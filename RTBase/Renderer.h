#pragma once

#include "Core.h"
#include "Sampling.h"
#include "Geometry.h"
#include "Imaging.h"
#include "Materials.h"
#include "Lights.h"
#include "Scene.h"
#include "GamesEngineeringBase.h"
#include <thread>
#include <functional>

class RayTracer
{
public:
	Scene* scene;
	GamesEngineeringBase::Window* canvas;
	Film* film;
	MTRandom *samplers;
	std::thread **threads;
	int numProcs;
	void init(Scene* _scene, GamesEngineeringBase::Window* _canvas)
	{
		scene = _scene;
		canvas = _canvas;
		film = new Film();
		film->init((unsigned int)scene->camera.width, (unsigned int)scene->camera.height, new BoxFilter());
		SYSTEM_INFO sysInfo;
		GetSystemInfo(&sysInfo);
		numProcs = sysInfo.dwNumberOfProcessors;
		threads = new std::thread*[numProcs];
		samplers = new MTRandom[numProcs];
		clear();
	}
	void clear()
	{
		film->clear();
	}
	//计算传入的空间中一点在随机一入射方向上的渲染结果，但是只计算能射到光源的入射方向，
	//也就是说只计算已有光源，且能直接照射到此点的一条光线
	//不包括反射效果，只计算了一个方向，取样时也是直接取已有光源，Next Event Estimation(nee)
	Colour computeDirect111(ShadingData shadingData, Sampler* sampler)
	{
		//直接跳过镜面类型
		if (shadingData.bsdf->isPureSpecular() == true)
		{
			return Colour(0.0f, 0.0f, 0.0f);
		}
		// Sample a light
		float pmf;
		Light* light = scene->sampleLight(sampler, pmf);
		// Sample a point on the light，返回光源上一点方向
		float pdf;
		Colour emitted;
		Vec3 p = light->sample(shadingData, sampler, emitted, pdf);
		//判断类型
		if (light->isArea())
		{
			// Calculate GTerm
			//wi就是入射方向
			Vec3 wi = p - shadingData.x;
			float l = wi.lengthSq();
			wi = wi.normalize();
			//计算几何项gterm，用于后续光照计算
			float GTerm = (max(Dot(wi, shadingData.sNormal), 0.0f) * max(-Dot(wi, light->normal(shadingData, wi)), 0.0f)) / l;
			if (GTerm > 0)
			{
				// Trace，可见性判断，是否被阻挡
				if (scene->visible(shadingData.x, p))
				{
					// Shade，计算这一方向的渲染效果
					return shadingData.bsdf->evaluate(shadingData, wi) * emitted * GTerm / (pmf * pdf);
				}
			}
		}
		else
	    //说明此光照是方向光，等没有固定源头,需要进行近似计算处理
		{
			// Calculate GTerm
			Vec3 wi = p;
			float GTerm = max(Dot(wi, shadingData.sNormal), 0.0f);
			if (GTerm > 0)
			{
				// Trace
				if (scene->visible(shadingData.x, shadingData.x + (p * 10000.0f)))
				{
					// Shade
					return shadingData.bsdf->evaluate(shadingData, wi) * emitted * GTerm / (pmf * pdf);
				}
			}
		}
		return Colour(0.0f, 0.0f, 0.0f);
	}

	Colour computeDirect(ShadingData shadingData, Sampler* sampler)
	{
		//直接跳过镜面类型
		if (shadingData.bsdf->isPureSpecular() == true)
		{
			return Colour(0.0f, 0.0f, 0.0f);
		}
		//用来存储最终光源
		Colour L(0.0f, 0.0f, 0.0f);
		//光源采样
		{
			// Sample a light
			float pmf;
			Light* light = scene->sampleLight(sampler, pmf);
			// Sample a point on the light，返回光源上一点方向
			float pdf;
			Colour emitted;
			Vec3 p = light->sample(shadingData, sampler, emitted, pdf);
			//判断类型
			if (light->isArea())
			{
				// Calculate GTerm
				//wi就是入射方向
				Vec3 wi = p - shadingData.x;
				float l = wi.lengthSq();
				wi = wi.normalize();
				//计算几何项gterm，用于后续光照计算
				float NoL = max(Dot(wi, shadingData.sNormal), 0.0f);
				float NoL_light = max(-Dot(wi, light->normal(shadingData, wi)), 0.0f);
				float GTerm = (NoL * NoL_light) / l;
				if (GTerm > 0)
				{
					// Trace，可见性判断，是否被阻挡
					if (scene->visible(shadingData.x, p))
					{
						Colour bsdfVal = shadingData.bsdf->evaluate(shadingData, wi);
						if (!bsdfVal.isBlack())
						{
							// ----------------------
							// PDF 转换（面积 → 方向）
							// ----------------------
							// p_light = pmf * pdfLightArea * (distance^2 / |n·wi|)
							float pdfLight = (pmf * pdf*l) / NoL_light;

							// BSDF pdf for same direction
							float pdfBSDF = shadingData.bsdf->PDF(shadingData, wi);

							// MIS 权重：balance heuristic
							float weight = pdfLight / (pdfLight + pdfBSDF + 1e-6f);

							// 最终值：emit * bsdf * G / pdfLight * weight
							//L = L + bsdfVal* emitted* GTerm* weight / pdfLight;
							Colour c = bsdfVal * emitted;
							c =c* GTerm;
							c =c* weight;
							c =c/ pdfLight;
							L =L+ c;
						}
					}
				}
			}
			else
				//说明此光照是方向光，等没有固定源头,需要进行近似计算处理
			{
				// Calculate GTerm
				Vec3 wi = p;
				float NoL = max(Dot(wi, shadingData.sNormal), 0.0f);
				if (NoL > 0)
				{
					// Trace
					if (scene->visible(shadingData.x, shadingData.x + (p * 10000.0f)))
					{
						Colour bsdfVal = shadingData.bsdf->evaluate(shadingData, wi);
						if (!bsdfVal.isBlack())
						{
							// 假设 pdfLight 为方向光的方向PDF
							float pdfLight = pmf * pdf;
							float pdfBSDF = shadingData.bsdf->PDF(shadingData, wi);
							float weight = pdfLight / (pdfLight + pdfBSDF + 1e-6f);

							//L = L+bsdfVal * emitted * NoL * weight / pdfLight;
							Colour c = bsdfVal * emitted;
							c =c* NoL;
							c =c* weight;
							c =c/ pdfLight;
							L = L+c;
						}
					}
				}
			}
		}
		//BSDF采样
		{
			float pdfBSDF;
			Colour bsdfVal;
			Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
			pdfBSDF = SamplingDistributions::cosineHemispherePDF(wi);
			wi = shadingData.frame.toWorld(wi);
			bsdfVal = shadingData.bsdf->evaluate(shadingData, wi);

			if (pdfBSDF > 1e-6f && !bsdfVal.isBlack())
			{
				Ray shadowRay(shadingData.x + wi * EPSILON, wi);
				IntersectionData hit = scene->traverse(shadowRay);
				if (hit.t < FLT_MAX)
				{
					ShadingData hitShading = scene->calculateShadingData(hit, shadowRay);
					if (hitShading.bsdf && hitShading.bsdf->isLight())
					{
						Colour emitted = hitShading.bsdf->emit(hitShading, -wi);
						if (!emitted.isBlack())
						{
							float NoL = max(Dot(shadingData.sNormal, wi), 0.0f);
							float pdfLight = scene->lightPdf(shadingData, wi);  // 你可能需要实现这个
							float weight = pdfBSDF / (pdfBSDF + pdfLight + 1e-6f);

							//L = L+bsdfVal * emitted * NoL * weight / pdfBSDF;
							Colour c = bsdfVal * emitted;
							c =c* NoL;
							c =c* weight;
							c =c/ pdfBSDF;
							L =L+ c;
						}
					}
				}
			}
		}
	}


	//canhitlight用于判断能否命中光源
	//这个是直接光加上间接光
	Colour pathTrace(Ray& r, Colour& pathThroughput, int depth, Sampler* sampler, bool canHitLight = true)
	{
		//摄像机发出的光线r遍历场景，得到相交点
		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		//如果确实命中
		if (shadingData.t < FLT_MAX)
		{
			//如果是光源
			if (shadingData.bsdf->isLight())
			{
				if (canHitLight == true)
				{
					//返回光的颜色，这里的path这个值代表光的反射衰减率，保留原有效果的多少
					return pathThroughput * shadingData.bsdf->emit(shadingData, shadingData.wo);
				}
				else
				{
					return Colour(0.0f, 0.0f, 0.0f);
				}
			}
			//不是光，说明是物体，返回一个方向的颜色，一样乘衰减率
			Colour direct = pathThroughput * computeDirect(shadingData, sampler);
			//判断光的弹射次数，达标就返回
			if (depth > 5)
			{
				return direct;
			}
			//若是光的能量过低，也停止计算，返回
			float russianRouletteProbability = min(pathThroughput.Lum(), 0.9f);
			if (sampler->next() < russianRouletteProbability)
			{
				//并且更新衰减率
				pathThroughput = pathThroughput / russianRouletteProbability;
			}
			else
			{
				return direct;
			}
			Colour bsdf;
			float pdf;
			//bsdf采样，不同于光源采样，随机取样
			Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
			pdf = SamplingDistributions::cosineHemispherePDF(wi);
			wi = shadingData.frame.toWorld(wi);
			bsdf = shadingData.bsdf->evaluate(shadingData, wi);
			//这里递归调用了这个方法，以这个交点当相机再次射线，模拟一次光线反射
			pathThroughput = pathThroughput * bsdf * fabsf(Dot(wi, shadingData.sNormal)) / pdf;
			r.init(shadingData.x + (wi * EPSILON), wi);
			return (direct + pathTrace(r, pathThroughput, depth + 1, sampler, shadingData.bsdf->isPureSpecular()));
		}
		//没打到，返回背景材质
		return scene->background->evaluate(shadingData, r.dir);
	}


	//从相机发出的一条射线返回一个颜色
	//这个方法使用后就没有间接光了
	Colour direct(Ray& r, Sampler* sampler)
	{
		//摄像机发出的光线r遍历场景，得到相交点
		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		//判断确实打到物体
		if (shadingData.t < FLT_MAX)
		{
			//判断打到的物体是否是光源
			if (shadingData.bsdf->isLight())
			{
				//是光源就返回光的颜色
				return shadingData.bsdf->emit(shadingData, shadingData.wo);
			}
			//不是光则计算一次compute，得到此点朝向某方向的颜色
			return computeDirect(shadingData, sampler);
		}
		//未打到物体则返回黑色
		return Colour(0.0f, 0.0f, 0.0f);
	}
	
	Colour albedo(Ray& r)
	{
		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		if (shadingData.t < FLT_MAX)
		{
			if (shadingData.bsdf->isLight())
			{
				return shadingData.bsdf->emit(shadingData, shadingData.wo);
			}
			return shadingData.bsdf->evaluate(shadingData, Vec3(0, 1, 0));
		}
		return scene->background->evaluate(shadingData, r.dir);
	}
	Colour viewNormals(Ray& r)
	{
		IntersectionData intersection = scene->traverse(r);
		if (intersection.t < FLT_MAX)
		{
			ShadingData shadingData = scene->calculateShadingData(intersection, r);
			return Colour(fabsf(shadingData.sNormal.x), fabsf(shadingData.sNormal.y), fabsf(shadingData.sNormal.z));
		}
		return Colour(0.0f, 0.0f, 0.0f);
	}

	 //单线程版
	//void render()
	//{
	//	film->incrementSPP();
	//	for (unsigned int y = 0; y < film->height; y++)
	//	{
	//		for (unsigned int x = 0; x < film->width; x++)
	//		{//这里就是光线追踪（遍历每一个像素点，生成光线，返回一个颜色绘制）
	//		    //这里取每个像素点的中心来生成光线
	//			//float px = x + 0.5f;
	//			//float py = y + 0.5f;
	//			float px = x + samplers->next();
	//			float py = y + samplers->next();
	//			Ray ray = scene->camera.generateRay(px, py);
	//			//Colour col = viewNormals(ray);
	//			//Colour col = albedo(ray);
	//			Colour temp(1.f, 1.f, 1.f);
	//			Colour col = pathTrace(ray, temp, 0, samplers);
	//			film->splat(px, py, col);
	//			//把法线的颜色变成255版本
	//			unsigned char r = (unsigned char)(col.r * 255);
	//			unsigned char g = (unsigned char)(col.g * 255);
	//			unsigned char b = (unsigned char)(col.b * 255);
	//			film->tonemap(x, y, r, g, b);
	//			canvas->draw(x, y, r, g, b);
	//		}
	//	}
	//}

		//多线程
	void render()
	{
		static const int TILE_SIZE = 32;
		film->incrementSPP();
		std::vector<std::thread> workers;
		int numThreads = numProcs;
		int numTilesX = (film->width + TILE_SIZE - 1) / TILE_SIZE;
		int numTilesY = (film->height + TILE_SIZE - 1) / TILE_SIZE;


		auto renderTile = [&](int tileX, int tileY, int threadId)
			{
				int startX = tileX * TILE_SIZE;
				int startY = tileY * TILE_SIZE;
				int endX = min(startX + TILE_SIZE, (int)film->width);
				int endY = min(startY + TILE_SIZE, (int)film->height);

				for (int y = startY; y < endY; y++)
				{
					for (int x = startX; x < endX; x++)
					{
						//float px = x + 0.5f;
						//float py = y + 0.5f;

						float px = x + samplers->next();
						float py = y + samplers->next();
						Ray ray = scene->camera.generateRay(px, py);

						//Colour col = viewNormals(ray);
			            //Colour col = albedo(ray);

						Colour pathThroughput(1.0f, 1.0f, 1.0f);
						//Colour col = direct(ray, &samplers[threadId]);
						Colour col = pathTrace(ray, pathThroughput, 0, &samplers[threadId]);

						film->splat(px, py, col);

						unsigned char r = (unsigned char)(col.r * 255);
						unsigned char g = (unsigned char)(col.g * 255);
						unsigned char b = (unsigned char)(col.b * 255);
						film->tonemap(x, y, r, g, b);
						canvas->draw(x, y, r, g, b);
					}
				}
			};

		auto workerFunc = [&](int threadId)
			{
				for (int tileY = 0; tileY < numTilesY; tileY++)
				{
					for (int tileX = 0; tileX < numTilesX; tileX++)
					{
						if (((tileY * numTilesX) + tileX) % numThreads == threadId)
						{
							renderTile(tileX, tileY, threadId);
						}
					}
				}
			};

		for (int i = 0; i < numThreads; i++)
		{
			workers.emplace_back(workerFunc, i);
		}

		for (auto& worker : workers)
		{
			worker.join();
		}

	}

	int getSPP()
	{
		return film->SPP;
	}
	void saveHDR(std::string filename)
	{
		film->save(filename);
	}
	void savePNG(std::string filename)
	{
		stbi_write_png(filename.c_str(), canvas->getWidth(), canvas->getHeight(), 3, canvas->getBackBuffer(), canvas->getWidth() * 3);
	}
};