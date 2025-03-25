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
	Colour computeDirect1(ShadingData shadingData, Sampler* sampler)
	{
		if (shadingData.bsdf->isPureSpecular() == true)
		{
			return Colour(0.0f, 0.0f, 0.0f);
		}
		// Sample a light
		float pmf;
		Light* light = scene->sampleLight(sampler, pmf);
		// Sample a point on the light
		float pdf;
		Colour emitted;
		Vec3 p = light->sample(shadingData, sampler, emitted, pdf);
		if (light->isArea())
		{
			// Calculate GTerm
			Vec3 wi = p - shadingData.x;
			float l = wi.lengthSq();
			wi = wi.normalize();
			float GTerm = (max(Dot(wi, shadingData.sNormal), 0.0f) * max(-Dot(wi, light->normal(shadingData, wi)), 0.0f)) / l;
			if (GTerm > 0)
			{
				// Trace
				if (scene->visible(shadingData.x, p))
				{
					// Shade
					return shadingData.bsdf->evaluate(shadingData, wi) * emitted * GTerm / (pmf * pdf);
				}
			}
		}
		else
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

	Colour computeDirect11(ShadingData shadingData, Sampler* sampler)
	{
		if (shadingData.bsdf->isPureSpecular()) return Colour(0.0f, 0.0f, 0.0f);

		// Sample a light
		float pmf;
		Light* light = scene->sampleLight(sampler, pmf);
		if (!light || pmf < 1e-6f) return Colour(0.0f, 0.0f, 0.0f);

		// Sample a point/direction on the light
		float pdf;
		Colour emitted;
		Vec3 sample = light->sample(shadingData, sampler, emitted, pdf);

		if (pdf < 1e-6f || !std::isfinite(pdf)) return Colour(0.0f, 0.0f, 0.0f);

		Vec3 wi;
		float G = 1.0f;

		if (light->isArea())
		{
			Vec3 p = sample;
			wi = p - shadingData.x;
			float dist2 = wi.lengthSq();
			wi = wi.normalize();

			float cosSurface = max(Dot(shadingData.sNormal, wi), 0.0f);
			float cosLight = max(-Dot(wi, light->normal(shadingData, wi)), 1e-4f); // avoid div0

			if (cosSurface <= 0.0f || cosLight <= 0.0f) return Colour(0.0f, 0.0f, 0.0f);
			if (!scene->visible(shadingData.x, p)) return Colour(0.0f, 0.0f, 0.0f);

			G = (cosSurface * cosLight) / dist2;

			// Convert area pdf to solid angle domain
			pdf = (pdf * dist2) / cosLight;

			if (pdf < 1e-6f || !std::isfinite(pdf)) return Colour(0.0f, 0.0f, 0.0f);
		}
		else
		{
			// Environment / infinite light
			wi = sample;
			float cosSurface = max(Dot(shadingData.sNormal, wi), 0.0f);
			G = cosSurface;

			if (cosSurface <= 0.0f) return Colour(0.0f, 0.0f, 0.0f);
			if (!scene->visible(shadingData.x, shadingData.x + wi * 1e4f)) return Colour(0.0f, 0.0f, 0.0f);
		}

		Colour bsdfVal = shadingData.bsdf->evaluate(shadingData, wi);

		if (!std::isfinite(bsdfVal.r) || !std::isfinite(emitted.r)) return Colour(0.0f, 0.0f, 0.0f);

		Colour result = bsdfVal * emitted * G / (pdf * pmf);

		// 额外保险：强制 clamp 太亮的异常值（可选）
		result.r = min(result.r, 50.0f);
		result.g = min(result.g, 50.0f);
		result.b = min(result.b, 50.0f);

		return result;
	}

	Colour computeDirect(ShadingData shadingData, Sampler* sampler)
	{
		if (shadingData.bsdf->isPureSpecular())
			return Colour(0.0f, 0.0f, 0.0f);

		// 1. Sample a light
		float pmf;
		Light* light = scene->sampleLight(sampler, pmf);
		if (!light || pmf < 1e-6f)
			return Colour(0.0f, 0.0f, 0.0f);

		// 2. Sample a point or direction on the light
		float pdf;
		Colour emitted;
		Vec3 sample = light->sample(shadingData, sampler, emitted, pdf);

		if (pdf < 1e-6f || !std::isfinite(pdf))
			return Colour(0.0f, 0.0f, 0.0f);

		Vec3 wi;
		float G = 1.0f;

		if (light->isArea())
		{
			Vec3 p = sample;
			wi = p - shadingData.x;
			float dist2 = wi.lengthSq();
			float dist = sqrtf(dist2);
			wi = wi.normalize();

			float cosSurface = max(Dot(shadingData.sNormal, wi), 0.0f);
			float cosLight = max(-Dot(wi, light->normal(shadingData, wi)), 1e-2f); // ⚠️ 更高的下限，抗白边

			if (cosSurface <= 0.0f || cosLight <= 0.0f)
				return Colour(0.0f, 0.0f, 0.0f);

			G = (cosSurface * cosLight) / dist2;

			// Convert area PDF to solid angle domain
			pdf = (pdf * dist2) / cosLight;
			if (pdf < 1e-6f || !std::isfinite(pdf))
				return Colour(0.0f, 0.0f, 0.0f);

			// ✅ 加 offset，防止自交导致 shadow ray 撞到自己
			Vec3 shadowOrigin = shadingData.x + wi * 1e-4f;
			if (!scene->visible(shadowOrigin, p))
				return Colour(0.0f, 0.0f, 0.0f);
		}
		else
		{
			// Environment light
			wi = sample;
			float cosSurface = max(Dot(shadingData.sNormal, wi), 0.0f);
			G = cosSurface;

			if (cosSurface <= 0.0f)
				return Colour(0.0f, 0.0f, 0.0f);

			Vec3 shadowOrigin = shadingData.x + wi * 1e-4f;
			if (!scene->visible(shadowOrigin, shadowOrigin + wi * 1e4f))
				return Colour(0.0f, 0.0f, 0.0f);
		}

		// 3. BSDF eval
		Colour bsdfVal = shadingData.bsdf->evaluate(shadingData, wi);
		if (!std::isfinite(bsdfVal.r) || !std::isfinite(emitted.r))
			return Colour(0.0f, 0.0f, 0.0f);

		// 4. Final result
		Colour result = bsdfVal * emitted * G / (pdf * pmf);

		// ✅ Soft clamp to avoid overbright edges
		float maxIntensity = 20.0f;
		result.r = min(result.r, maxIntensity);
		result.g = min(result.g, maxIntensity);
		result.b = min(result.b, maxIntensity);

		return result;
	}

	Colour pathTrace(Ray& r, Colour& pathThroughput, int depth, Sampler* sampler, bool canHitLight = true)
	{
		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		if (shadingData.t < FLT_MAX)
		{
			if (shadingData.bsdf->isLight())
			{
				if (canHitLight == true)
				{
					return pathThroughput * shadingData.bsdf->emit(shadingData, shadingData.wo);
				}
				else
				{
					return Colour(0.0f, 0.0f, 0.0f);
				}
			}
			Colour direct = pathThroughput * computeDirect(shadingData, sampler);
			if (depth > 3)
			{
				return direct;
			}
			float russianRouletteProbability = min(pathThroughput.Lum(), 0.9f);
			if (sampler->next() < russianRouletteProbability)
			{
				pathThroughput = pathThroughput / russianRouletteProbability;
			}
			else
			{
				return direct;
			}
			Colour bsdf;
			float pdf;
			Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
			pdf = SamplingDistributions::cosineHemispherePDF(wi);
			wi = shadingData.frame.toWorld(wi);
			bsdf = shadingData.bsdf->evaluate(shadingData, wi);
			pathThroughput = pathThroughput * bsdf * fabsf(Dot(wi, shadingData.sNormal)) / pdf;
			r.init(shadingData.x + (wi * EPSILON), wi);
			return (direct + pathTrace(r, pathThroughput, depth + 1, sampler, shadingData.bsdf->isPureSpecular()));
		}
		return scene->background->evaluate(shadingData, r.dir);
	}
	Colour direct(Ray& r, Sampler* sampler)
	{
		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		if (shadingData.t < FLT_MAX)
		{
			if (shadingData.bsdf->isLight())
			{
				return shadingData.bsdf->emit(shadingData, shadingData.wo);
			}
			return computeDirect(shadingData, sampler);
		}
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
	//			float px = x + 0.5f;
	//			float py = y + 0.5f;
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
						float px = x + 0.5f;
						float py = y + 0.5f;

						Ray ray = scene->camera.generateRay(px, py);

						//Colour col = viewNormals(ray);
			            //Colour col = albedo(ray);

						Colour pathThroughput(1.0f, 1.0f, 1.0f);
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