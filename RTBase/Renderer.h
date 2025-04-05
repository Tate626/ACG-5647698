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

//for adaptive sampler
struct Tile {
	//position
	int tileX, tileY;
	float variance = 0.0f;
	//weight
	float weight = 1.0f;
	//the number of sample
	int currentSPP = 0;
	int targetSPP = 1;
};

struct VPL {
	ShadingData shadingData;
	Colour Le;
	bool fromLight = false;
	bool isEnvironment = false; //from environment light or not
	Vec3 envDir;               //if from environmeny light, memory the direction
};

class RayTracer
{
public:
	Scene* scene;
	GamesEngineeringBase::Window* canvas;
	Film* film;
	MTRandom *samplers;
	std::thread **threads;
	int numProcs;
	std::vector<VPL> vpls;
	void init(Scene* _scene, GamesEngineeringBase::Window* _canvas)
	{
		scene = _scene;
		canvas = _canvas;
		film = new Film();
		//改成高斯滤波器了
		//film->init((unsigned int)scene->camera.width, (unsigned int)scene->camera.height, new GaussianFilter());
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

	//Used for path tracing: computeDirect  pathTrace  direct.
	//计算传入的空间中一点在随机一入射方向上的渲染结果，但是只计算能射到光源的入射方向，
	//也就是说只计算已有光源，且能直接照射到此点的一条光线
	//不包括反射效果，只计算了一个方向，取样时也是直接取已有光源，Next Event Estimation(nee)
	Colour computeDirect(ShadingData shadingData, Sampler* sampler)
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

		Vec3 wi;
		float GTerm = 0.0f;
		//判断类型
		if (light->isArea())
		{
			// Calculate GTerm
			//wi就是入射方向
			wi = p - shadingData.x;
			float l = wi.lengthSq();
			wi = wi.normalize();
			//计算几何项gterm，用于后续光照计算
			GTerm = (max(Dot(wi, shadingData.sNormal), 0.0f) * max(-Dot(wi, light->normal(wi)), 0.0f)) / l;
			if (GTerm > 0)
			{
				// Trace，可见性判断，是否被阻挡
				if (scene->visible(shadingData.x, p))
				{
					// 计算 BSDF 对这个方向的 pdf
					float bsdf_pdf = shadingData.bsdf->PDF(shadingData, wi);
					// 光源采样的 pdf：结合采样光源和采样点的概率
					float light_pdf = pmf * pdf;
					// 计算 MIS 权重（加个小量防止除0）
					float misWeight = light_pdf / (light_pdf + bsdf_pdf + 1e-6f);
					// Shade，计算这一方向的渲染效果
					//return shadingData.bsdf->evaluate(shadingData, wi) * emitted * GTerm / (pmf * pdf);
					return shadingData.bsdf->evaluate(shadingData, wi) * emitted * GTerm * misWeight / light_pdf;
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
					float bsdf_pdf = shadingData.bsdf->PDF(shadingData, wi);
					float light_pdf = pmf * pdf;
					float misWeight = light_pdf / (light_pdf + bsdf_pdf + 1e-6f);
					return shadingData.bsdf->evaluate(shadingData, wi) * emitted * GTerm * misWeight / light_pdf;
					// Shade
					//return shadingData.bsdf->evaluate(shadingData, wi) * emitted * GTerm / (pmf * pdf);
				}
			}
		}
		return Colour(0.0f, 0.0f, 0.0f);
	}
	//canhitlight用于判断能否命中光源
	//这个是直接光加上间接光
	Colour pathTrace(Ray& r, Colour& pathThroughput, int depth, Sampler* sampler,Colour& a, Vec3& N, bool canHitLight = true)
	{
		//摄像机发出的光线r遍历场景，得到相交点
		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		N = shadingData.sNormal.normalize();

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
			if (depth >8)
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
			//Colour bsdf;
			//float pdf;
			//Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
			//pdf = SamplingDistributions::cosineHemispherePDF(wi);
			//wi = shadingData.frame.toWorld(wi);
			//bsdf = shadingData.bsdf->evaluate(shadingData, wi);
			Colour indirect;
			float pdf;
			Vec3 wi = shadingData.bsdf->sample(shadingData, sampler, indirect, pdf,a);
			//这里递归调用了这个方法，以这个交点当相机再次射线，模拟一次光线反射
			pathThroughput = pathThroughput * indirect * fabsf(Dot(wi, shadingData.sNormal)) / pdf;
			r.init(shadingData.x + (wi * EPSILON), wi);
			// 递归调用 pathTrace 来计算间接光照（反射/折射后采样的贡献）
			return (direct + pathTrace(r, pathThroughput, depth + 1, sampler, a, N,shadingData.bsdf->isPureSpecular()));
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


	//Used for Instant Radiosity: VPLTracePath  traceVPLs  computeForPixel.
	const size_t MAX_VPLS = 100;//Limit the max number of VPLs
	//This method simulates the path of a single light ray and generates a VPL at each intersection point.
	void VPLTracePath(Ray& r, Colour pathThroughput, Colour Le, Sampler* sampler, int depth) {
		//stop if the number reach the max
		if (vpls.size() >= MAX_VPLS) return;

		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		if (shadingData.t < FLT_MAX) {
			if (depth > 5) {
				return;
			}
			// if it is a mirror,return.I use path tracing to render mirror. it can be seen in computeForPixel()
			if (shadingData.bsdf->isPureSpecular() == true) {
				return;
			}
			//create a new VPL
			if (vpls.size() < MAX_VPLS) {
				VPL vpl;
				vpl.shadingData = shadingData;
				vpl.Le = pathThroughput * Le;
				vpl.fromLight = false;
				vpl.isEnvironment = false; 
				vpls.push_back(vpl);
			}
			else {
				return;
			}
			//Recursively call to simulate light transport.
			float pdfBSDF = 0.0f;
			Colour bsdfVal;
			Colour A;
			Vec3 newDir = shadingData.bsdf->sample(shadingData, sampler, bsdfVal, pdfBSDF, A);
			pathThroughput = pathThroughput * bsdfVal
				* fabsf(Dot(newDir, shadingData.sNormal)) / pdfBSDF;
			float russianRouletteProbability = min(pathThroughput.Lum(), 0.9f);
			if (sampler->next() < russianRouletteProbability) {
				pathThroughput = pathThroughput / russianRouletteProbability;
			}
			else {
				return;
			}
			r.init(shadingData.x + (newDir * EPSILON), newDir);
			VPLTracePath(r, pathThroughput, Le, sampler, depth + 1);
		}
	}
	//This method randomly samples and generates N_VPLs initial light source VPLs. 
	//If the light is not an environment light, calls VPLTracePath to generate all VPLs along the path.
	void traceVPLs(Sampler* sampler, int N_VPLs) {
		vpls.clear();
		vpls.reserve(MAX_VPLS);

		//generate N_VPLs initial VPLs
		for (int i = 0; i < N_VPLs; ++i) {
			float pmf;
			Light* light = scene->sampleLight(sampler, pmf);
			if (pmf <= 0.0f) {
				continue;
			}
			//distinguish environment light and other lights
			if (light->isENV()) {
				//sample a direction
				float pdfDir;
				Vec3 envDir = light->sampleDirectionFromLight(sampler, pdfDir);
				if (pdfDir <= 0.0f) continue;

				//build the VPL
				Colour envEmit = light->getEmitted(envDir);
				Colour LE = envEmit / (pmf * pdfDir * float(N_VPLs));
				VPL vpl;
				// environment light does not require a real position or normal.
				vpl.shadingData.x = Vec3(0.0f, 0.0f, 0.0f);
				vpl.shadingData.sNormal = Vec3(0.0f, 0.0f, 0.0f);
				vpl.Le = LE;
				vpl.fromLight = true;
				vpl.isEnvironment = true;
				vpl.envDir = envDir; 
				vpls.push_back(vpl);
			}
			//for other lights
			else {
				//sample a position
				float pdfPosition;
				Vec3 p = light->samplePositionFromLight(sampler, pdfPosition);
				if (pdfPosition <= 0.0f) continue;

				ShadingData shadingData(p, light->normal(p));
				Colour emitted = light->getEmitted(p);
				Colour LE = emitted / (pmf * pdfPosition * float(N_VPLs));

				//the first VPL
				VPL vpl;
				vpl.shadingData = shadingData;
				vpl.Le = LE;
				vpl.fromLight = true;
				vpl.isEnvironment = false;
				vpls.push_back(vpl);

				//generate ray from the first VPL and use VPLTracePath to get more VPLs.
				float pdfDir;
				Vec3 wi = light->sampleDirectionFromLight(sampler, pdfDir);
				if (pdfDir <= 0.0f) {
					continue;
				}
				float cosLight = Dot(wi, light->normal(p));
				if (cosLight < 0.0f) cosLight = 0.0f;
				Colour le = emitted * cosLight / (pmf * pdfPosition * pdfDir * float(N_VPLs));
				Ray r;
				r.init(shadingData.x + (wi * EPSILON), wi);
				Colour pathThroughput(1.f, 1.f, 1.f);
				VPLTracePath(r, pathThroughput, le, sampler, 0);
			}
		}
	}
	//For a pixel on the screen, compute the intersection to get its position in the world. 
	//Then,traverse all VPLs and compute their color contribution to this position.
	//The sum of these contributions is the final color seen by the camera.
	Colour computeForPixel(Ray& r, Sampler* sampler, Colour& a, Vec3& N) {
		Colour result(0.0f, 0.0f, 0.0f);
		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		N = shadingData.sNormal.normalize();

		//if hitting a light return its colour.
		if(shadingData.bsdf!=NULL)
		if (shadingData.bsdf->isLight()) {
			return shadingData.bsdf->emission;
		}

		if (shadingData.t < FLT_MAX) {
			//if it is a mirror, using path tracing
			if (shadingData.bsdf->isPureSpecular() == true) {
				Colour throughput(1.0f,1.0f,1.0f);
				return pathTrace(r, throughput, 0, sampler, a, N);
			}

			// traverse all VPLs to get the colour of this pixel
			for (const VPL& vpl : vpls) {
				if (vpl.isEnvironment) {
					//path 1: environment light
					Vec3 wi = vpl.envDir;
					float cos = max(0.0f, Dot(wi, shadingData.sNormal));
					if (cos < 1e-4f) {
						continue;
					}
					Ray testRay(shadingData.x + wi * EPSILON, wi);
					IntersectionData testIsect = scene->traverse(testRay);
					bool visibleEnv = (testIsect.t >= FLT_MAX);
					if (visibleEnv) {
						Colour bsdfVal = shadingData.bsdf->evaluate(shadingData, wi);
						Colour contrib = vpl.Le * bsdfVal * cos;
						result = result + contrib;
					}
				}
				else {
					//path 1: other lights
					Vec3 vplPos = vpl.shadingData.x;
					Vec3 toVPL = vplPos - shadingData.x;
					float dist2 = Dot(toVPL, toVPL);
					if (dist2 < 1e-8f) {
						continue;
					}
					Vec3 wi = toVPL.normalize();

					float cosh = max(0.0f, Dot(wi, shadingData.sNormal));
					float cost = max(0.0f, -Dot(wi, vpl.shadingData.sNormal));
					if (cosh < 0.1f || cost < 0.1f) {
						continue;
					}
					float G = (cosh * cost) / (dist2 + EPSILON);
					if (G < 1e-4f) {
						continue;
					}
					if (!scene->visible(shadingData.x, vplPos)) {
						continue;
					}
					Colour temp(0.0f, 0.0f, 0.0f);
					if (vpl.fromLight) {
						// if the VPL is a light
						temp = vpl.Le
							* shadingData.bsdf->evaluate(shadingData, wi)
							* G;
					}
					else {
						// if the VPL is a bsdf
						temp = vpl.Le
							* vpl.shadingData.bsdf->evaluate(vpl.shadingData, -wi)
							* shadingData.bsdf->evaluate(shadingData, wi)
							* G;
					}
					result = result + temp;
				}
			}
			return result;
		}
		return scene->background->evaluate(shadingData, r.dir);
	}


	//Used for light tracing: connectToCamera  lightTracePath  lightTrace.
	//用于检测输入点p是否在摄像机视线范围内
	void connectToCamera(const Vec3& p, const Vec3& n, const Colour& col) {
		float x, y;
		//将 3D 点 p 投影到图像平面上
		if (scene->camera.projectOntoCamera(p, x, y)) {
			
			// 计算最终颜色
			Vec3 toCamera = (scene->camera.origin - p).normalize();
			float cosTheta = Dot(scene->camera.viewDirection, toCamera);

			//cosTheta = max(cosTheta, 0.01f);
			//if (cosTheta <= 0.0f) return; // 背对相机，不可见
			//float distanceSquared = (scene->camera.origin - p).lengthSq();
			//float We = cosTheta / (scene->camera.Afilm * distanceSquared * pow(cosTheta, 3));
			float We = 1.0f / (scene->camera.Afilm * pow(cosTheta, 4));
			Colour finalContribution = col*We;
			
			//调用splat累计贡献
			film->splat(static_cast<int>(x), static_cast<int>(y), finalContribution);
		}
	}
    //传入光线遍历场景，相交调用connect，再递归反射
	void lightTracePath(Ray& r, Colour pathThroughput, int depth, Colour Le, Sampler* sampler, Colour& a, Vec3& N) {
		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		N = shadingData.sNormal.normalize();
		//如果确实命中
		if (shadingData.t < FLT_MAX)
		{
			float cos = Dot(shadingData.sNormal, shadingData.wo);
			if (cos <= 0.0f) return;

	/*		Colour f_cam = shadingData.bsdf->evaluate(shadingData, shadingData.wo);
			Colour finalCol = pathThroughput * f_cam * Le * cos;*/
			Vec3 toCamera = (scene->camera.origin - shadingData.x).normalize();
			Colour f_cam = shadingData.bsdf->evaluate(shadingData, toCamera);
			Colour finalCol = pathThroughput * f_cam * Le;
			// 1. 连接相机
			connectToCamera(shadingData.x, shadingData.sNormal, finalCol);

			// 2. 采样 BSDF 新方向
			float pdf;
			Colour f;
			Vec3 wi = shadingData.bsdf->sample(shadingData, sampler, f, pdf, a);

			float cosTheta = Dot(wi, shadingData.sNormal);
			if (pdf <= 0.0f || cosTheta <= 0.0f) return;

			//判断光的弹射次数，达标就返回
			if (depth > 8)
			{
				return;
			}
			// 4. 更新 pathThroughput（PPT公式：BSDF × cosθ / pdf）
			pathThroughput = pathThroughput * f * cosTheta / pdf;

			// 5. 构建新 Ray 并递归（PPT的“Goto 1”）
			r.init(shadingData.x + (wi * EPSILON), wi);
			lightTracePath(r, pathThroughput, depth + 1, Le, sampler, a, N);
			return;
		}
		//没打到
		return;
	}
	// 入口函数：从光源发射光线开始追踪
	void lightTrace(Sampler* sampler, Colour& a, Vec3& N) {
		//采样一个光源（Sample a light source）
		float pmf;
		Light* light = scene->sampleLight(sampler,pmf);

		//从光源采样一个位置和一个出射方向（Sample position and direction）
		float pdfPos, pdfDir;
		Vec3 p, wi;
		ShadingData temp;
		Colour Le;

		p = light->samplePositionFromLight(sampler, pdfPos); // 采样光源表面一个点
		wi = light->sampleDirectionFromLight(sampler, pdfDir); // 从该点采样一个发射方向

		// [PPT 第3步] 计算从该方向出射的光强（Le），再除以位置采样概率密度
		// PPT公式：col = Le(x0 → x1) / pA(x0)
		Le = light->evaluate(temp,-wi) / pdfPos; // 注意方向是 -wi，因为 evaluate 期望入射方向

		//初始化路径能量为1（Path throughput）
		Colour pathThroughput(1.0f,1.0f,1.0f);

		//构造初始光线，从光源出发
		Ray r(p + wi * EPSILON, wi);

		//调用路径追踪函数，从光源起始点开始追踪光线
		lightTracePath(r, pathThroughput, 0, Le, sampler, a, N);
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

	//This render is a basic single-threaded version
	//Some interfaces have been modified, so it may no longer be used.
	 
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

	//This render is a basic multiple-threaded version
	//Some interfaces have been modified, so it may no longer be used.
	
	//void render()
	//{
	//	static const int TILE_SIZE = 32;
	//	film->incrementSPP();
	//	std::vector<std::thread> workers;
	//	int numThreads = numProcs;
	//	int numTilesX = (film->width + TILE_SIZE - 1) / TILE_SIZE;
	//	int numTilesY = (film->height + TILE_SIZE - 1) / TILE_SIZE;
	//	auto renderTile = [&](int tileX, int tileY, int threadId)
	//		{
	//			int startX = tileX * TILE_SIZE;
	//			int startY = tileY * TILE_SIZE;
	//			int endX = min(startX + TILE_SIZE, (int)film->width);
	//			int endY = min(startY + TILE_SIZE, (int)film->height);
	//			for (int y = startY; y < endY; y++)
	//			{
	//				for (int x = startX; x < endX; x++)
	//				{
	//					float px = x + samplers->next();
	//					float py = y + samplers->next();
	//					Ray ray = scene->camera.generateRay(px, py);
	//					//Colour col = viewNormals(ray);
	//		            //Colour col = albedo(ray);
	//					Colour pathThroughput(1.0f, 1.0f, 1.0f);
	//					//Colour col = direct(ray, &samplers[threadId]);
	//					Colour col = pathTrace(ray, pathThroughput, 0, &samplers[threadId]);
	//					film->splat(px, py, col);
	//					unsigned char r = (unsigned char)(col.r * 255);
	//					unsigned char g = (unsigned char)(col.g * 255);
	//					unsigned char b = (unsigned char)(col.b * 255);
	//					film->tonemap(x, y, r, g, b);
	//					canvas->draw(x, y, r, g, b);
	//				}
	//			}
	//		};
	//	auto workerFunc = [&](int threadId)
	//		{
	//			for (int tileY = 0; tileY < numTilesY; tileY++)
	//			{
	//				for (int tileX = 0; tileX < numTilesX; tileX++)
	//				{
	//					if (((tileY * numTilesX) + tileX) % numThreads == threadId)
	//					{
	//						renderTile(tileX, tileY, threadId);
	//					}
	//				}
	//			}
	//		};
	//	for (int i = 0; i < numThreads; i++)
	//	{
	//		workers.emplace_back(workerFunc, i);
	//	}
	//	for (auto& worker : workers)
	//	{
	//		worker.join();
	//	}
	//}

	//This render is a multiple-threaded version with adaptive sampler and denoise
	//Using path tracing
	//The denoiser is automatically applied after the 10th accumulated frame of the same image.

	void render()
	{
		static const int TILE_SIZE = 32;
		//memory the number of render
		film->incrementSPP();
		std::vector<std::thread> workers;
		std::vector<Tile> tiles;
		int numThreads = numProcs;
		int numTilesX = (film->width + TILE_SIZE - 1) / TILE_SIZE;
		int numTilesY = (film->height + TILE_SIZE - 1) / TILE_SIZE;
		//initialize all tiles
		for (int tileY = 0; tileY < numTilesY; ++tileY)
		{
			for (int tileX = 0; tileX < numTilesX; ++tileX)
			{
				Tile tile;
				tile.tileX = tileX;
				tile.tileY = tileY;
				tiles.push_back(tile);
			}
		}
		//compute one tile
		auto renderTile = [&](Tile& tile, int threadId)
			{
				int startX = tile.tileX * TILE_SIZE;
				int startY = tile.tileY * TILE_SIZE;
				int endX = min(startX + TILE_SIZE, (int)film->width);
				int endY = min(startY + TILE_SIZE, (int)film->height);
				for (int s = 0; s < tile.targetSPP; ++s)
				{
					for (int y = startY; y < endY; y++)
					{
						for (int x = startX; x < endX; x++)
						{
							float px = x + samplers->next();
							float py = y + samplers->next();
							Ray ray = scene->camera.generateRay(px, py);
							Colour pathThroughput(1.0f, 1.0f, 1.0f);
							Colour A;
							Vec3 N;
							Colour col = pathTrace(ray, pathThroughput, 0, &samplers[threadId],A,N);
							film->splat(px, py, col);
							//use for denoising
							film->AOV(int(px),int(py),A,N);
						}
					}
				}
				tile.currentSPP += tile.targetSPP;
			};
		//traverse all tiles
		auto workerFunc = [&](int threadId)
			{
				for (size_t i = 0; i < tiles.size(); ++i)
				{
					if (i % numThreads == threadId)
					{
						renderTile(tiles[i], threadId);
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
		//using denoiser
		DenoiseFilm(film);
		//draw the scene
		for (int y = 0; y < film->height; ++y)
		{
			for (int x = 0; x < film->width; ++x)
			{
				unsigned char r, g, b;
				film->tonemap(x, y, r, g, b);
				canvas->draw(x, y, r, g, b);
			}
		}
	    //every 4 renders refresh
		if (getSPP() % 4 == 0)
		{
			float totalVariance = 0.0f;
			//comnpute variance
			for (auto& tile : tiles)
			{
				int startX = tile.tileX * TILE_SIZE;
				int startY = tile.tileY * TILE_SIZE;
				int endX = min(startX + TILE_SIZE, (int)film->width);
				int endY = min(startY + TILE_SIZE, (int)film->height);
				int pixelCount = 0;
				float variance = 0.0f;
				for (int y = startY; y < endY; ++y)
				{
					for (int x = startX; x < endX; ++x)
					{
						int idx = y * film->width + x;
						int spp = max(1, film->sppBuffer[idx]);
						Colour E = film->film[y * film->width + x] * (1.0f / (float)spp);
						Colour I = Colour(0.0f, 0.0f, 0.0f);
						float weightSum = 0.0f;
						int size = film->filter->size();
						for (int dy = -size; dy <= size; dy++)
						{
							for (int dx = -size; dx <= size; dx++)
							{
								int nx = x + dx;
								int ny = y + dy;
								if (nx >= 0 && nx < film->width && ny >= 0 && ny < film->height)
								{
									float w = film->filter->filter(dx, dy);
									Colour neighbor = film->film[ny * film->width + nx] * (1.0f / (float)spp);
									I = I + neighbor * w;
									weightSum += w;
								}
							}
						}
						float E_lum = 0.2126f * E.r + 0.7152f * E.g + 0.0722f * E.b;
						float I_lum = 0.2126f * I.r + 0.7152f * I.g + 0.0722f * I.b;
						float diff = E_lum - I_lum;
						variance += diff * diff;
						pixelCount++;
					}
				}
				if (pixelCount > 1)
					variance /= (pixelCount - 1);
				else
					variance = 0.0f;
				tile.variance = variance;
				totalVariance += variance;
			}
			//compute weight based on variance
			for (auto& tile : tiles)
			{
				if (totalVariance > 0.0f)
					tile.weight = tile.variance / totalVariance;
				else
					tile.weight = 1.0f / tiles.size();
			}
			//refresh all tiles
			int totalSamplesPerFrame = film->width * film->height; 
			for (auto& tile : tiles)
			{
				float tileSampleBudget = tile.weight * totalSamplesPerFrame;
				int tileSampleCount = (int)(tileSampleBudget); 
				tile.targetSPP = max(1, tileSampleCount / (TILE_SIZE * TILE_SIZE));
			}
		}
	}

	//This render using Instant Radiosity, other methods are same as the previous one

	//void render()
	//{
	//	static const int TILE_SIZE = 32;
	//	//memory the number of render
	//	film->incrementSPP();
	//	std::vector<std::thread> workers;
	//	std::vector<Tile> tiles;
	//	int numThreads = numProcs;
	//	int numTilesX = (film->width + TILE_SIZE - 1) / TILE_SIZE;
	//	int numTilesY = (film->height + TILE_SIZE - 1) / TILE_SIZE;
	//	//initialize all tiles
	//	for (int tileY = 0; tileY < numTilesY; ++tileY)
	//	{
	//		for (int tileX = 0; tileX < numTilesX; ++tileX)
	//		{
	//			Tile tile;
	//			tile.tileX = tileX;
	//			tile.tileY = tileY;
	//			tiles.push_back(tile);
	//		}
	//	}
	//	//precompute VPLs
	//	traceVPLs(samplers, 16);
	//	//compute one tile
	//	auto renderTile = [&](Tile& tile, int threadId)
	//		{
	//			int startX = tile.tileX * TILE_SIZE;
	//			int startY = tile.tileY * TILE_SIZE;
	//			int endX = min(startX + TILE_SIZE, (int)film->width);
	//			int endY = min(startY + TILE_SIZE, (int)film->height);
	//			for (int s = 0; s < tile.targetSPP; ++s)
	//			{
	//				for (int y = startY; y < endY; y++)
	//				{
	//					for (int x = startX; x < endX; x++)
	//					{
	//						float px = x + samplers->next();
	//						float py = y + samplers->next();
	//						Ray ray = scene->camera.generateRay(px, py);
	//						Colour A;
	//						Vec3 N;
	//						//change to ir methods
	//						Colour col = computeForPixel(ray, &samplers[threadId],A,N);
	//						film->splat(px, py, col);
	//						//use for denoising
	//						film->AOV(int(px), int(py), A, N);
	//					}
	//				}
	//			}
	//			tile.currentSPP += tile.targetSPP;
	//		};
	//	//traverse all tiles
	//	auto workerFunc = [&](int threadId)
	//		{
	//			for (size_t i = 0; i < tiles.size(); ++i)
	//			{
	//				if (i % numThreads == threadId)
	//				{
	//					renderTile(tiles[i], threadId);
	//				}
	//			}
	//		};
	//	for (int i = 0; i < numThreads; i++)
	//	{
	//		workers.emplace_back(workerFunc, i);
	//	}
	//	for (auto& worker : workers)
	//	{
	//		worker.join();
	//	}
	//	//using denoiser
	//	DenoiseFilm(film);
	//	//draw the scene
	//	for (int y = 0; y < film->height; ++y)
	//	{
	//		for (int x = 0; x < film->width; ++x)
	//		{
	//			unsigned char r, g, b;
	//			film->tonemap(x, y, r, g, b);
	//			canvas->draw(x, y, r, g, b);
	//		}
	//	}
	//	//every 4 renders refresh
	//	if (getSPP() % 4 == 0)
	//	{
	//		float totalVariance = 0.0f;
	//		//comnpute variance
	//		for (auto& tile : tiles)
	//		{
	//			int startX = tile.tileX * TILE_SIZE;
	//			int startY = tile.tileY * TILE_SIZE;
	//			int endX = min(startX + TILE_SIZE, (int)film->width);
	//			int endY = min(startY + TILE_SIZE, (int)film->height);
	//			int pixelCount = 0;
	//			float variance = 0.0f;
	//			for (int y = startY; y < endY; ++y)
	//			{
	//				for (int x = startX; x < endX; ++x)
	//				{
	//					int idx = y * film->width + x;
	//					int spp = max(1, film->sppBuffer[idx]); 
	//					Colour E = film->film[y * film->width + x] * (1.0f / (float)spp);
	//					Colour I = Colour(0.0f, 0.0f, 0.0f);
	//					float weightSum = 0.0f;
	//					int size = film->filter->size();
	//					for (int dy = -size; dy <= size; dy++)
	//					{
	//						for (int dx = -size; dx <= size; dx++)
	//						{
	//							int nx = x + dx;
	//							int ny = y + dy;
	//							if (nx >= 0 && nx < film->width && ny >= 0 && ny < film->height)
	//							{
	//								float w = film->filter->filter(dx, dy);
	//								Colour neighbor = film->film[ny * film->width + nx] * (1.0f / (float)spp);
	//								I = I + neighbor * w;
	//								weightSum += w;
	//							}
	//						}
	//					}
	//					float E_lum = 0.2126f * E.r + 0.7152f * E.g + 0.0722f * E.b;
	//					float I_lum = 0.2126f * I.r + 0.7152f * I.g + 0.0722f * I.b;
	//					float diff = E_lum - I_lum;
	//					variance += diff * diff;
	//					pixelCount++;
	//				}
	//			}
	//			if (pixelCount > 1)
	//				variance /= (pixelCount - 1);
	//			else
	//				variance = 0.0f;
	//			tile.variance = variance;
	//			totalVariance += variance;
	//		}
	//		//compute weight based on variance
	//		for (auto& tile : tiles)
	//		{
	//			if (totalVariance > 0.0f)
	//				tile.weight = tile.variance / totalVariance;
	//			else
	//				tile.weight = 1.0f / tiles.size(); 
	//		}
	//		//refresh all tiles
	//		int totalSamplesPerFrame = film->width * film->height;
	//		for (auto& tile : tiles)
	//		{
	//			float tileSampleBudget = tile.weight * totalSamplesPerFrame;
	//			int tileSampleCount = (int)(tileSampleBudget);
	//			tile.targetSPP = max(1, tileSampleCount / (TILE_SIZE * TILE_SIZE));
	//		}
	//	}
	//}

    //	// 多线程
//void render()
//{
//	film->incrementSPP();
//	std::vector<std::thread> workers;
//	int numThreads = numProcs;
//
//	// 定义每帧总路径数：可设置为图像像素数量的倍数（比如每像素发射 1 条路径）
//	const int totalPaths = film->width * film->height;
//	const int pathsPerThread = totalPaths / numThreads;
//
//	// 每个线程执行的逻辑
//	auto workerFunc = [&](int threadId)
//		{
//			for (int i = 0; i < pathsPerThread; ++i)
//			{
//				// 使用线程 id + 局部 i 生成 sampler（或用共享采样器数组）
//				Colour A;
//				Vec3 N;
//
//				// 发射一条从光源出发的路径
//				lightTrace(&samplers[threadId], A, N);
//			}
//		};
//
//	// 启动线程并加入
//	for (int i = 0; i < numThreads; i++)
//	{
//		workers.emplace_back(workerFunc, i);
//	}
//	for (auto& worker : workers)
//	{
//		worker.join();
//	}
//
//	// 后处理：降噪 + 色调映射 + 绘制
//	for (int y = 0; y < film->height; ++y)
//	{
//		for (int x = 0; x < film->width; ++x)
//		{
//			unsigned char r, g, b;
//			film->tonemap(x, y, r, g, b);
//			canvas->draw(x, y, r, g, b);
//		}
//	}
//}

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