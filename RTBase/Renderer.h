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

struct Tile {
	int tileX, tileY;
	float variance = 0.0f;
	//本轮权重
	float weight = 1.0f;
	int currentSPP = 0;
	//本轮采样数
	int targetSPP = 1;
};

struct VPL {
	ShadingData shadingData;
	Colour Le;
	bool fromLight = false;
	bool isEnvironment = false; // 该 VPL 是否来自环境光
	Vec3 envDir;               // 若 isEnvironment = true，则此处记录环境光方向
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

	//下面连续三个方法用于instant
	const size_t MAX_VPLS = 100;

	void VPLTracePath(Ray& r, Colour pathThroughput, Colour Le, Sampler* sampler, int depth) {
		// 如果已有 VPL 数量足够，就直接返回
		if (vpls.size() >= MAX_VPLS) return;

		// 如果这是环境光的 VPL，就不追踪了（可按需选择是否写这条判断）
		// （本函数一般是追踪物体表面反弹产生的 VPL，环境光可以跳过）
		// 你也可以把这句逻辑放在 traceVPLs 里来控制
		// -------------------------------------
		// if (/* 某些情况表示这是环境光 */) { return; }
		// -------------------------------------

		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);

		if (shadingData.t < FLT_MAX) {
			if (depth > 5) {
				return;
			}

			// 镜面直接返回，不用在此生成 VPL
			if (shadingData.bsdf->isPureSpecular() == true) {
				return;
			}

			// 生成一个 VPL
			if (vpls.size() < MAX_VPLS) {
				VPL vpl;
				vpl.shadingData = shadingData;
				vpl.Le = pathThroughput * Le;
				vpl.fromLight = false;      // 此处是物体反射 VPL
				vpl.isEnvironment = false;  // 普通表面，不是环境
				vpls.push_back(vpl);
			}
			else {
				return;
			}

			// 以下是继续向场景追踪，采样新的方向
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

	void traceVPLs(Sampler* sampler, int N_VPLs) {
		vpls.clear();
		vpls.reserve(MAX_VPLS);

		for (int i = 0; i < N_VPLs; ++i) {
			float pmf;
			Light* light = scene->sampleLight(sampler, pmf);
			if (pmf <= 0.0f) {
				continue;
			}

			// ---------------------------
			// 如果是环境光：
			// ---------------------------
			if (light->isENV()) {
				float pdfDir;
				Vec3 envDir = light->sampleDirectionFromLight(sampler, pdfDir);
				if (pdfDir <= 0.0f) {
					continue;
				}

				// 该方向的环境辐射
				Colour envEmit = light->getEmitted(envDir);
				// 归一化
				Colour LE = envEmit / (pmf * pdfDir * float(N_VPLs));

				// 构建一个环境 VPL
				VPL vpl;
				// 由于是环境光，不需要真实位置和法线
				vpl.shadingData.x = Vec3(0.0f, 0.0f, 0.0f);
				vpl.shadingData.sNormal = Vec3(0.0f, 0.0f, 0.0f);
				vpl.Le = LE;
				vpl.fromLight = true;
				vpl.isEnvironment = true; // 重点：标记为环境光
				vpl.envDir = envDir;      // 记录该方向

				vpls.push_back(vpl);

				// 通常对环境光就到此为止，不再用 VPLTracePath 继续追踪
				// 因为它不是一个“表面”可以再反弹
				// 如果你要做奇巧淫技，也可再发射射线，但这里先不做。
			}
			// ---------------------------
			// 如果是普通面积光或其他光：
			// ---------------------------
			else {
				float pdfPosition;
				Vec3 p = light->samplePositionFromLight(sampler, pdfPosition);
				if (pdfPosition <= 0.0f) {
					continue;
				}

				ShadingData shadingData(p, light->normal(p));
				Colour emitted = light->getEmitted(p);

				Colour LE = emitted / (pmf * pdfPosition * float(N_VPLs));

				// 第一个 VPL
				{
					VPL vpl;
					vpl.shadingData = shadingData;
					vpl.Le = LE;
					vpl.fromLight = true;
					vpl.isEnvironment = false; // 普通表面
					vpls.push_back(vpl);
				}

				// 从该光源点再采样方向
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
				// 这里才进入 VPLTracePath，去追踪后续反射产生的 VPL
				VPLTracePath(r, pathThroughput, le, sampler, 0);
			}
		}
	}

	Colour computeForPixel(Ray& r, Sampler* sampler, Colour& a, Vec3& N) {
		Colour result(0.0f, 0.0f, 0.0f);

		// 主射线求交
		IntersectionData intersection = scene->traverse(r);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		N = shadingData.sNormal.normalize();

		if(shadingData.bsdf!=NULL)
		if (shadingData.bsdf->isLight()) {
			return shadingData.bsdf->emission;
		}

		if (shadingData.t < FLT_MAX) {
			// 如果是纯镜面，仍然用你现有的 pathTrace
			if (shadingData.bsdf->isPureSpecular() == true) {
				Colour throughput(1.0f,1.0f,1.0f);
				return pathTrace(r, throughput, 0, sampler, a, N);
			}

			// 遍历所有 VPL
			for (const VPL& vpl : vpls) {
				if (vpl.isEnvironment) {
					// --------- 环境光 VPL 的计算 ---------
					Vec3 wi = vpl.envDir;
					float cosTheta_hit = max(0.0f, Dot(wi, shadingData.sNormal));
					if (cosTheta_hit < 1e-4f) {
						continue; // 太背面，不贡献
					}

					// 无穷远可见性检测：从当前点沿 wi 发射射线，看是否击中场景
					Ray testRay(shadingData.x + wi * EPSILON, wi);
					IntersectionData testIsect = scene->traverse(testRay);
					// 若 testIsect.t 很大或没碰到物体，则视作可见环境
					bool visibleEnv = (testIsect.t >= FLT_MAX);

					if (visibleEnv) {
						// BSDF 求值
						Colour bsdfVal = shadingData.bsdf->evaluate(shadingData, wi);
						// 一般乘一次 cosTheta_hit
						Colour contrib = vpl.Le * bsdfVal * cosTheta_hit;
						result = result + contrib;
					}
				}
				else {
					// --------- 普通 VPL (面积光 or 物体反射) 的计算 ---------
					Vec3 vplPos = vpl.shadingData.x;
					Vec3 toVPL = vplPos - shadingData.x;
					float dist2 = Dot(toVPL, toVPL);
					if (dist2 < 1e-8f) {
						continue;
					}
					Vec3 wi = toVPL.normalize();

					float cosTheta_hit = max(0.0f, Dot(wi, shadingData.sNormal));
					float cosTheta_vpl = max(0.0f, -Dot(wi, vpl.shadingData.sNormal));
					if (cosTheta_hit < 0.1f || cosTheta_vpl < 0.1f) {
						continue;
					}

					float G = (cosTheta_hit * cosTheta_vpl) / (dist2 + EPSILON);
					if (G < 1e-4f) {
						continue;
					}

					// 点到点可见性测试
					if (!scene->visible(shadingData.x, vplPos)) {
						continue;
					}

					Colour temp(0.0f, 0.0f, 0.0f);
					if (vpl.fromLight) {
						// 若是从光源直接采样的 VPL
						temp = vpl.Le
							* shadingData.bsdf->evaluate(shadingData, wi)
							* G;
					}
					else {
						// 若是物体反射生成的 VPL
						temp = vpl.Le
							* vpl.shadingData.bsdf->evaluate(vpl.shadingData, -wi)
							* shadingData.bsdf->evaluate(shadingData, wi)
							* G;
					}

					// 你原本可能还有对 shadingData.bsdf->sample(...) 的调用
					// 是否继续保留取决于你的需求，这里省略。

					result = result + temp;
				}
			}

			return result;
		}

		// 未击中任何几何体则返回背景
		return scene->background->evaluate(shadingData, r.dir);
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

	// //多线程，适应性采样版
	//void render()
	//{
	//	static const int TILE_SIZE = 32;
	//	film->incrementSPP();
	//	std::vector<std::thread> workers;
	//	std::vector<Tile> tiles;
	//	int numThreads = numProcs;
	//	int numTilesX = (film->width + TILE_SIZE - 1) / TILE_SIZE;
	//	int numTilesY = (film->height + TILE_SIZE - 1) / TILE_SIZE;
	//	// 初始化所有 tile 信息
	//	for (int tileY = 0; tileY < numTilesY; ++tileY)
	//	{
	//		for (int tileX = 0; tileX < numTilesX; ++tileX)
	//		{
	//			Tile tile;
	//			tile.tileX = tileX;
	//			tile.tileY = tileY;
	//			// 其他字段（如 variance、weight、currentSPP、targetSPP）可以在 Tile 结构体中预置默认值
	//			tiles.push_back(tile);
	//		}
	//	}
	//	// 渲染每个 tile 的函数：每个像素采 INIT_SPP 个样本（此处为 1 个）
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
	//						Colour pathThroughput(1.0f, 1.0f, 1.0f);
	//						Colour A;
	//						Vec3 N;
	//						Colour col = pathTrace(ray, pathThroughput, 0, &samplers[threadId],A,N);
	//						film->splat(px, py, col);
	//						film->AOV(int(px),int(py),A,N);
	//					}
	//				}
	//			}
	//			tile.currentSPP += tile.targetSPP;
	//		};
	//	// 多线程遍历 tiles 数组（按 round-robin 分配）
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
	//	DenoiseFilm(film);
	//	for (int y = 0; y < film->height; ++y)
	//	{
	//		for (int x = 0; x < film->width; ++x)
	//		{
	//			unsigned char r, g, b;
	//			film->tonemap(x, y, r, g, b);
	//			canvas->draw(x, y, r, g, b);
	//		}
	//	}
	//	// 每 4 帧更新一次方差（仅计算并输出调试信息）
	//	if (getSPP() % 4 == 0)
	//	{
	//		float totalVariance = 0.0f;
	//		//计算方差
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
	//					//历史累计值
	//					int idx = y * film->width + x;
	//					int spp = max(1, film->sppBuffer[idx]); // 防止除 0
	//					Colour E = film->film[y * film->width + x] * (1.0f / (float)spp);
	//					//高斯近似计算当前帧值
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
	//		//由方差计算权重
	//		for (auto& tile : tiles)
	//		{
	//			if (totalVariance > 0.0f)
	//				tile.weight = tile.variance / totalVariance;
	//			else
	//				tile.weight = 1.0f / tiles.size(); // fallback：平均采样
	//		}
	//		//调整tile结构体
	//		int totalSamplesPerFrame = film->width * film->height; 
	//		for (auto& tile : tiles)
	//		{
	//			float tileSampleBudget = tile.weight * totalSamplesPerFrame;
	//			int tileSampleCount = (int)(tileSampleBudget); // 这是总共多少次采样
	//			tile.targetSPP = max(1, tileSampleCount / (TILE_SIZE * TILE_SIZE));
	//		}
	//	}
	//}

	//多线程，instant版本
	void render()
	{
		static const int TILE_SIZE = 32;
		film->incrementSPP();
		std::vector<std::thread> workers;
		std::vector<Tile> tiles;
		int numThreads = numProcs;
		int numTilesX = (film->width + TILE_SIZE - 1) / TILE_SIZE;
		int numTilesY = (film->height + TILE_SIZE - 1) / TILE_SIZE;
		// 初始化所有 tile 信息
		for (int tileY = 0; tileY < numTilesY; ++tileY)
		{
			for (int tileX = 0; tileX < numTilesX; ++tileX)
			{
				Tile tile;
				tile.tileX = tileX;
				tile.tileY = tileY;
				// 其他字段（如 variance、weight、currentSPP、targetSPP）可以在 Tile 结构体中预置默认值
				tiles.push_back(tile);
			}
		}
		traceVPLs(samplers, 32);
		//std::cout << vpls.size() << std::endl;
		// 渲染每个 tile 的函数：每个像素采 INIT_SPP 个样本（此处为 1 个）
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
							Colour A;
							Vec3 N;
							Colour col = computeForPixel(ray, &samplers[threadId],A,N);
							film->splat(px, py, col);
							film->AOV(int(px), int(py), A, N);
						}
					}
				}
				tile.currentSPP += tile.targetSPP;
			};
		// 多线程遍历 tiles 数组（按 round-robin 分配）
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
		DenoiseFilm(film);
		for (int y = 0; y < film->height; ++y)
		{
			for (int x = 0; x < film->width; ++x)
			{
				unsigned char r, g, b;
				film->tonemap(x, y, r, g, b);
				canvas->draw(x, y, r, g, b);
			}
		}
		// 每 4 帧更新一次方差（仅计算并输出调试信息）
		if (getSPP() % 4 == 0)
		{
			float totalVariance = 0.0f;
			//计算方差
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
						//历史累计值
						int idx = y * film->width + x;
						int spp = max(1, film->sppBuffer[idx]); // 防止除 0
						Colour E = film->film[y * film->width + x] * (1.0f / (float)spp);
						//高斯近似计算当前帧值
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
			//由方差计算权重
			for (auto& tile : tiles)
			{
				if (totalVariance > 0.0f)
					tile.weight = tile.variance / totalVariance;
				else
					tile.weight = 1.0f / tiles.size(); // fallback：平均采样
			}
			//调整tile结构体
			int totalSamplesPerFrame = film->width * film->height;
			for (auto& tile : tiles)
			{
				float tileSampleBudget = tile.weight * totalSamplesPerFrame;
				int tileSampleCount = (int)(tileSampleBudget); // 这是总共多少次采样
				tile.targetSPP = max(1, tileSampleCount / (TILE_SIZE * TILE_SIZE));
			}
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