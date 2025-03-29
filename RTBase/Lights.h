#pragma once

#include "Core.h"
#include "Geometry.h"
#include "Materials.h"
#include "Sampling.h"

#pragma warning( disable : 4244)

class SceneBounds
{
public:
	Vec3 sceneCentre;
	float sceneRadius;
};

//抽象成两种，区域光源和方向光，在光照计算中只分成两类，是否能得到光上一点
class Light
{
public:
	//从光源区域中随机取样一点
	virtual Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& emittedColour, float& pdf) = 0;
	//判断传入入射角度下，能否看到光源的正面，是否要返回光的颜色
	virtual Colour evaluate(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual float PDF(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual bool isArea() = 0;
	virtual Vec3 normal(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual float totalIntegratedPower() = 0;
	virtual Vec3 samplePositionFromLight(Sampler* sampler, float& pdf) = 0;
	virtual Vec3 sampleDirectionFromLight(Sampler* sampler, float& pdf) = 0;
};

class AreaLight : public Light
{
public:
	Triangle* triangle = NULL;
	Colour emission;
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& emittedColour, float& pdf)
	{
		emittedColour = emission;
		return triangle->sample(sampler, pdf);
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		if (Dot(wi, triangle->gNormal()) < 0)
		{
			return emission;
		}
		return Colour(0.0f, 0.0f, 0.0f);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		return 1.0f / triangle->area;
	}
	bool isArea()
	{
		return true;
	}
	Vec3 normal(const ShadingData& shadingData, const Vec3& wi)
	{
		return triangle->gNormal();
	}
	float totalIntegratedPower()
	{
		return (triangle->area * emission.Lum());
	}
	Vec3 samplePositionFromLight(Sampler* sampler, float& pdf)
	{
		return triangle->sample(sampler, pdf);
	}
	Vec3 sampleDirectionFromLight(Sampler* sampler, float& pdf)
	{
		// Add code to sample a direction from the light
		Vec3 wi = Vec3(0, 0, 1);
		pdf = 1.0f;
		Frame frame;
		frame.fromVector(triangle->gNormal());
		return frame.toWorld(wi);
	}
};

class BackgroundColour : public Light
{
public:
	Colour emission;
	BackgroundColour(Colour _emission)
	{
		emission = _emission;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		Vec3 wi = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
		pdf = SamplingDistributions::uniformSpherePDF(wi);
		reflectedColour = emission;
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		return emission;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		return SamplingDistributions::uniformSpherePDF(wi);
	}
	bool isArea()
	{
		return false;
	}
	Vec3 normal(const ShadingData& shadingData, const Vec3& wi)
	{
		return -wi;
	}
	float totalIntegratedPower()
	{
		return emission.Lum() * 4.0f * M_PI;
	}
	Vec3 samplePositionFromLight(Sampler* sampler, float& pdf)
	{
		Vec3 p = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
		p = p * use<SceneBounds>().sceneRadius;
		p = p + use<SceneBounds>().sceneCentre;
		pdf = 4 * M_PI * use<SceneBounds>().sceneRadius * use<SceneBounds>().sceneRadius;
		return p;
	}
	Vec3 sampleDirectionFromLight(Sampler* sampler, float& pdf)
	{
		Vec3 wi = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
		pdf = SamplingDistributions::uniformSpherePDF(wi);
		return wi;
	}
};

#define ENV_SCALE 0.18f  // ✅ 控制环境图亮度缩放（你可以随时调）

class EnvironmentMap : public Light
{
public:
	Texture* env;
	std::vector<float> cdf;
	float totalSum;

	EnvironmentMap(Texture* _env)
	{
		env = _env;
		int w = env->width;
		int h = env->height;
		cdf.resize(w * h, 0.0f);
		totalSum = 0.0f;

		float accum = 0.0f;
		for (int i = 0; i < h; i++)
		{
			float theta = M_PI * ((float)i + 0.5f) / (float)h;
			float sinTheta = sinf(theta);
			for (int j = 0; j < w; j++)
			{
				float lum = env->texels[i * w + j].Lum();
				float weight = lum * sinTheta;
				accum += weight;
				cdf[i * w + j] = accum;
			}
		}
		totalSum = cdf.back();
	}

	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		int w = env->width;
		int h = env->height;

		float r = sampler->next() * totalSum;
		auto it = std::lower_bound(cdf.begin(), cdf.end(), r);
		int idx = int(it - cdf.begin());
		if (idx >= w * h) idx = w * h - 1;

		int i = idx / w;
		int j = idx % w;

		float u = (j + 0.5f) / float(w);
		float v = (i + 0.5f) / float(h);
		float theta = v * M_PI;
		float phi = u * 2.0f * M_PI;

		float sinTheta = sinf(theta);
		Vec3 wi(
			sinTheta * cosf(phi),
			cosf(theta),
			sinTheta * sinf(phi)
		);

		float lum = env->texels[i * w + j].Lum();
		float texelSolidAngle = (2.0f * M_PI / w) * (M_PI / h) * sinTheta;
		pdf = (lum * sinTheta / totalSum) / texelSolidAngle;

		reflectedColour = evaluate(shadingData, wi);
		return wi;
	}

	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		float u = atan2f(wi.z, wi.x);
		u = (u < 0.0f) ? u + (2.0f * M_PI) : u;
		u = u / (2.0f * M_PI);
		float v = acosf(wi.y) / M_PI;

		// ✅ 加入亮度缩放，防止HDR爆白
		return env->sample(u, v) * ENV_SCALE;
	}

	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		int w = env->width;
		int h = env->height;

		float phi = atan2f(wi.z, wi.x);
		if (phi < 0.0f) phi += 2.0f * M_PI;
		float theta = acosf(wi.y);
		float sinTheta = sinf(theta);

		float u = phi / (2.0f * M_PI);
		float v = theta / M_PI;
		int j = std::min(std::max(int(u * w), 0), w - 1);
		int i = std::min(std::max(int(v * h), 0), h - 1);

		float lum = env->texels[i * w + j].Lum();
		float texelSolidAngle = (2.0f * M_PI / w) * (M_PI / h) * sinTheta;
		if (texelSolidAngle < 1e-6f) return 0.0f;

		return (lum * sinTheta / totalSum) / texelSolidAngle;
	}

	bool isArea() { return false; }

	Vec3 normal(const ShadingData& shadingData, const Vec3& wi)
	{
		return -wi;
	}

	float totalIntegratedPower()
	{
		float total = 0;
		for (int i = 0; i < env->height; i++)
		{
			float theta = M_PI * ((float)i + 0.5f) / (float)env->height;
			float st = sinf(theta);
			for (int n = 0; n < env->width; n++)
			{
				total += (env->texels[(i * env->width) + n].Lum() * st);
			}
		}
		total = total / (float)(env->width * env->height);
		return total * 4.0f * M_PI * ENV_SCALE;  // ✅ 也乘上亮度缩放
	}

	Vec3 samplePositionFromLight(Sampler* sampler, float& pdf)
	{
		Vec3 p = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
		p = p * use<SceneBounds>().sceneRadius;
		p = p + use<SceneBounds>().sceneCentre;
		pdf = 1.0f / (4 * M_PI * SQ(use<SceneBounds>().sceneRadius));
		return p;
	}

	Vec3 sampleDirectionFromLight(Sampler* sampler, float& pdf)
	{
		Colour dummy;
		return sample(ShadingData(), sampler, dummy, pdf);
	}
};