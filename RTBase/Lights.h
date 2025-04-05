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
	virtual bool isENV() = 0;
	virtual Vec3 normal(const Vec3& wi) = 0;
	virtual Colour getEmitted(const Vec3& wi) = 0;
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
	bool isENV()
	{
		return false;
	}
	Vec3 normal(const Vec3& wi)
	{
		return triangle->gNormal();
	}
	Colour getEmitted(const Vec3& wi)
	{
		return emission;
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
		// 使用余弦加权半球采样
		float r1 = sampler->next();
		float r2 = sampler->next();
		float phi = 2.0f * M_PI * r1;
		float theta = acos(sqrt(1.0f - r2)); // 余弦加权采样

		float x = sin(theta) * cos(phi);
		float y = sin(theta) * sin(phi);
		float z = cos(theta);

		// 构建局部坐标系，以光源法线为 Z 轴
		Frame frame;
		frame.fromVector(triangle->gNormal());
		pdf = z / M_PI;  // 余弦半球采样的 pdf

		return frame.toWorld(Vec3(x, y, z));
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
	bool isENV()
	{
		return false;
	}
	Vec3 normal( const Vec3& wi)
	{
		return -wi;
	}
	Colour getEmitted(const Vec3& wi)
	{
		return emission;
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

//Old environmentmap

//class EnvironmentMap : public Light
//{
//public:
//	Texture* env;
//	EnvironmentMap(Texture* _env)
//	{
//		env = _env;
//	}
//	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
//	{
//		// Assignment: Update this code to importance sampling lighting based on luminance of each pixel
//		Vec3 wi = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
//		pdf = SamplingDistributions::uniformSpherePDF(wi);
//		reflectedColour = evaluate(shadingData, wi);
//		return wi;
//	}
//	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
//	{
//		float u = atan2f(wi.z, wi.x);
//		u = (u < 0.0f) ? u + (2.0f * M_PI) : u;
//		u = u / (2.0f * M_PI);
//		float v = acosf(wi.y) / M_PI;
//		return env->sample(u, v);
//	}
//	float PDF(const ShadingData& shadingData, const Vec3& wi)
//	{
//		// Assignment: Update this code to return the correct PDF of luminance weighted importance sampling
//		return SamplingDistributions::uniformSpherePDF(wi);
//	}
//	bool isArea()
//	{
//		return false;
//	}
//	bool isENV()
//	{
//		return true;
//	}
//	Vec3 normal(const Vec3& wi)
//	{
//		return -wi;
//	}
//	Colour getEmitted(const Vec3& wi)
//	{
//		float u = atan2f(wi.z, wi.x);
//		u = (u < 0.0f) ? u + (2.0f * M_PI) : u;
//		u = u / (2.0f * M_PI);
//		float v = acosf(wi.y) / M_PI;
//		return env->sample(u, v);
//	}
//	float totalIntegratedPower()
//	{
//		float total = 0;
//		for (int i = 0; i < env->height; i++)
//		{
//			float st = sinf(((float)i / (float)env->height) * M_PI);
//			for (int n = 0; n < env->width; n++)
//			{
//				total += (env->texels[(i * env->width) + n].Lum() * st);
//			}
//		}
//		total = total / (float)(env->width * env->height);
//		return total * 4.0f * M_PI;
//	}
//	Vec3 samplePositionFromLight(Sampler* sampler, float& pdf)
//	{
//		// Samples a point on the bounding sphere of the scene. Feel free to improve this.
//		Vec3 p = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
//		p = p * use<SceneBounds>().sceneRadius;
//		p = p + use<SceneBounds>().sceneCentre;
//		pdf = 1.0f / (4 * M_PI * SQ(use<SceneBounds>().sceneRadius));
//		return p;
//	}
//	Vec3 sampleDirectionFromLight(Sampler* sampler, float& pdf)
//	{
//		// Replace this tabulated sampling of environment maps
//		Vec3 wi = SamplingDistributions::uniformSampleSphere(sampler->next(), sampler->next());
//		pdf = SamplingDistributions::uniformSpherePDF(wi);
//		return wi;
//	}
//};

//New environmentmap with importance sampling based on pixel luminance
#define ENV_SCALE 0.26f//used to control evaluate. without it, the image becomes overexposed.(Not sure why)
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
		//precompute CDF
		float accum = 0.0f;
		for (int i = 0; i < h; i++)
		{
			float theta = M_PI * ((float)i + 0.5f) / (float)h;
			float sinTheta = sinf(theta);
			for (int j = 0; j < w; j++)
			{
				float lum = env->texels[i * w + j].Lum();
				accum += lum * sinTheta;
				cdf[i * w + j] = accum;
			}
		}
		totalSum = accum;
	}

	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		int w = env->width;
		int h = env->height;

		//get a random pixel
		float r = sampler->next() * totalSum;
		auto it = std::lower_bound(cdf.begin(), cdf.end(), r);
		int idx = int(it - cdf.begin());
		if (idx >= w * h) idx = w * h - 1;
		int i = idx / w;
		int j = idx % w;
		//convert to spherical coordinates
		float u = (j + 0.5f) / float(w);
		float v = (i + 0.5f) / float(h);
		float theta = v * M_PI;
		float phi = u * 2.0f * M_PI;
		float sin = sinf(theta);
		Vec3 wi(
			sin * cosf(phi),
			cosf(theta),
			sin * sinf(phi)
		);
		//compute pdf and colour
		float lum = env->texels[i * w + j].Lum();
		float texelSolidAngle = (2.0f * M_PI / w) * (M_PI / h) * sin;
		pdf = (lum * sin / totalSum) / texelSolidAngle;
		reflectedColour = evaluate(shadingData, wi);
		return wi;
	}

	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		float u = atan2f(wi.z, wi.x);
		u = (u < 0.0f) ? u + (2.0f * M_PI) : u;
		u = u / (2.0f * M_PI);
		float v = acosf(wi.y) / M_PI;
		return env->sample(u, v) * ENV_SCALE;
	}

	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		int w = env->width;
		int h = env->height;
		//compute spherical coordinates 
		float phi = atan2f(wi.z, wi.x);
		if (phi < 0.0f) phi += 2.0f * M_PI;
		float theta = acosf(wi.y);
		float sin = sinf(theta);
		//convert to uv coordinates and get xy position
		float u = phi / (2.0f * M_PI);
		float v = theta / M_PI;
		int j = std::min(std::max(int(u * w), 0), w - 1);
		int i = std::min(std::max(int(v * h), 0), h - 1);
		//compute pdf
		float lum = env->texels[i * w + j].Lum();
		float texelSolidAngle = (2.0f * M_PI / w) * (M_PI / h) * sin;
		if (texelSolidAngle < 1e-6f) return 0.0f;
		return (lum * sin / totalSum) / texelSolidAngle;
	}

	bool isArea()
	{
		return false;
	}
	bool isENV()
	{
		return true;
	}
	Vec3 normal(const Vec3& wi)
	{
		return -wi;
	}
	Colour getEmitted(const Vec3& wi)
	{
		float u = atan2f(wi.z, wi.x);
		u = (u < 0.0f) ? u + (2.0f * M_PI) : u;
		u = u / (2.0f * M_PI);
		float v = acosf(wi.y) / M_PI;
		return env->sample(u, v) * ENV_SCALE;
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
		return total * 4.0f * M_PI * ENV_SCALE;
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
		//use sample
		Colour temp;
		return sample(ShadingData(), sampler, temp, pdf);
	}
};