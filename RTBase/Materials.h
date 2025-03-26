#pragma once

#include "Core.h"
#include "Imaging.h"
#include "Sampling.h"

#pragma warning( disable : 4244)

class BSDF;

class ShadingData
{
public:
	Vec3 x;
	Vec3 wo;
	Vec3 sNormal;
	Vec3 gNormal;
	float tu;
	float tv;
	Frame frame;
	BSDF* bsdf;
	float t;
	ShadingData() {}
	ShadingData(Vec3 _x, Vec3 n)
	{
		x = _x;
		gNormal = n;
		sNormal = n;
		bsdf = NULL;
	}
};

class ShadingHelper
{
public:
	static float fresnelDielectric(float cosTheta, float iorInt, float iorExt)
	{
		if (cosTheta < -1.0f) cosTheta = -1.0f;
		if (cosTheta > 1.0f) cosTheta = 1.0f;

		bool entering = cosTheta > 0.0f;

		// 根据入射方向选择 etaI / etaT
		float etaI = entering ? iorExt : iorInt;
		float etaT = entering ? iorInt : iorExt;

		// 修正 cosThetaI
		float cosThetaI = std::abs(cosTheta);

		// 计算 sin²θT
		float sinThetaTSq = (etaI / etaT) * (etaI / etaT) * (1.0f - cosThetaI * cosThetaI);

		// 全反射判定
		if (sinThetaTSq >= 1.0f)
			return 1.0f;

		float cosThetaT = std::sqrt(std::max(0.0f, 1.0f - sinThetaTSq));

		// 反射系数 Rs / Rp
		float Rs = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
			((etaT * cosThetaI) + (etaI * cosThetaT));
		float Rp = ((etaI * cosThetaI) - (etaT * cosThetaT)) /
			((etaI * cosThetaI) + (etaT * cosThetaT));

		return (Rs * Rs + Rp * Rp) * 0.5f;
	}
	static Colour fresnelConductor(float cosTheta, Colour ior, Colour k)
	{
		// Add code here
		Colour cos2 = Colour(cosTheta * cosTheta, cosTheta * cosTheta, cosTheta * cosTheta);
		Colour sin2 = Colour(1.0f, 1.0f, 1.0f) - cos2;

		Colour ior2 = ior * ior;
		Colour k2 = k * k;

		Colour twoEtaCos = ior * (2.0f * cosTheta);

		Colour t0 = ior2 + k2;

		Colour one = Colour(1.0f, 1.0f, 1.0f);

		Colour rParallel2 = ((t0 * cos2 - twoEtaCos + one) /
			(t0 * cos2 + twoEtaCos + one));

		Colour rPerpendicular2 = ((t0 - twoEtaCos + cos2) /
			(t0 + twoEtaCos + cos2));

		return (rParallel2 + rPerpendicular2) * 0.5f;
	}
	//这是 Smith 的 λ 函数，用于简化遮挡因子的计算
	static float lambdaGGX(Vec3 wi, float alpha)
	{
		// Add code here
		float cosTheta = wi.z;
		float sinTheta = std::sqrt(std::max(0.0f, 1.0f - cosTheta * cosTheta));
		float tanTheta = sinTheta / std::max(1e-4f, cosTheta);
		float a = alpha * tanTheta;

		return (-1.0f + std::sqrt(1.0f + a * a)) * 0.5f;
	}
	static float Gggx(Vec3 wi, Vec3 wo, float alpha)
	{
		// Add code here
		return 1.0f / (1.0f + lambdaGGX(wi, alpha) + lambdaGGX(wo, alpha));
	}
	//GGX 的法线分布函数，用于控制微表面朝向的统计密度
	static float Dggx(Vec3 h, float alpha)
	{
		// Add code here
		float cosThetaH = h.z;
		if (cosThetaH <= 0.0f) return 0.0f;

		float alpha2 = alpha * alpha;
		float cosThetaH2 = cosThetaH * cosThetaH;

		float denom = cosThetaH2 * (alpha2 - 1.0f) + 1.0f;
		denom = M_PI * denom * denom;

		return alpha2 / denom;
	}
};

class BSDF
{
public:
	Colour emission;
	virtual Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf) = 0;
	virtual Colour evaluate(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual float PDF(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual bool isPureSpecular() = 0;
	virtual bool isTwoSided() = 0;
	bool isLight()
	{
		return emission.Lum() > 0 ? true : false;
	}
	void addLight(Colour _emission)
	{
		emission = _emission;
	}
	Colour emit(const ShadingData& shadingData, const Vec3& wi)
	{
		return emission;
	}
	virtual float mask(const ShadingData& shadingData) = 0;
};

//漫反射
class DiffuseBSDF : public BSDF
{
public:
	Texture* albedo;
	DiffuseBSDF() = default;
	DiffuseBSDF(Texture* _albedo)
	{
		albedo = _albedo;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// Add correct sampling code here
		Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(),sampler->next());
		pdf = wi.z/M_PI;
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		wi = shadingData.frame.toWorld(wi);
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		return albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Add correct PDF code here
		Vec3 temp = shadingData.frame.toLocal(wi);
		return SamplingDistributions::cosineHemispherePDF(temp);
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};
//镜面反射，完美反射
class MirrorBSDF : public BSDF
{
public:
	Texture* albedo;
	MirrorBSDF() = default;
	MirrorBSDF(Texture* _albedo)
	{
		albedo = _albedo;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// Replace this with Mirror sampling code
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
		Vec3 wi = shadingData.frame.toWorld(wiLocal);
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv);
		pdf = 1.0f;
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Mirror evaluation code
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		Vec3 idealReflection = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
		if (Dot(wiLocal, idealReflection) > 0.999f)
		{
			float cosTheta = std::abs(wiLocal.z);
			return albedo->sample(shadingData.tu, shadingData.tv) / cosTheta;
		}
		return Colour(0.0f, 0.0f, 0.0f);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Mirror PDF
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		Vec3 idealReflection = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
		if (Dot(wiLocal, idealReflection) > 0.999f)
			return 1.0f;
		return 0.0f;
	}
	bool isPureSpecular()
	{
		return true;
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};
//金属反射
class ConductorBSDF : public BSDF
{
public:
	Texture* albedo;
	Colour eta;
	Colour k;
	float alpha;
	ConductorBSDF() = default;
	ConductorBSDF(Texture* _albedo, Colour _eta, Colour _k, float roughness)
	{
		albedo = _albedo;
		eta = _eta;
		k = _k;
		alpha = 1.62142f * sqrtf(roughness);
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// Replace this with Conductor sampling code
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		if (woLocal.z <= 0.0f)
		{
			pdf = 0.0f;
			reflectedColour = Colour(0.0f, 0.0f, 0.0f);
			return Vec3(0.0f, 0.0f, 0.0f);
		}

		// 采样半角向量（微表面法线）
		float u1 = sampler->next();
		float u2 = sampler->next();
		float phi = 2.0f * M_PI * u1;
		float cosTheta = std::sqrt((1.0f - u2) / (1.0f + (alpha * alpha - 1.0f) * u2));
		float sinTheta = std::sqrt(std::max(0.0f, 1.0f - cosTheta * cosTheta));
		Vec3 hLocal = Vec3(std::cos(phi) * sinTheta, std::sin(phi) * sinTheta, cosTheta).normalize();

		Vec3 wiLocal = Vec3(
			-woLocal.x + 2.0f * hLocal.x * Dot(woLocal, hLocal),
			-woLocal.y + 2.0f * hLocal.y * Dot(woLocal, hLocal),
			-woLocal.z + 2.0f * hLocal.z * Dot(woLocal, hLocal)
		);
		if (wiLocal.z <= 0.0f)
		{
			pdf = 0.0f;
			reflectedColour = Colour(0.0f, 0.0f, 0.0f);
			return Vec3(0.0f, 0.0f, 0.0f);
		}

		// 变换到世界空间
		Vec3 wi = shadingData.frame.toWorld(wiLocal);

		// Fresnel
		float cosThetaH = std::max(0.0f, Dot(wiLocal, hLocal));
		Colour F = ShadingHelper::fresnelConductor(cosThetaH, eta, k);

		// D、G
		float D = ShadingHelper::Dggx(hLocal, alpha);
		float G = ShadingHelper::Gggx(wiLocal, woLocal, alpha);

		// 最终 BSDF 计算
		Colour base = albedo->sample(shadingData.tu, shadingData.tv);
		float denom = 4.0f * std::max(1e-4f, wiLocal.z * woLocal.z);
		reflectedColour = base * F * (D * G / denom);

		// PDF 计算
		float dotWH = std::max(Dot(wiLocal, hLocal), 1e-4f);
		float dotNH = std::max(hLocal.z, 1e-4f);
		pdf = D * dotNH / (4.0f * dotWH);

		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Conductor evaluation code
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		if (woLocal.z <= 0.0f || wiLocal.z <= 0.0f)
			return Colour(0.0f, 0.0f, 0.0f);

		Vec3 h = (wiLocal + woLocal).normalize();
		if (Dot(woLocal, h) <= 0.0f) return Colour(0.0f, 0.0f, 0.0f);

		// Fresnel
		float cosThetaH = std::max(0.0f, Dot(wiLocal, h));
		Colour F = ShadingHelper::fresnelConductor(cosThetaH, eta, k);

		// D、G
		float D = ShadingHelper::Dggx(h, alpha);
		float G = ShadingHelper::Gggx(wiLocal, woLocal, alpha);

		// 最终 BSDF
		Colour base = albedo->sample(shadingData.tu, shadingData.tv);
		float denom = 4.0f * std::max(1e-4f, wiLocal.z * woLocal.z);
		return base * F * (D * G / denom);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Conductor PDF
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		if (woLocal.z <= 0.0f || wiLocal.z <= 0.0f)
			return 0.0f;

		Vec3 h = (wiLocal + woLocal).normalize();
		if (Dot(woLocal, h) <= 0.0f) return 0.0f;

		float D = ShadingHelper::Dggx(h, alpha);
		float dotWH = std::max(Dot(wiLocal, h), 1e-4f);
		float dotNH = std::max(h.z, 1e-4f);

		return D * dotNH / (4.0f * dotWH);
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};
//玻璃球
class GlassBSDF : public BSDF
{
public:
	Texture* albedo;
	float intIOR;
	float extIOR;
	GlassBSDF() = default;
	GlassBSDF(Texture* _albedo, float _intIOR, float _extIOR)
	{
		albedo = _albedo;
		intIOR = _intIOR;
		extIOR = _extIOR;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// Replace this with Glass sampling code
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		bool entering = woLocal.z > 0.0f;
		float etaI = entering ? extIOR : intIOR;
		float etaT = entering ? intIOR : extIOR;
		float eta = etaI / etaT;

		float cosThetaI = std::abs(woLocal.z);
		float fresnel = ShadingHelper::fresnelDielectric(cosThetaI, etaI, etaT);

		bool reflectEvent = sampler->next() < fresnel;

		Vec3 wiLocal;
		if (reflectEvent)
		{
			// Perfect mirror reflection
			wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv);
			pdf = fresnel;
		}
		else
		{
			// Perfect transmission (Snell's law)
			float sinTheta2 = eta * eta * (1.0f - cosThetaI * cosThetaI);
			if (sinTheta2 > 1.0f)
			{
				// Total internal reflection fallback
				wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
				reflectedColour = albedo->sample(shadingData.tu, shadingData.tv);
				pdf = 1.0f;
				return shadingData.frame.toWorld(wiLocal);
			}

			float cosThetaT = std::sqrt(1.0f - sinTheta2);
			cosThetaT = entering ? -cosThetaT : cosThetaT;

			wiLocal = Vec3(
				-eta * woLocal.x,
				-eta * woLocal.y,
				cosThetaT
			);

			// Radiance scaling
			float scale = (eta * eta);
			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * scale * (1.0f - fresnel);
			pdf = 1.0f - fresnel;
		}

		return shadingData.frame.toWorld(wiLocal);
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Glass evaluation code
		return Colour(0.0f, 0.0f, 0.0f);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with GlassPDF
		return 0.0f;
	}
	bool isPureSpecular()
	{
		return true;
	}
	bool isTwoSided()
	{
		return false;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class DielectricBSDF : public BSDF
{
public:
	Texture* albedo;
	float intIOR;
	float extIOR;
	float alpha;
	DielectricBSDF() = default;
	DielectricBSDF(Texture* _albedo, float _intIOR, float _extIOR, float roughness)
	{
		albedo = _albedo;
		intIOR = _intIOR;
		extIOR = _extIOR;
		alpha = 1.62142f * sqrtf(roughness);
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// Replace this with Dielectric sampling code
		Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		pdf = wi.z / M_PI;
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		wi = shadingData.frame.toWorld(wi);
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Dielectric evaluation code
		return albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Dielectric PDF
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		return SamplingDistributions::cosineHemispherePDF(wiLocal);
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return false;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class OrenNayarBSDF : public BSDF
{
public:
	Texture* albedo;
	float sigma;
	OrenNayarBSDF() = default;
	OrenNayarBSDF(Texture* _albedo, float _sigma)
	{
		albedo = _albedo;
		sigma = _sigma;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// Replace this with OrenNayar sampling code
		Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		pdf = wi.z / M_PI;
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		wi = shadingData.frame.toWorld(wi);
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with OrenNayar evaluation code
		return albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with OrenNayar PDF
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		return SamplingDistributions::cosineHemispherePDF(wiLocal);
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class PlasticBSDF : public BSDF
{
public:
	Texture* albedo;
	float intIOR;
	float extIOR;
	float alpha;
	PlasticBSDF() = default;
	PlasticBSDF(Texture* _albedo, float _intIOR, float _extIOR, float roughness)
	{
		albedo = _albedo;
		intIOR = _intIOR;
		extIOR = _extIOR;
		alpha = 1.62142f * sqrtf(roughness);
	}
	float alphaToPhongExponent()
	{
		return (2.0f / SQ(std::max(alpha, 0.001f))) - 2.0f;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// Replace this with Plastic sampling code
		Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		pdf = wi.z / M_PI;
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		wi = shadingData.frame.toWorld(wi);
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Plastic evaluation code
		return albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Plastic PDF
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		return SamplingDistributions::cosineHemispherePDF(wiLocal);
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class LayeredBSDF : public BSDF
{
public:
	BSDF* base;
	Colour sigmaa;
	float thickness;
	float intIOR;
	float extIOR;
	LayeredBSDF() = default;
	LayeredBSDF(BSDF* _base, Colour _sigmaa, float _thickness, float _intIOR, float _extIOR)
	{
		base = _base;
		sigmaa = _sigmaa;
		thickness = _thickness;
		intIOR = _intIOR;
		extIOR = _extIOR;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// Add code to include layered sampling
		return base->sample(shadingData, sampler, reflectedColour, pdf);
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Add code for evaluation of layer
		return base->evaluate(shadingData, wi);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Add code to include PDF for sampling layered BSDF
		return base->PDF(shadingData, wi);
	}
	bool isPureSpecular()
	{
		return base->isPureSpecular();
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return base->mask(shadingData);
	}
};