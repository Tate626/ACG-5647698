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
	//返回折射方向wi（局部），传入的也是局部入射wo及折射率比
	//全反射则返回false
	bool refract(const Vec3& wo, float eta, Vec3& wi)
	{
		// cosθi：入射角余弦（在局部空间中 = wo.z）
		float cosThetaI = wo.z;

		// sin²θi = 1 - cos²θi
		float sin2ThetaI = std::max(0.0f, 1.0f - cosThetaI * cosThetaI);

		// Snell's Law：sin²θt = η² * sin²θi
		float sin2ThetaT = eta * eta * sin2ThetaI;

		// 如果大于 1，表示发生全反射（TIR）
		if (sin2ThetaT >= 1.0f)
			return false;

		// 否则计算 cosθt = sqrt(1 - sin²θt)
		float cosThetaT = std::sqrt(1.0f - sin2ThetaT);

		// 构造折射方向 wi（局部坐标系中：z-up）
		wi = Vec3(
			-eta * wo.x,
			-eta * wo.y,
			wo.z > 0 ? -cosThetaT : cosThetaT
		);
		return true;
	}
	//返回的值代表多少比例的光是反射走的（剩下就是折射）
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
	//使用schlick公式近似计算反射率，完美代替上面那个diele方法
	static float fresnelSchlick(float cosTheta, float etaI, float etaT)
	{
		float f0 = (etaI - etaT) / (etaI + etaT);
		f0 = f0 * f0;

		float oneMinusCos = 1.0f - cosTheta;
		return f0 + (1.0f - f0) * std::pow(oneMinusCos, 5.0f);
	}
	//返回金属的反射率
	static Colour fresnelConductor(float cosTheta, Colour ior, Colour k)
	{
		// Add code here
		Colour cos2 = Colour(cosTheta * cosTheta, cosTheta * cosTheta, cosTheta * cosTheta);
		Colour sin2 = Colour(1.0f, 1.0f, 1.0f) - cos2;

		Colour ior2 = ior * ior;
		Colour k2 = k * k;
		Colour t0 = ior2 + k2;
		Colour twoEtaCos = ior * (2.0f * cosTheta);

		Colour rParallel2 = ((t0 * cos2 - twoEtaCos + sin2) /
			(t0 * cos2 + twoEtaCos + sin2));

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
		//返回材质图对应位置的颜色
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
//镜面反射，完美反射,有问题，pdf到底是几，贡献除不除cos
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
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wi = Vec3(-wo.x, -wo.y, wo.z);
		float cosTheta = std::abs(Dot(wi,shadingData.gNormal));
		wi = shadingData.frame.toWorld(wi);
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv);
		pdf = 1.0f;
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Mirror evaluation code
		return Colour(0.0f, 0.0f, 0.0f);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Mirror PDF
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
//金属反射,没有折射，但是不是完美反射，部分反射
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
//电介质材料，水，冰，等
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
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		bool entering = woLocal.z > 0.0f;
		float etaI = entering ? extIOR : intIOR;
		float etaT = entering ? intIOR : extIOR;
		float eta = etaI / etaT;

		float cosThetaI = woLocal.z;
		if (cosThetaI < -1.0f) cosThetaI = -1.0f;
		else if (cosThetaI > 1.0f) cosThetaI = 1.0f;
		float fresnel = ShadingHelper::fresnelDielectric(cosThetaI, intIOR, extIOR);

		bool reflect = sampler->next() < fresnel;

		Vec3 wiLocal;
		if (reflect)
		{
			// Perfect reflection
			wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * fresnel;
			pdf = fresnel;
		}
		else
		{
			// Refraction using Snell's law
			float sin2ThetaI = std::max(0.0f, 1.0f - cosThetaI * cosThetaI);
			float sin2ThetaT = eta * eta * sin2ThetaI;

			if (sin2ThetaT >= 1.0f)
			{
				// Total internal reflection
				wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
				reflectedColour = albedo->sample(shadingData.tu, shadingData.tv);
				pdf = 1.0f;
				return shadingData.frame.toWorld(wiLocal);
			}

			float cosThetaT = std::sqrt(1.0f - sin2ThetaT);
			cosThetaT = entering ? -cosThetaT : cosThetaT;

			wiLocal = Vec3(
				-eta * woLocal.x,
				-eta * woLocal.y,
				cosThetaT
			);

			float scale = (eta * eta);
			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * scale * (1.0f - fresnel);
			pdf = 1.0f - fresnel;
		}

		return shadingData.frame.toWorld(wiLocal);
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Dielectric evaluation code
		return Colour(0.0f, 0.0f, 0.0f);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Dielectric PDF
		return 0.0f;
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
//塑料材质
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
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		float exponent = alphaToPhongExponent();

		// 随机决定是采样漫反射还是高光反射
		float choice = sampler->next();
		Vec3 wiLocal;

		if (choice < 0.5f)
		{
			// Diffuse 半球采样
			wiLocal = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		}
		else
		{
			// Phong 高光采样：采样一个围绕理想反射方向的 lobe
			Vec3 reflectDir = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
			float phi = 2.0f * M_PI * sampler->next();
			float cosTheta = powf(sampler->next(), 1.0f / (exponent + 1.0f));
			float sinTheta = std::sqrt(std::max(0.0f, 1.0f - cosTheta * cosTheta));
			Vec3 lobe = Vec3(std::cos(phi) * sinTheta, std::sin(phi) * sinTheta, cosTheta);
			wiLocal = shadingData.frame.toWorld(wiLocal);

		}

		if (wiLocal.z <= 0.0f)
		{
			pdf = 0.0f;
			reflectedColour = Colour(0.0f, 0.0f, 0.0f);
			return Vec3(0.0f, 0.0f, 0.0f);
		}

		Vec3 wi = shadingData.frame.toWorld(wiLocal);

		// Fresnel
		float cosTheta = std::abs(woLocal.z);
		float F = ShadingHelper::fresnelDielectric(cosTheta, intIOR, extIOR);

		// 颜色和 PDF
		Colour diffuse = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		float phongNorm = (exponent + 2.0f) / (2.0f * M_PI);
		Vec3 reflectDir = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
		float dotRH = std::max(0.0f, Dot(reflectDir, wiLocal));
		Colour glossy = albedo->sample(shadingData.tu, shadingData.tv) * phongNorm * std::pow(dotRH, exponent);

		reflectedColour = diffuse*(1.0f - F) * 2.0f + glossy *F * 2.0f;
		pdf = 0.5f * SamplingDistributions::cosineHemispherePDF(wiLocal) + 0.5f * phongNorm * std::pow(dotRH, exponent);

		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Plastic evaluation code
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);

		if (wiLocal.z <= 0.0f || woLocal.z <= 0.0f)
			return Colour(0.0f, 0.0f, 0.0f);

		// Fresnel
		float cosTheta = std::abs(woLocal.z);
		float F = ShadingHelper::fresnelDielectric(cosTheta, intIOR, extIOR);

		// Diffuse term
		Colour diffuse = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;

		// Phong specular
		Vec3 reflectDir = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
		float exponent = alphaToPhongExponent();
		float dotRH = std::max(0.0f, Dot(reflectDir, wiLocal));
		float phongNorm = (exponent + 2.0f) / (2.0f * M_PI);
		Colour glossy = albedo->sample(shadingData.tu, shadingData.tv) * phongNorm * std::pow(dotRH, exponent);

		return diffuse * (1.0f - F) + glossy * F;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Plastic PDF
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		if (wiLocal.z <= 0.0f || woLocal.z <= 0.0f)
			return 0.0f;

		Vec3 reflectDir = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
		float exponent = alphaToPhongExponent();
		float dotRH = std::max(0.0f, Dot(reflectDir, wiLocal));
		float phongPDF = (exponent + 1.0f) / (2.0f * M_PI) * std::pow(dotRH, exponent);

		return 0.5f * SamplingDistributions::cosineHemispherePDF(wiLocal) + 0.5f * phongPDF;
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