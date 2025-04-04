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
	//compute fresnel using physical formula
	static float fresnelDielectric(float cosTheta, float etaI, float etaT)
	{
		float cos = std::abs(cosTheta);
		float sinSQ = (etaI / etaT) * (etaI / etaT) * (1.0f - cos * cos);
		if (sinSQ >= 1.0f) return 1.0f;
		cos = std::sqrt(1.0f - sinSQ);

		float Rs = ((etaT * cos) - (etaI * cos)) /
			((etaT * cos) + (etaI * cos));
		float Rp = ((etaI * cos) - (etaT * cos)) /
			((etaI * cos) + (etaT * cos));

		return (Rs * Rs + Rp * Rp) * 0.5f;
	}
	//compute fresnel using schlick approximation
	static float fresnelSchlick(float cosTheta, float etaI, float etaT)
	{
		float f0 = (etaI - etaT) / (etaI + etaT);
		f0 = f0 * f0;

		float temp = 1.0f - cosTheta;
		return f0 + (1.0f - f0) * std::pow(temp, 5.0f);
	}
	//compute fresnel for conductor
	static Colour fresnelConductor(float cosTheta, Colour ior, Colour k)
	{
		Colour cos2 = Colour(cosTheta * cosTheta, cosTheta * cosTheta, cosTheta * cosTheta);
		Colour sin2 = Colour(1.0f, 1.0f, 1.0f) - cos2;

		Colour ior2 = ior * ior;
		Colour k2 = k * k;
		Colour t0 = ior2 + k2;
		Colour twoEtaCos = ior * (2.0f * cosTheta);

		Colour r1 = ((t0 * cos2 - twoEtaCos + sin2) /
			(t0 * cos2 + twoEtaCos + sin2));

		Colour r2 = ((t0 - twoEtaCos + cos2) /
			(t0 + twoEtaCos + cos2));

		return (r1 + r2) * 0.5f;
	}
	
	//compute the lambda for GGX
	static float lambdaGGX(Vec3 wi, float alpha)
	{
		float cos = wi.z;
		float sin = std::sqrt(std::max(0.0f, 1.0f - cos * cos));
		float tan = sin / std::max(1e-4f, cos);
		float a = alpha * tan;

		return (-1.0f + std::sqrt(1.0f + a * a)) * 0.5f;
	}
	//compute the G
	static float Gggx(Vec3 wi, Vec3 wo, float alpha)
	{
		float lambdaI = lambdaGGX(wi, alpha);
		float lambdaO = lambdaGGX(wo, alpha);
		return 1.0f / ((1.0f + lambdaI) * (1.0f + lambdaO));
	}
	//compute the D
	static float Dggx(Vec3 h, float alpha)
	{
		float cos = std::max(0.0f, std::min(1.0f, h.z));
		if (cos <= 0.0f) return 0.0f;
		float alpha2 = alpha * alpha;
		float cos2 = cos * cos;

		float denom = cos2 * (alpha2 - 1.0f) + 1.0f;
		denom = M_PI * denom * denom;

		return alpha2 / denom;
	}
};

class BSDF
{
public:
	Colour emission;
	virtual Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf, Colour& a) = 0;
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
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf, Colour& a)
	{
		// Add correct sampling code here
		Vec3 wi = SamplingDistributions::cosineSampleHemisphere(sampler->next(),sampler->next());
		pdf = wi.z/M_PI;
		a = albedo->sample(shadingData.tu, shadingData.tv);
		reflectedColour = a/ M_PI;
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
//Perfect Specular Reflection(Only reflection)
class MirrorBSDF : public BSDF
{
public:
	Texture* albedo;
	MirrorBSDF() = default;
	MirrorBSDF(Texture* _albedo)
	{
		albedo = _albedo;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf, Colour& a)
	{
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);

		//reflection direction
		Vec3 wi = Vec3(-wo.x, -wo.y, wo.z);

		//relection colour
		//didn't divide it by the cos term 
		//after dividing, the specular sphere in the material ball scene looked strange
		//not sure if this is due to the cos division already being applied in the UV texture.
		a= albedo->sample(shadingData.tu, shadingData.tv);
		reflectedColour = a;

		wi = shadingData.frame.toWorld(wi);
		pdf = 1.0f;
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		//specular will be handled specially, so set it to 0
		return Colour(0.0f, 0.0f, 0.0f);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// similar to evaluate
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
//Conductor(Only reflection with microfacet)
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
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf, Colour& a)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);

		//compute the direction
		float u1 = sampler->next();
		float u2 = sampler->next();
		float cos = std::sqrt((1.0f - u1) / (u1 * (alpha * alpha - 1.0f) + 1.0f));
		float sin = std::sqrt(fabs(1.0f - cos * cos));
		float phi = 2.0f * M_PI * u2;
		Vec3 h = Vec3(std::cos(phi) * sin, std::sin(phi) * sin, cos).normalize();
		if (h.z * woLocal.z < 0.f)
			h = -h;
		Vec3 wiLocal = h * 2 * Dot(woLocal, h) - woLocal;	

		//compute the D,G,F
		float D = ShadingHelper::Dggx(h, alpha);
		float G = ShadingHelper::Gggx(wiLocal, woLocal, alpha);
		float cosH = fabs(Dot(wiLocal, h));
		Colour F = ShadingHelper::fresnelConductor(cosH, eta, k);

		//comnpute the colour
		a = albedo->sample(shadingData.tu, shadingData.tv);
		float denom = fabs(wiLocal.z * woLocal.z);
		reflectedColour = F * a * D * G /(4 * denom);

		//comnpute the pdf
		float temp = fabs(Dot(wiLocal, h));
		pdf = (D * fabs(h.z)) /(4 * temp);

		return shadingData.frame.toWorld(wiLocal);
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		
		//compute fresnel
		Vec3 h = (wiLocal + woLocal).normalize();
		float cosThetaH = std::max(0.0f, Dot(wiLocal, h));
		Colour F = ShadingHelper::fresnelConductor(cosThetaH, eta, k);

		//compute D、G
		float D = ShadingHelper::Dggx(h, alpha);
		float G = ShadingHelper::Gggx(wiLocal, woLocal, alpha);

		//compute colour
		Colour base = albedo->sample(shadingData.tu, shadingData.tv);
		float denom = 4.0f * std::max(1e-4f, wiLocal.z * woLocal.z);
		return base * F * (D * G / denom);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);

		Vec3 h = (wiLocal + woLocal).normalize();
		float dotWoH = Dot(woLocal, h);
		
		//compute D
		float D = ShadingHelper::Dggx(h, alpha);
		//compute pdf
		float pdf = (D * fabs(h.z)) / (4.0f * fabs(dotWoH));

		return pdf;
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
//Glass（Only reflections and refractions in two directions）
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
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf, Colour& a)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		
		//adjust the incident direction
		bool entering = woLocal.z > 0.0f;
		float etaI = entering ? extIOR : intIOR;
		float etaT = entering ? intIOR : extIOR;
		float eta = etaI / etaT;

		//compute fresnel
		float cosI = std::abs(woLocal.z);
		float fresnel = ShadingHelper::fresnelDielectric(cosI, etaI, etaT);


		Vec3 wiLocal;
		if (sampler->next() < fresnel) {
			//path 1: reflection
			wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
			a = albedo->sample(shadingData.tu, shadingData.tv);
			reflectedColour = a * fresnel / std::abs(woLocal.z);
			pdf = fresnel;
		}
		else {
			//path 2: refraction
			Vec3 n = (woLocal.z > 0.0f) ? Vec3(0, 0, 1) : Vec3(0, 0, -1);
			float cosTT = Dot(woLocal, n);
			float k=1.0f- eta * eta * (1.0f - cosTT * cosTT);

			//check k
			if (k < 0.0f) {
				// use reflection
				wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
				a = albedo->sample(shadingData.tu, shadingData.tv);
				reflectedColour = a * fresnel / std::abs(woLocal.z);
				pdf = 1.0f;
			}
			else {
				//use refraction
				wiLocal = Vec3(-eta * woLocal.x,-eta * woLocal.y,(entering ? -1 : 1)* std::sqrt(k));
				a = albedo->sample(shadingData.tu, shadingData.tv);
				float scale = (eta * eta);
				reflectedColour = a * scale * (1.0f - fresnel) / std::abs(wiLocal.z);
				pdf = 1.0f - fresnel;
			}
		}
		return shadingData.frame.toWorld(wiLocal);
	}

	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// similar to mirror
		return Colour(0.0f, 0.0f, 0.0f);
	}

	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// similar to mirror
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
//dielectric(reglections and refractions with microfacet)
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
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf, Colour& a)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);

		//compute half-vector
		float u1 = sampler->next();
		float u2 = sampler->next();
		float cos = std::sqrt((1.0f - u1) / (u1 * (alpha * alpha - 1.0f) + 1.0f));
		float sin = std::sqrt(fabs(1.0f - cos * cos));
		float phi = 2.0f * M_PI * u2;
		Vec3 h = Vec3(std::cos(phi) * sin, std::sin(phi) * sin, cos).normalize();
		if (h.z * woLocal.z < 0.f)
			h = -h;

		//compute fresnel
		float cosH = fabs(Dot(woLocal, h));
		float fresnel = ShadingHelper::fresnelSchlick(cosH, intIOR, extIOR);

		if (sampler->next() < fresnel)
		{
			//path 1: reflection
			//direction
			Vec3 wi = -woLocal + h * 2.0f * cosH;
			//pdf
			float pdf_m = ShadingHelper::Dggx(h, alpha) * fabs(h.z);
			pdf = fresnel * pdf_m / (4.0f * cosH);
			//colour
			float D = ShadingHelper::Dggx(h, alpha);
			float G = ShadingHelper::Gggx(woLocal, wi, alpha);
			float brdf = fresnel * D * G / (4.0f * fabs(woLocal.z) * fabs(wi.z));
			a = albedo->sample(shadingData.tu, shadingData.tv);
			reflectedColour = a * brdf;
			return wi;
		}
		else
		{
			//path 2: refraction
			float eta = (woLocal.z > 0.0f) ? (extIOR / intIOR) : (intIOR / extIOR);
			float sin2O = std::max(0.0f, 1.0f - cosH * cosH);
			float sin2I = eta * eta * sin2O;

			if (sin2I >= 1.0f)
			{
				Vec3 wi = -woLocal + h * 2.0f * cosH;
				float pdf_m = ShadingHelper::Dggx(h, alpha) * fabs(h.z);
				pdf = pdf_m / (4.0f * fabs(Dot(woLocal, h)));
				float D = ShadingHelper::Dggx(h, alpha);
				float G = ShadingHelper::Gggx(woLocal, wi, alpha);
				float brdf = fresnel * D * G / (4.0f * fabs(woLocal.z) * fabs(wi.z));
				a = albedo->sample(shadingData.tu, shadingData.tv);
				reflectedColour = a * brdf;
				return wi;
			}
			float cosI = std::sqrt(1.0f - sin2I);
			float sign = (woLocal.z > 0.0f) ? 1.0f : -1.0f;
			Vec3 wi = -woLocal * eta + h * (eta * cosH - sign * cosI);
			float wi_dot_m = Dot(wi, h);
			float denom = eta * cosH + wi_dot_m;
			float pdf_m = ShadingHelper::Dggx(h, alpha) * fabs(h.z);
			pdf = (1.0f - fresnel) * pdf_m * (eta * eta * fabs(wi_dot_m)) / (denom * denom);
			float D = ShadingHelper::Dggx(h, alpha);
			float G = ShadingHelper::Gggx(woLocal, wi, alpha);
			float brdf = ((1.0f - fresnel) * D * G * eta * eta * fabs(wi_dot_m)) /
				(fabs(woLocal.z) * fabs(wi.z) * (denom * denom));
			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * brdf;
			return wi;
		}
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		
		if (woLocal.z * wiLocal.z > 0.0f) {
			//path 1: reflection
			Vec3 h = (woLocal + wiLocal).normalize();
			float F = ShadingHelper::fresnelSchlick(fabs(Dot(woLocal, h)), intIOR, extIOR);
			float D = ShadingHelper::Dggx(h, alpha);
			float G = ShadingHelper::Gggx(woLocal, wiLocal, alpha);
			float brdf = F * D * G / (4.0f * fabs(woLocal.z) * fabs(wiLocal.z));
			return albedo->sample(shadingData.tu, shadingData.tv) * brdf;
		}
		else {
			//path 2: refraction
			float eta = (woLocal.z > 0.0f) ? (extIOR / intIOR) : (intIOR / extIOR);
			Vec3 h = (woLocal + wiLocal * eta).normalize();
			float F = ShadingHelper::fresnelSchlick(fabs(Dot(woLocal, h)), intIOR, extIOR);
			float D = ShadingHelper::Dggx(h, alpha);
			float G = ShadingHelper::Gggx(woLocal, wiLocal, alpha);
			float wi_dot_m = Dot(wiLocal, h);
			float denom = eta * Dot(woLocal, h) + wi_dot_m;
			float brdf = ((1.0f - F) * D * G * eta * eta * fabs(wi_dot_m)) /
				(fabs(woLocal.z) * fabs(wiLocal.z) * (denom * denom));
			return albedo->sample(shadingData.tu, shadingData.tv) * brdf;
		}
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);

		if (woLocal.z * wiLocal.z > 0.0f) {
			//path 1: reflection
			Vec3 h = (woLocal + wiLocal).normalize();
			float F = ShadingHelper::fresnelSchlick(fabs(Dot(woLocal, h)), intIOR, extIOR);
			float pdf_m = ShadingHelper::Dggx(h, alpha) * fabs(h.z);
			return F * pdf_m / (4.0f * fabs(Dot(woLocal, h)));
		}
		else {
			//path 2: refraction
			float eta = (woLocal.z > 0.0f) ? (extIOR / intIOR) : (intIOR / extIOR);
			Vec3 h = (woLocal + wiLocal * eta).normalize();
			float F = ShadingHelper::fresnelSchlick(fabs(Dot(woLocal, h)), intIOR, extIOR);
			float pdf_m = ShadingHelper::Dggx(h, alpha) * fabs(h.z);
			float wi_dot_m = Dot(wiLocal, h);
			float denom = eta * Dot(woLocal, h) + wi_dot_m;
			return (1.0f - F) * pdf_m * (eta * eta * fabs(wi_dot_m)) / (denom * denom);
		}
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
//Not implemented
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
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf, Colour& a)
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
//Plastic (Phong glossy and diffuse)
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
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf, Colour& a)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);

		//adjust the incident direction
		bool entering = woLocal.z > 0.0f;
		float etaI = entering ? extIOR : intIOR;
		float etaT = entering ? intIOR : extIOR;
		float eta = etaI / etaT;

		//compute fresnel
		float cosI = std::abs(woLocal.z);
		float fresnel = ShadingHelper::fresnelDielectric(cosI, etaI, etaT);

		Vec3 wiLocal;
		if (sampler->next() < fresnel) {
			//path 1: Phong glossy
			//compute direction
			float exponent = alphaToPhongExponent();
			float xi1 = sampler->next();
			float xi2 = sampler->next();
			float theta = acos(pow(xi1, 1.0f / (exponent + 1.0f)));
			float phi = 2.0f * M_PI * xi2;
			float sinT = sin(theta);
			Vec3 dir(
				sinT * cos(phi),
				sinT * sin(phi),
				cos(theta)
			);
			Vec3 wr = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
			Frame glossyFrame;
			glossyFrame.fromVector(wr);
			wiLocal = glossyFrame.toWorld(dir);

			// compute colour
			float dotWRWi = std::max(0.0f, Dot(wr, wiLocal));
			pdf = fresnel * (exponent + 1.0f) / (2.0f * M_PI) * pow(dotWRWi, exponent);
			float phongEval = fresnel * (exponent + 2.0f) / (2.0f * M_PI) * pow(dotWRWi, exponent);
			a = albedo->sample(shadingData.tu, shadingData.tv);
			reflectedColour = a * phongEval;
		}
		else {
			//path 1: diffuse
			wiLocal = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
			pdf = (1- fresnel)* wiLocal.z / M_PI;
			a = albedo->sample(shadingData.tu, shadingData.tv);
			reflectedColour = a * (1 - fresnel) / (M_PI);
		}
		return shadingData.frame.toWorld(wiLocal);
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);

		//adjust the incident direction
		bool entering = woLocal.z > 0.0f;
		float etaI = entering ? extIOR : intIOR;
		float etaT = entering ? intIOR : extIOR;
		float eta = etaI / etaT;

		//compute fresnel
		float cosI = std::abs(woLocal.z);
		float fresnel = ShadingHelper::fresnelDielectric(cosI, etaI, etaT);

		//diffuse part
		Colour diffuse = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;

		//glossy part
		Vec3 wr = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
		float dotWRWi = std::max(0.0f, Dot(wr, wiLocal));
		float exponent = alphaToPhongExponent();
		Colour glossy = albedo->sample(shadingData.tu, shadingData.tv)
			*((exponent + 2.0f) / (2.0f * M_PI)) *pow(dotWRWi, exponent);

		return glossy * fresnel + diffuse * (1.0f - fresnel);
	}

	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);

		//adjust the incident direction
		bool entering = woLocal.z > 0.0f;
		float etaI = entering ? extIOR : intIOR;
		float etaT = entering ? intIOR : extIOR;
		float eta = etaI / etaT;

		//compute fresnel
		float cosI = std::abs(woLocal.z);
		float fresnel = ShadingHelper::fresnelDielectric(cosI, etaI, etaT);

		//diffuse part
		float diffuse = (wiLocal.z > 0.0f) ? (wiLocal.z / M_PI) : 0.0f;

		//glossy part
		Vec3 wr = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
		float exponent = alphaToPhongExponent();
		float dotWRWi = std::max(0.0f, Dot(wr, wiLocal));
		float glossy = ((exponent + 1.0f) / (2.0f * M_PI)) * pow(dotWRWi, exponent);

		return fresnel * glossy + (1.0f - fresnel) * diffuse;
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
//Not implemented
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
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf, Colour& a)
	{
		// Add code to include layered sampling
		return base->sample(shadingData, sampler, reflectedColour, pdf,a);
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