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
	static bool refract(const Vec3& wo, float eta, Vec3& wi)
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
	//这是 Smith 的 λ 函数，用于简化遮挡因子的计算,.就是给Gggx提供一个值
	//wi要是单位向量
	static float lambdaGGX(Vec3 wi, float alpha)
	{
		// Add code here
		float cosTheta = wi.z;
		float sinTheta = std::sqrt(std::max(0.0f, 1.0f - cosTheta * cosTheta));
		float tanTheta = sinTheta / std::max(1e-4f, cosTheta);
		float a = alpha * tanTheta;

		return (-1.0f + std::sqrt(1.0f + a * a)) * 0.5f;
	}
	//GGX版本的遮蔽函数，计算微表面的G项，调用了lambdaggx方法
	static float Gggx(Vec3 wi, Vec3 wo, float alpha)
	{
		// Add code here
		float lambdaI = lambdaGGX(wi, alpha);
		float lambdaO = lambdaGGX(wo, alpha);
		return 1.0f / ((1.0f + lambdaI) * (1.0f + lambdaO));
	}
	//GGX 的法线分布函数，用于控制微表面朝向的统计密度,计算微表面的D值
	//h要是单位向量
	static float Dggx(Vec3 h, float alpha)
	{
		// Add code here
		float cosThetaH = std::max(0.0f, std::min(1.0f, h.z));
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
		if (woLocal.z < 1e-6f)
		{
			pdf = 0.0f;
			reflectedColour = Colour(0.0f, 0.0f, 0.0f);
			return Vec3(0.0f, 0.0f, 0.0f);
		}

		// 采样半角向量（微表面法线）
		float u1 = sampler->next();
		float u2 = sampler->next();
		float cosTheta = std::sqrt((1.0f - u1) / (u1 * (alpha * alpha - 1.0f) + 1.0f));
		float sinTheta = std::sqrt(fabs(1.0f - cosTheta * cosTheta));
		float phi = 2.0f * M_PI * u2;
		Vec3 h = Vec3(std::cos(phi) * sinTheta, std::sin(phi) * sinTheta, cosTheta).normalize();
		if (h.z * woLocal.z < 0.f)
			h = -h;
		Vec3 wiLocal = h * 2 * Dot(woLocal, h) - woLocal;	
		if (wiLocal.z < 1e-6f)
		{
			pdf = 0.0f;
			reflectedColour = Colour(0.0f, 0.0f, 0.0f);
			return Vec3(0.0f, 0.0f, 0.0f);
		}
		// D,G,F
		float D = ShadingHelper::Dggx(h, alpha);
		float G = ShadingHelper::Gggx(wiLocal, woLocal, alpha);
		float cosThetaH = fabs(Dot(wiLocal, h));
		Colour F = ShadingHelper::fresnelConductor(cosThetaH, eta, k);

		//光照贡献
		Colour base = albedo->sample(shadingData.tu, shadingData.tv);
		float denom = fabs(wiLocal.z * woLocal.z);
		reflectedColour = F * base * D * G /(4 * denom);

		// PDF 计算
		float temp = fabs(Dot(wiLocal, h));
		pdf = (D * fabs(h.z)) /(4 * temp);
		//float cosWO_WH = fabs(Dot(wiLocal, h));
		//pdf = (D * fabs(h.z)) / (4.f * cosWO_WH);
		//if (pdf < 1e-8f) {
		//	std::cout << "[DEBUG] pdf too small: " << pdf << "\n";
		//}

		return shadingData.frame.toWorld(wiLocal);
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Conductor evaluation code
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		if (woLocal.z <= 0.0f || wiLocal.z <= 0.0f)
			return Colour(0.0f, 0.0f, 0.0f);
		//半角向量
		Vec3 h = (wiLocal + woLocal).normalize();

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
		// Step 1️⃣: 变换到局部空间（切空间）
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		float cosThetaO = woLocal.z;
		float cosThetaI = wiLocal.z;
		// Step 2️⃣: 上半球方向合法性检查
		if (cosThetaO <= 0.0f || cosThetaI <= 0.0f)
			return 0.0f;
		// Step 3️⃣: 构造半角向量（微表面法线）
		Vec3 h = (wiLocal + woLocal).normalize();
		// 保证 h 指向上半球（与 GGX 分布方向一致）
		if (h.z < 0.0f)
			h = -h;
		// Step 4️⃣: dot(wo, h) 是变换微表面采样到方向采样的重要 Jacobian
		float dotWoH = Dot(woLocal, h);
		if (dotWoH <= 1e-6f)
			return 0.0f;
		// Step 5️⃣: 查 GGX 分布 D(h)
		float D = ShadingHelper::Dggx(h, alpha);
		// Step 6️⃣: PDF 计算（PPT 给出）
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
//理想玻璃（只考虑两个方向上的反射和折射）
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
		float fresnel = ShadingHelper::fresnelSchlick(cosThetaI, etaI, etaT);
		//用来判断发生哪个事件
		bool reflectEvent = sampler->next() < fresnel;

		Vec3 wiLocal;
		if (reflectEvent)
		{
			// 反射
			wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * fresnel;
			pdf = fresnel;
		}
		else
		{
			//折射
			bool temp = ShadingHelper::refract(woLocal, eta, wiLocal);
			if (!temp)
			{
				// 如果是全反射
				wiLocal = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
				reflectedColour = albedo->sample(shadingData.tu, shadingData.tv);
				pdf = 1.0f;
				return shadingData.frame.toWorld(wiLocal);
			}

			//
			float scale = (eta * eta);
			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * scale * (1.0f - fresnel);
			pdf = 1.0f - fresnel;
		}

		return shadingData.frame.toWorld(wiLocal);
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Glass evaluation code
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);

		bool entering = woLocal.z > 0.0f;
		float etaI = entering ? extIOR : intIOR;
		float etaT = entering ? intIOR : extIOR;
		float eta = etaI / etaT;

		float cosThetaI = std::abs(woLocal.z);
		float fresnel = ShadingHelper::fresnelSchlick(cosThetaI, etaI, etaT);

		// === 检查是否是镜面反射方向 ===
		Vec3 expectedReflect = Vec3(-woLocal.x, -woLocal.y, woLocal.z);
		if (Dot(expectedReflect, wiLocal)>0.999f) {
			return albedo->sample(shadingData.tu, shadingData.tv) * fresnel / std::abs(expectedReflect.z);
		}

		// === 检查是否是折射方向 ===
		Vec3 expectedRefract;
		if (ShadingHelper::refract(woLocal, eta, expectedRefract)) {
			if (Dot(expectedRefract, wiLocal) > 0.999f) {
				float scale = eta * eta;
				return albedo->sample(shadingData.tu, shadingData.tv) * scale * (1.0f - fresnel) / std::abs(expectedRefract.z);
			}
		}

		// 不是镜面反射或折射方向，返回 0
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
//微表面版本的glass，不单是一个方向上的
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
	Vec3 sampleGGXNormal(Sampler* sampler) const
	{
		float u1 = sampler->next();
		float u2 = sampler->next();
		// 采用常用采样方法：tan²θ = α² u1/(1-u1)
		float tan2Theta = alpha * alpha * u1 / (1.0f - u1);
		float cosTheta = 1.0f / std::sqrt(1.0f + tan2Theta);
		float sinTheta = std::sqrt(std::max(0.0f, 1.0f - cosTheta * cosTheta));
		float phi = 2.0f * M_PI * u2;
		Vec3 m(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
		// 确保 m 与表面法线（局部空间中 n=(0,0,1)）在同一侧
		if (m.z < 0.0f) m = -m;
		return m;
	}

	// 采样函数：同时采样反射和透射分支
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		Vec3 wo = shadingData.wo; // 出射方向（局部空间）
		// 从 GGX 分布采样半角向量 m
		Vec3 m = sampleGGXNormal(sampler);
		// 计算 wo 与 m 的点积
		float cosThetaO = Dot(wo, m);
		// 计算菲涅尔反射率（注意取绝对值）
		float F = ShadingHelper::fresnelDielectric(fabs(cosThetaO), intIOR, extIOR);

		// 随机决定走反射还是透射分支
		float r = sampler->next();
		if (r < F)
		{
			// --- 反射分支 ---
			Vec3 wi = -wo + m * 2.0f * cosThetaO;
			// 计算 PDF（注意 m.z 为 |m·n|，n = (0,0,1)）
			float pdf_m = ShadingHelper::Dggx(m, alpha) * fabs(m.z);
			pdf = F * pdf_m / (4.0f * fabs(Dot(wo, m)));
			// 计算反射 BSDF 值
			float D = ShadingHelper::Dggx(m, alpha);
			float G = ShadingHelper::Gggx(wo, wi, alpha);
			float brdf = F * D * G / (4.0f * fabs(wo.z) * fabs(wi.z));
			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * brdf;
			return wi;
		}
		else
		{
			// --- 透射分支 ---
			// 根据 wo 的正负确定相对折射率
			float eta = (wo.z > 0.0f) ? (extIOR / intIOR) : (intIOR / extIOR);
			// 计算 Snell 定律所需量
			float sin2ThetaO = std::max(0.0f, 1.0f - cosThetaO * cosThetaO);
			float sin2ThetaI = eta * eta * sin2ThetaO;
			// 如果发生全内反射，退化为反射分支
			if (sin2ThetaI >= 1.0f)
			{
				Vec3 wi = -wo + m * 2.0f * cosThetaO;
				float pdf_m = ShadingHelper::Dggx(m, alpha) * fabs(m.z);
				pdf = pdf_m / (4.0f * fabs(Dot(wo, m))); // 此时 F=1
				float D = ShadingHelper::Dggx(m, alpha);
				float G = ShadingHelper::Gggx(wo, wi, alpha);
				float brdf = D * G / (4.0f * fabs(wo.z) * fabs(wi.z));
				reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * brdf;
				return wi;
			}
			float cosThetaI = std::sqrt(1.0f - sin2ThetaI);
			// 根据 wo.z 的正负确定符号，确保透射方向位于另一侧
			float sign = (wo.z > 0.0f) ? 1.0f : -1.0f;
			Vec3 wi = -wo * eta + m * (eta * cosThetaO - sign * cosThetaI);
			// 计算透射分支 PDF
			float wi_dot_m = Dot(wi, m);
			float denom = eta * cosThetaO + wi_dot_m;
			float pdf_m = ShadingHelper::Dggx(m, alpha) * fabs(m.z);
			pdf = (1.0f - F) * pdf_m * (eta * eta * fabs(wi_dot_m)) / (denom * denom);
			// 计算透射 BSDF 值（参见 Walter 的模型）
			float D = ShadingHelper::Dggx(m, alpha);
			float G = ShadingHelper::Gggx(wo, wi, alpha);
			float brdf = ((1.0f - F) * D * G * eta * eta * fabs(cosThetaO) * fabs(wi_dot_m)) /
				(fabs(wo.z) * fabs(wi.z) * (denom * denom));
			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * brdf;
			return wi;
		}
	}

	// BSDF 评价函数：给定 wo 和 wi，根据反射或透射分别计算
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 wo = shadingData.wo;
		// 判断是反射还是透射：反射时 wo 与 wi 同侧，透射则异侧
		if (wo.z * wi.z > 0.0f) {
			// --- 反射分支 ---
			Vec3 m = (wo + wi).normalize();
			float F = ShadingHelper::fresnelDielectric(fabs(Dot(wo, m)), intIOR, extIOR);
			float D = ShadingHelper::Dggx(m, alpha);
			float G = ShadingHelper::Gggx(wo, wi, alpha);
			float brdf = F * D * G / (4.0f * fabs(wo.z) * fabs(wi.z));
			return albedo->sample(shadingData.tu, shadingData.tv) * brdf;
		}
		else {
			// --- 透射分支 ---
			float eta = (wo.z > 0.0f) ? (extIOR / intIOR) : (intIOR / extIOR);
			// 透射时半角向量按如下方式计算：
			Vec3 m = (wo + wi * eta).normalize();
			float F = ShadingHelper::fresnelDielectric(fabs(Dot(wo, m)), intIOR, extIOR);
			float D = ShadingHelper::Dggx(m, alpha);
			float G = ShadingHelper::Gggx(wo, wi, alpha);
			float wi_dot_m = Dot(wi, m);
			float denom = eta * Dot(wo, m) + wi_dot_m;
			float brdf = ((1.0f - F) * D * G * eta * eta * fabs(Dot(wo, m)) * fabs(wi_dot_m)) /
				(fabs(wo.z) * fabs(wi.z) * (denom * denom));
			return albedo->sample(shadingData.tu, shadingData.tv) * brdf;
		}
	}

	// PDF 计算，同样分为反射与透射两支
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 wo = shadingData.wo;
		if (wo.z * wi.z > 0.0f) {
			// --- 反射分支 ---
			Vec3 m = (wo + wi).normalize();
			float F = ShadingHelper::fresnelDielectric(fabs(Dot(wo, m)), intIOR, extIOR);
			float pdf_m = ShadingHelper::Dggx(m, alpha) * fabs(m.z);
			return F * pdf_m / (4.0f * fabs(Dot(wo, m)));
		}
		else {
			// --- 透射分支 ---
			float eta = (wo.z > 0.0f) ? (extIOR / intIOR) : (intIOR / extIOR);
			Vec3 m = (wo + wi * eta).normalize();
			float F = ShadingHelper::fresnelDielectric(fabs(Dot(wo, m)), intIOR, extIOR);
			float pdf_m = ShadingHelper::Dggx(m, alpha) * fabs(m.z);
			float wi_dot_m = Dot(wi, m);
			float denom = eta * Dot(wo, m) + wi_dot_m;
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
//漫反射升级版
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
//塑料材质（高光+漫反射）
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
		//// Replace this with Plastic sampling code
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		bool entering = woLocal.z > 0.0f;
		float etaI = entering ? extIOR : intIOR;
		float etaT = entering ? intIOR : extIOR;
		float eta = etaI / etaT;

		float cosThetaI = std::abs(woLocal.z);
		float fresnel = ShadingHelper::fresnelSchlick(cosThetaI, etaI, etaT);
		bool reflectEvent = sampler->next() < fresnel;
		Vec3 wiLocal;
		if (reflectEvent) {
			//glossy
			float exponent = alphaToPhongExponent();
			// 采样两个随机数
			float xi1 = sampler->next();
			float xi2 = sampler->next();
			// 根据公式采样 theta 和 phi（z-up lobe）
			float theta = acos(pow(xi1, 1.0f / (exponent + 1.0f)));
			float phi = 2.0f * M_PI * xi2;

			// lobe 中的方向（局部空间 z-up）
			float sinTheta = sin(theta);
			Vec3 lobeDir(
				sinTheta * cos(phi),
				sinTheta * sin(phi),
				cos(theta)
			);
			// 计算镜面反射方向 wr
			Vec3 wr = Vec3(-woLocal.x, -woLocal.y, woLocal.z);

			// 创建以 wr 为主轴的局部坐标系
			Frame glossyFrame;
			glossyFrame.fromVector(wr);

			// 将 z-up lobe 方向旋转到 wr 方向上
			wiLocal = glossyFrame.toWorld(lobeDir);

			// 反射颜色
			float dotWRWi = std::max(0.0f, Dot(wr, wiLocal));
			pdf = fresnel * (exponent + 1.0f) / (2.0f * M_PI) * pow(dotWRWi, exponent);
			float phongEval = (exponent + 2.0f) / (2.0f * M_PI) * pow(dotWRWi, exponent);
			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * phongEval;
		}
		else {
			//漫反射
			wiLocal = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
			pdf = (1- fresnel)*wiLocal.z / M_PI;
			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / (M_PI);
		}
		return shadingData.frame.toWorld(wiLocal);
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Plastic evaluation code
			// === Step 1: 把方向从世界空间变到局部空间 ===
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		return Colour(0.0f, 0.0f, 0.0f);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Replace this with Plastic PDF
		return 0.f;
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