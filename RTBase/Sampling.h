#pragma once

#include "Core.h"
#include <random>
#include <algorithm>

class Sampler
{
public:
	virtual float next() = 0;
};

class MTRandom : public Sampler
{
public:
	std::mt19937 generator;
	std::uniform_real_distribution<float> dist;
	MTRandom(unsigned int seed = 1) : dist(0.0f, 1.0f)
	{
		generator.seed(seed);
	}
	float next()
	{
		return dist(generator);
	}
};

// Note all of these distributions assume z-up coordinate system
class SamplingDistributions
{
public:
	static Vec3 uniformSampleHemisphere(float r1, float r2)
	{
		// Add code here
		float phi = 2.0f * M_PI * r1;
		float cosTheta = r2;
		float sinTheta = sqrt(1 - cosTheta * cosTheta);
		return Vec3(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);
	}
	static float uniformHemispherePDF(const Vec3 wi)
	{
		// Add code here
		return 1.0f / (2.0f * M_PI);
	}
	static Vec3 cosineSampleHemisphere(float r1, float r2)
	{
		// Add code here
		float r = sqrt(r1);
		float theta = 2.0f * M_PI * r2;
		float x = r * cos(theta);
		float y = r * sin(theta);
		float z = sqrt(std::max(0.0f, 1.0f - x * x - y * y));
		return Vec3(x, y, z);
	}
	static float cosineHemispherePDF(const Vec3 wi)
	{
		// Add code here
		return wi.z / M_PI;
	}
	static Vec3 uniformSampleSphere(float r1, float r2)
	{
		// Add code here
		float phi = 2.0f * M_PI * r1;
		float cosTheta = 1.0f - 2.0f * r2;
		float sinTheta = sqrt(1.0f - cosTheta * cosTheta);
		return Vec3(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);
	}
	static float uniformSpherePDF(const Vec3& wi)
	{
		// Add code here
		return 1.0f / (4.0f * M_PI);
	}
};