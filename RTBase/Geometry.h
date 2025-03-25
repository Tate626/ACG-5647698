#pragma once

#include "Core.h"
#include "Sampling.h"

class Ray
{
public:
	Vec3 o;
	Vec3 dir;
	Vec3 invDir;
	Ray()
	{
	}
	Ray(Vec3 _o, Vec3 _d)
	{
		init(_o, _d);
	}
	void init(Vec3 _o, Vec3 _d)
	{
		o = _o;
		dir = _d;
		invDir = Vec3(1.0f / dir.x, 1.0f / dir.y, 1.0f / dir.z);
	}
	Vec3 at(const float t) const
	{
		return (o + (dir * t));
	}
};

class Plane
{
public:
	Vec3 n;
	float d;
	void init(Vec3& _n, float _d)
	{
		n = _n;
		d = _d;
	}
	// Add code here
	bool rayIntersect(Ray& r, float& t)
	{
		float dt = Dot(n, r.dir);
		float temp = d - Dot(n, r.o);
		if (dt != 0) {
			temp = temp / dt;
			if (temp > 0) {
				return true;
			}
		}
		return false;
	}
};

#define EPSILON 0.001f

class Triangle
{
public:
	Vertex vertices[3];
	Vec3 e1; // Edge 1
	Vec3 e2; // Edge 2
	Vec3 n; // Geometric Normal
	float area; // Triangle area
	float d; // For ray triangle if needed
	unsigned int materialIndex;
	void init(Vertex v0, Vertex v1, Vertex v2, unsigned int _materialIndex)
	{
		materialIndex = _materialIndex;
		vertices[0] = v0;
		vertices[1] = v1;
		vertices[2] = v2;
		//e1 = vertices[2].p - vertices[1].p;
		//e2 = vertices[0].p - vertices[2].p;
		e1 = vertices[1].p - vertices[0].p;
		e2 = vertices[2].p - vertices[0].p;
		n = e1.cross(e2).normalize();
		area = e1.cross(e2).length() * 0.5f;
		d = Dot(n, vertices[0].p);
	}
	Vec3 centre() const
	{
		return (vertices[0].p + vertices[1].p + vertices[2].p) / 3.0f;
	}
	// Add code here
	//m法，注意需要修改init
	bool rayIntersect(const Ray& r, float& t, float& u, float& v) const
	{
		Vec3 P = r.dir.cross(e2);
		float det = Dot(e1, P);

		if (fabs(det) < 1e-6) return false;

		float invdet = 1.0 / det;
		Vec3 T = r.o - vertices[0].p;
		u = Dot(T, P) * invdet;
		if (u < 0 || u > 1) return false;

		Vec3 Q = T.cross(e1);
		v = Dot(r.dir, Q) * invdet;
		if (v < 0 || u + v > 1) return false;

		t = Dot(e2, Q) * invdet;
		if (t < 1e-5) return false;
		return true;
	}
	void interpolateAttributes(const float alpha, const float beta, const float gamma, Vec3& interpolatedNormal, float& interpolatedU, float& interpolatedV) const
	{
		interpolatedNormal = vertices[0].normal * alpha + vertices[1].normal * beta + vertices[2].normal * gamma;
		interpolatedNormal = interpolatedNormal.normalize();
		interpolatedU = vertices[0].u * alpha + vertices[1].u * beta + vertices[2].u * gamma;
		interpolatedV = vertices[0].v * alpha + vertices[1].v * beta + vertices[2].v * gamma;
	}
	// Add code here
	Vec3 sample(Sampler* sampler, float& pdf)
	{
		float r1 = sampler->next();
		float r2 = sampler->next();

		// 根据 Barycentric 坐标对三角形进行采样
		float sqrtr1 = sqrt(r1);
		float dV1 = 1.0f - sqrtr1;
		float dV2 = r2 * sqrtr1;
		float dV3 = 1 - dV1 - dV2;
		Vec3 pos = vertices[0].p * dV1 + vertices[1].p * dV2 + vertices[2].p * dV3;
		return pos;
	}
	Vec3 gNormal()
	{
		return (n * (Dot(vertices[0].normal, n) > 0 ? 1.0f : -1.0f));
	}
};

class AABB
{
public:
	Vec3 max;
	Vec3 min;
	AABB()
	{
		reset();
	}
	void reset()
	{
		max = Vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
		min = Vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	}
	void extend(const Vec3 p)
	{
		max = Max(max, p);
		min = Min(min, p);
	}

	void extend(const Triangle t) {
		this->extend(t.vertices[0].p);
		this->extend(t.vertices[1].p);
		this->extend(t.vertices[2].p);
	}

	// Add code here
	bool rayAABB(const Ray& r, float& t)
	{
		Vec3 Tmin = (min - r.o) * r.invDir;
		Vec3 Tmax = (max - r.o) * r.invDir;
		if (Tmin.x > Tmax.x)
		{
			float temp = Tmax.x;
			Tmax.x = Tmin.x;
			Tmin.x = temp;
		}
		if (Tmin.y > Tmax.y)
		{
			float temp = Tmax.y;
			Tmax.y = Tmin.y;
			Tmin.y = temp;
		}
		if (Tmin.z > Tmax.z)
		{
			float temp = Tmax.z;
			Tmax.z = Tmin.z;
			Tmin.z = temp;
		}
		float tmin = std::max(Tmin.x, std::max(Tmin.y, Tmin.z));
		float tmax = std::min(Tmax.x, std::min(Tmax.y, Tmax.z));
		if (tmin > tmax) return false;
		if (tmax < 0) return false;
		t = tmin;
		return true;
	}
	// Add code here
	bool rayAABB(const Ray& r)
	{
		Vec3 Tmin = (min - r.o) * r.invDir;
		Vec3 Tmax = (max - r.o) * r.invDir;
		if (Tmin.x > Tmax.x)
		{
			float temp = Tmax.x;
			Tmax.x = Tmin.x;
			Tmin.x = temp;
		}
		if (Tmin.y > Tmax.y)
		{
			float temp = Tmax.y;
			Tmax.y = Tmin.y;
			Tmin.y = temp;
		}
		if (Tmin.z > Tmax.z)
		{
			float temp = Tmax.z;
			Tmax.z = Tmin.z;
			Tmin.z = temp;
		}
		float tmin = std::max(Tmin.x, std::max(Tmin.y, Tmin.z));
		float tmax = std::min(Tmax.x, std::min(Tmax.y, Tmax.z));
		if (tmin > tmax)return false;
		if (tmax < 0)return false;
		return true;
	}
	// Add code here
	float area()
	{
		Vec3 size = max - min;
		return ((size.x * size.y) + (size.y * size.z) + (size.x * size.z)) * 2.0f;
	}
};

class Sphere
{
public:
	Vec3 centre;
	float radius;
	void init(Vec3& _centre, float _radius)
	{
		centre = _centre;
		radius = _radius;
	}
	// Add code here
	bool rayIntersect(Ray& r, float& t)
	{
		Vec3 L = r.o - centre;
		float A = Dot(r.dir, r.dir);
		if (A == 0)return false;
		float B = 2 * Dot(L, r.dir);
		float C = Dot(L, L) - radius * radius;
		float temp = B * B - 4 * A * C;
		if (temp > 0) {
			float D = std::sqrt(temp);
			float inv = 0.5f / A;
			float t1 = (-B - D) * inv;
			float t2 = (-B + D) * inv;
			float tmin = std::min(t1, t2);
			float tmax = std::max(t1, t2);
			if (tmin >= 0) {
				t = tmin;
				return true;
			}
			else if (tmax < 0) {
				return false;
			}
			t = tmax;
			return true;
		}

		return false;
	}
};

//用于储存光线发射后与最近三角形相交得到的交点信息
struct IntersectionData
{
	unsigned int ID;//三角形编号
	float t;//交点
	float alpha;//颜色
	float beta;
	float gamma;
};

#define MAXNODE_TRIANGLES 8
#define TRAVERSE_COST 1.0f
#define TRIANGLE_COST 2.0f
#define BUILD_BINS 32

class BVHNode
{
public:
	AABB bounds;
	BVHNode* r;
	BVHNode* l;
	// This can store an offset and number of triangles in a global triangle list for example
	// But you can store this however you want!
	// unsigned int offset;
	// unsigned char num;
	BVHNode()
	{
		r = NULL;
		l = NULL;
	}
	// Note there are several options for how to implement the build method. Update this as required
	void build(std::vector<Triangle>& inputTriangles)
	{
		// Add BVH building code here
	}
	void traverse(const Ray& ray, const std::vector<Triangle>& triangles, IntersectionData& intersection)
	{
		// Add BVH Traversal code here
	}
	IntersectionData traverse(const Ray& ray, const std::vector<Triangle>& triangles)
	{
		IntersectionData intersection;
		intersection.t = FLT_MAX;
		traverse(ray, triangles, intersection);
		return intersection;
	}
	bool traverseVisible(const Ray& ray, const std::vector<Triangle>& triangles, const float maxT)
	{
		// Add visibility code here
		return true;
	}
};
