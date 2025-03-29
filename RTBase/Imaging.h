#pragma once

#include "Core.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define __STDC_LIB_EXT1__
#include "stb_image_write.h"
#include <OpenImageDenoise/oidn.hpp>

// Stop warnings about buffer overruns if size is zero. Size should never be zero and if it is the code handles it.
#pragma warning( disable : 6386)

constexpr float texelScale = 1.0f / 255.0f;

class Texture
{
public:
	Colour* texels;
	float* alpha;
	int width;
	int height;
	int channels;
	void loadDefault()
	{
		width = 1;
		height = 1;
		channels = 3;
		texels = new Colour[1];
		texels[0] = Colour(1.0f, 1.0f, 1.0f);
	}
	void load(std::string filename)
	{
		alpha = NULL;
		if (filename.find(".hdr") != std::string::npos)
		{
			float* textureData = stbi_loadf(filename.c_str(), &width, &height, &channels, 0);
			if (width == 0 || height == 0)
			{
				loadDefault();
				return;
			}
			texels = new Colour[width * height];
			for (int i = 0; i < (width * height); i++)
			{
				texels[i] = Colour(textureData[i * channels], textureData[(i * channels) + 1], textureData[(i * channels) + 2]);
			}
			stbi_image_free(textureData);
			return;
		}
		unsigned char* textureData = stbi_load(filename.c_str(), &width, &height, &channels, 0);
		if (width == 0 || height == 0)
		{
			loadDefault();
			return;
		}
		texels = new Colour[width * height];
		for (int i = 0; i < (width * height); i++)
		{
			texels[i] = Colour(textureData[i * channels] / 255.0f, textureData[(i * channels) + 1] / 255.0f, textureData[(i * channels) + 2] / 255.0f);
		}
		if (channels == 4)
		{
			alpha = new float[width * height];
			for (int i = 0; i < (width * height); i++)
			{
				alpha[i] = textureData[(i * channels) + 3] / 255.0f;
			}
		}
		stbi_image_free(textureData);
	}
	Colour sample(const float tu, const float tv) const
	{
		Colour tex;
		float u = std::max(0.0f, fabsf(tu)) * width;
		float v = std::max(0.0f, fabsf(tv)) * height;
		int x = (int)floorf(u);
		int y = (int)floorf(v);
		float frac_u = u - x;
		float frac_v = v - y;
		float w0 = (1.0f - frac_u) * (1.0f - frac_v);
		float w1 = frac_u * (1.0f - frac_v);
		float w2 = (1.0f - frac_u) * frac_v;
		float w3 = frac_u * frac_v;
		x = x % width;
		y = y % height;
		Colour s[4];
		s[0] = texels[y * width + x];
		s[1] = texels[y * width + ((x + 1) % width)];
		s[2] = texels[((y + 1) % height) * width + x];
		s[3] = texels[((y + 1) % height) * width + ((x + 1) % width)];
		tex = (s[0] * w0) + (s[1] * w1) + (s[2] * w2) + (s[3] * w3);
		return tex;
	}
	float sampleAlpha(const float tu, const float tv) const
	{
		if (alpha == NULL)
		{
			return 1.0f;
		}
		float tex;
		float u = std::max(0.0f, fabsf(tu)) * width;
		float v = std::max(0.0f, fabsf(tv)) * height;
		int x = (int)floorf(u);
		int y = (int)floorf(v);
		float frac_u = u - x;
		float frac_v = v - y;
		float w0 = (1.0f - frac_u) * (1.0f - frac_v);
		float w1 = frac_u * (1.0f - frac_v);
		float w2 = (1.0f - frac_u) * frac_v;
		float w3 = frac_u * frac_v;
		x = x % width;
		y = y % height;
		float s[4];
		s[0] = alpha[y * width + x];
		s[1] = alpha[y * width + ((x + 1) % width)];
		s[2] = alpha[((y + 1) % height) * width + x];
		s[3] = alpha[((y + 1) % height) * width + ((x + 1) % width)];
		tex = (s[0] * w0) + (s[1] * w1) + (s[2] * w2) + (s[3] * w3);
		return tex;
	}
	~Texture()
	{
		delete[] texels;
		if (alpha != NULL)
		{
			delete alpha;
		}
	}
};

class ImageFilter
{
public:
	virtual float filter(const float x, const float y) const = 0;
	virtual int size() const = 0;
};

class BoxFilter : public ImageFilter
{
public:
	float filter(float x, float y) const
	{
		if (fabsf(x) < 0.5f && fabs(y) < 0.5f)
		{
			return 1.0f;
		}
		return 0;
	}
	int size() const
	{
		return 0;
	}
};

class GaussianFilter : public ImageFilter {
public:
	float sigma; // 高斯分布的标准差
	int radius;  // 滤波器半径，即向四周扩展的像素数

	// 构造函数，可传入 sigma，默认值为 1.0f
	GaussianFilter(float sigma = 1.0f) : sigma(sigma) {
		// 使用 ceil(2*sigma) 作为半径，通常能覆盖约 95% 的能量
		radius = static_cast<int>(std::ceil(2 * sigma));
	}

	// 返回给定偏移处的高斯权重（未归一化，归一化在 splat 中完成）
	float filter(float x, float y) const override {
		return std::exp(-(x * x + y * y) / (2 * sigma * sigma));
	}

	// 返回滤波器半径。splat 中会遍历 [ -size, size ] 范围
	int size() const override {
		return radius;
	}
};

//旧版，不支持自适应采样
//class Film
//{
//public:
//	film这个数组管理着累计渲染结果，splat每帧把结果存入，tomap除以spp就是累计渲染结果
//	Colour* film;
//	unsigned int width;//图片尺寸（分辨率）
//	unsigned int height;
//	int SPP;//采样次数（每个像素）
//	ImageFilter* filter;
//
//	void splat(const float x, const float y, const Colour& L)
//		计算每束光线对周围像素的影响，传入的x，y为浮点数，代表这束光线的位置
//		使用各种filter来计算这个光线对周围像素的权重
//	{
//		float filterWeights[25]; // Storage to cache weights 
//		unsigned int indices[25]; // Store indices to minimize computations 
//		unsigned int used = 0;
//		float total = 0;
//		int size = filter->size();
//		for (int i = -size; i <= size; i++) {
//			for (int j = -size; j <= size; j++) {
//				int px = (int)x + j;
//				int py = (int)y + i;
//				if (px >= 0 && px < width && py >= 0 && py < height) {
//					indices[used] = (py * width) + px;
//					filterWeights[used] = filter->filter(j, i);
//					total += filterWeights[used];
//					used++;
//				}
//			}
//		}
//		for (int i = 0; i < used; i++) {
//			film[indices[i]] = film[indices[i]] + (L * filterWeights[i] / total);
//		}
//		 Code to splat a smaple with colour L into the image plane using an ImageFilter
//	}
//	void tonemap(int x, int y, unsigned char& r, unsigned char& g, unsigned char& b, float exposure = 1.0f)
//		普通伽马校正方法
//	{
//		Colour pixel = film[y * width + x] * (exposure / (float)SPP);
//		r = std::min(powf(std::max(pixel.r, 0.0f), 1.0f / 2.2f) * 255, 255.0f);
//		g = std::min(powf(std::max(pixel.g, 0.0f), 1.0f / 2.2f) * 255, 255.0f);
//		b = std::min(powf(std::max(pixel.b, 0.0f), 1.0f / 2.2f) * 255, 255.0f);
//		 Return a tonemapped pixel at coordinates x, y
//	}
//	 Do not change any code below this line
//	void init(int _width, int _height, ImageFilter* _filter)
//	{
//		width = _width;
//		height = _height;
//		film = new Colour[width * height];
//		clear();
//		filter = _filter;
//	}
//	void clear()
//	{
//		memset(film, 0, width * height * sizeof(Colour));
//		SPP = 0;
//	}
//	void incrementSPP()
//	{
//		SPP++;
//	}
//	void save(std::string filename)
//	{
//		Colour* hdrpixels = new Colour[width * height];
//		for (unsigned int i = 0; i < (width * height); i++)
//		{
//			hdrpixels[i] = film[i] / (float)SPP;
//		}
//		stbi_write_hdr(filename.c_str(), width, height, 3, (float*)hdrpixels);
//		delete[] hdrpixels;
//	}
//};

class Film;

void DenoiseFilm(Film* film);

class Film
{
public:
	//film这个数组管理着累计渲染结果，splat每帧把结果存入，tomap除以spp就是累计渲染结果
	Colour* film;
	unsigned int width;//图片尺寸（分辨率）
	unsigned int height;
	int* sppBuffer;//记录每像素的实际采样数
	float* weightBuffer;//优化sppbuffer
	//以下三个存储用于降噪的
	float* albedoBuffer;
	float* normalBuffer;
	float* colourBuffer;
	float* outputBuffer;
	int SPP;//循环帧数（每个像素）
	ImageFilter* filter;

	void splat(const float x, const float y, const Colour& L)
	{
		float filterWeights[25];      // 滤波器的权重
		unsigned int indices[25];     // 对应像素在 buffer 中的 index
		unsigned int used = 0;
		float total = 0;
		int size = filter->size();

		// 1. 计算滤波器覆盖的区域 & 权重
		for (int i = -size; i <= size; i++) {
			for (int j = -size; j <= size; j++) {
				int px = (int)x + j;
				int py = (int)y + i;
				if (px >= 0 && px < (int)width && py >= 0 && py < (int)height) {
					int index = py * width + px;
					indices[used] = index;
					filterWeights[used] = filter->filter(j, i);
					total += filterWeights[used];
					used++;
				}
			}
		}

		// 2. 写入 film 累计原始值、更新权重
		for (int i = 0; i < used; i++) {
			int index = indices[i];
			float normWeight = filterWeights[i] / total;

			// 👇 未平均的累计值
			film[index] = film[index] + (L * normWeight);

			// 👇 累计权重
			weightBuffer[index] += normWeight;

			sppBuffer[index]++;
		}

		// 3. 实时生成平均值 → 写入 colourBuffer（float*）
		for (int i = 0; i < used; i++) {
			int index = indices[i];
			int pixelIndex = index * 3;
			float w = std::max(0.0001f, weightBuffer[index]); // 防止除0

			// 👇 当前像素平均颜色（提供给 OIDN）
			colourBuffer[pixelIndex + 0] = film[index].r / w;
			colourBuffer[pixelIndex + 1] = film[index].g / w;
			colourBuffer[pixelIndex + 2] = film[index].b / w;
		}
	}

	void AOV(int x, int y, const Colour& albedo, const Vec3& normal) {
		if (x < 0 || x >= (int)width || y < 0 || y >= (int)height)
			return;
		int index = (y * width + x) * 3;

		albedoBuffer[index + 0] = albedo.r;
		albedoBuffer[index + 1] = albedo.g;
		albedoBuffer[index + 2] = albedo.b;

		normalBuffer[index + 0] = normal.x;
		normalBuffer[index + 1] = normal.y;
		normalBuffer[index + 2] = normal.z;
	}

	void tonemap(int x, int y, unsigned char& r, unsigned char& g, unsigned char& b, float exposure = 1.0f)
		//普通伽马校正方法
	{
		int idx = y * width + x;
		Colour pixel;
		if (SPP > 20) {
			pixel.r = outputBuffer[idx * 3 + 0];
			pixel.g = outputBuffer[idx * 3 + 1];
			pixel.b = outputBuffer[idx * 3 + 2];
		}
		else {
			float totalWeight = std::max(0.0001f, weightBuffer[idx]);
		    // 使用累计贡献除以累计权重得到平均颜色
		    pixel = film[idx] * (exposure / totalWeight);
		}
		//float totalWeight = std::max(0.0001f, weightBuffer[idx]);
		//// 使用累计贡献除以累计权重得到平均颜色
		//pixel = film[idx] * (exposure / totalWeight);

		r = std::min(powf(std::max(pixel.r, 0.0f), 1.0f / 2.2f) * 255, 255.0f);
		g = std::min(powf(std::max(pixel.g, 0.0f), 1.0f / 2.2f) * 255, 255.0f);
		b = std::min(powf(std::max(pixel.b, 0.0f), 1.0f / 2.2f) * 255, 255.0f);
		// Return a tonemapped pixel at coordinates x, y
	}
	// Do not change any code below this line
	void init(int _width, int _height, ImageFilter* _filter)
	{
		width = _width;
		height = _height;
		film = new Colour[width * height];
		weightBuffer = new float[width * height];
		sppBuffer = new int[width * height];
		albedoBuffer = new float[width * height * 3];
		normalBuffer = new float[width * height * 3];
		outputBuffer = new float[width * height * 3];
		colourBuffer = new float[width * height * 3];
		clear();
		filter = _filter;
	}
	void clear()
	{
		memset(film, 0, width * height * sizeof(Colour));
		memset(weightBuffer, 0, width * height * sizeof(float));
		memset(sppBuffer, 0, width * height * sizeof(int));
		memset(albedoBuffer, 0, width * height * 3 * sizeof(float));
		memset(normalBuffer, 0, width * height * 3 * sizeof(float));
		memset(outputBuffer, 0, width * height * 3 * sizeof(float));
		memset(colourBuffer, 0, width * height * 3 * sizeof(float));
		SPP = 0;
	}
	void incrementSPP()
	{
		SPP++;
	}
	void save(std::string filename)
	{
		Colour* hdrpixels = new Colour[width * height];
		for (unsigned int i = 0; i < (width * height); i++)
		{
			int spp = std::max(1, sppBuffer[i]); // 避免除以 0
			hdrpixels[i] = film[i] * (1.0f / (float)spp);
		}
		stbi_write_hdr(filename.c_str(), width, height, 3, (float*)hdrpixels);
		delete[] hdrpixels;
	}
};

void DenoiseFilm(Film* film)
{
	// 创建默认 CPU 设备（OIDN v2 默认不带参数）
	oidn::DeviceRef device = oidn::newDevice();
	device.commit();

	// 计算图像大小（单位：字节，RGB 3 通道）
	size_t imageSize = film->width * film->height * 3 * sizeof(float);

	// 分别为 color, albedo, normal, output 分配设备缓冲区
	oidn::BufferRef colorBuf = device.newBuffer(imageSize);
	oidn::BufferRef albedoBuf = device.newBuffer(imageSize);
	oidn::BufferRef normalBuf = device.newBuffer(imageSize);
	oidn::BufferRef outputBuf = device.newBuffer(imageSize);

	// 将 CPU 内存中的数据拷贝到设备缓冲区
	memcpy(colorBuf.getData(), film->colourBuffer, imageSize);
	memcpy(albedoBuf.getData(), film->albedoBuffer, imageSize);
	memcpy(normalBuf.getData(), film->normalBuffer, imageSize);

	// 创建 OIDN RT 类型滤镜
	oidn::FilterRef filter = device.newFilter("RT");
	filter.setImage("color", colorBuf, oidn::Format::Float3, film->width, film->height);
	filter.setImage("albedo", albedoBuf, oidn::Format::Float3, film->width, film->height);
	filter.setImage("normal", normalBuf, oidn::Format::Float3, film->width, film->height);
	filter.setImage("output", outputBuf, oidn::Format::Float3, film->width, film->height);
	filter.set("hdr", true);
	filter.commit();
	filter.execute();

	// 将设备中处理后的降噪结果拷贝回 CPU 的 outputBuffer
	memcpy(film->outputBuffer, outputBuf.getData(), imageSize);
}
