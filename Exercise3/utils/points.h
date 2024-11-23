#pragma once

#include <vector>
#include <fstream>
#include <istream>
#include <iostream>
#include <sstream>


/*
struct Weight
{
	Weight() = default;

	Weight(double w_) : 
		w{ w_ }
	{
	}

	static Weight from_string(const std::string& line)
	{
		std::istringstream stream(line);
		double w{};
		stream >> w;

		Weight weight{ w };
		return weight;
	}

	double w{};
}; */
struct Weight {
    Weight() = default;

    // 使用 vector<double> 初始化 w
    Weight(const std::vector<double>& weights) : w(weights) {}

    // 从字符串构建 Weight 对象，支持多个权重值
    static Weight from_string(const std::string& line) {
        std::istringstream stream(line);
        std::vector<double> weights;
        double value;

        // 逐个读取 double 值并添加到 weights 中
        while (stream >> value) {
            weights.push_back(value);
        }

        return Weight(weights);
    }

    std::vector<double> w; // 支持多个权重值
};

struct Point2D
{
	Point2D() = default;
	
	Point2D(double x_, double y_)
		: x{ x_ }, y{ y_ }
	{
	}

	static Point2D from_string(const std::string& line)
	{
		std::istringstream stream(line);
		float x{};
		float y{};
		stream >> x >> y;

		Point2D point{ x, y };
		return point;
	}

	double x{};
	double y{};
};

struct Point3D
{
	Point3D(double x_, double y_, double z_)
		: x{ x_ }, y{ y_ }, z{ z_ }
	{
	}

	static Point3D from_string(const std::string& line)
	{
		std::istringstream stream(line);
		float x{};
		float y{};
		float z{};
		stream >> x >> y >> z;

		Point3D point{ x, y, z };
		return point;
	}

	double x{};
	double y{};
	double z{};
};

using PointList2D = std::vector<Point2D>;
using PointList3D = std::vector<Point3D>;