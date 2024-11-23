#pragma once

#ifndef IMPLICIT_SURFACE_H
#define IMPLICIT_SURFACE_H
#include <iostream>
#include "Eigen.h"
#include "SimpleMesh.h"
#include <cmath>

class ImplicitSurface
{
public:
	virtual double Eval(const Eigen::Vector3d& x) = 0;
};


class Sphere : public ImplicitSurface
{
public:
	Sphere(const Eigen::Vector3d& center, double radius) : m_center(center), m_radius(radius)
	{
	}

	double Eval(const Eigen::Vector3d& _x)
	{
		// TODO: implement the implicit sphere formula using the member variables m_center and m_radius
		double SDF = pow((_x(0) - m_center(0)), 2) + pow((_x(1) - m_center(1)), 2) + pow((_x(2) - m_center(2)), 2) - pow(m_radius, 2);
  	    return SDF;
		
	}

private:
	Eigen::Vector3d m_center;
	double m_radius;
};


class Torus : public ImplicitSurface
{
public:
	Torus(const Eigen::Vector3d& center, double radius, double a) : m_center(center), m_radius(radius), m_a(a)
	{
	}

	double Eval(const Eigen::Vector3d& _x)
	{
		// TODO: implement the implicit torus formula using the  variables m_center, m_radius (radius of the ring) and the radius m_a (small radius)
		// 我看出来了，就是计算出来真正的在catisian坐标下的 对应的函数值，比如传入的坐标x离与以center为中心的原不alinment， 
		double transit = pow((_x(0) - m_center(0)), 2) +
		 				 pow((_x(1) - m_center(1)), 2) + 
		 				 pow((_x(2) - m_center(2)), 2) +
						 pow(m_radius,2) - 
						 pow(m_a,2);
		double SDF = pow(transit,2) - 
				   4*pow(m_radius,2)*(pow((_x(0) - m_center(0)), 2) + pow((_x(1) - m_center(1)), 2));
		// 这里的甜甜圈方程(x2 + y2 + z2 + R2 - a2)2 - 4*R2(x2 + y2)
		return SDF;
	}

private:
	Eigen::Vector3d m_center;
	double m_radius;
	double m_a;
};


class Hoppe : public ImplicitSurface
{
public:
	Hoppe(const std::string& filenamePC)
	{
		//std::cout << "start reading the big point cloud file..." << std::endl;
		m_pointcloud.ReadFromFile(filenamePC); // pointcloud的实现在simple mesh的第91行
	}
		
	//std::cout << "calculating the eval ..." << std::endl;
	double Eval(const Eigen::Vector3d& _x)
	{
	
    Eigen::Vector3f x = Eigen::Vector3f((float)_x.x(), (float)_x.y(), (float)_x.z());
    unsigned int idx = m_pointcloud.GetClosestPoint(x);
    if (idx == m_pointcloud.GetPoints().size()) {
        return std::numeric_limits<double>::max();
    }
	//std::cout << "done getting point cloud" << std::endl;
    Eigen::Vector3f p = m_pointcloud.GetPoints()[idx];
    Eigen::Vector3f n = m_pointcloud.GetNormals()[idx];
	//std::cout << "calculating result..." << std::endl;
    return (x - p).dot(n);  // dot 返回的是 float，自动转换为 double
	}

private:
	PointCloud m_pointcloud;
};

///////////////////////////////////////////

class FunctionSamples
{
public:

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	FunctionSamples() {}

	void insertSample(const Vector3d& pos, const double targetFunctionValue)
	{
		m_pos.push_back(pos);
		m_val.push_back(targetFunctionValue);
	}

	std::vector<Vector3d> m_pos;
	std::vector<double> m_val;
};

#ifdef ORIGINAL
class RBF : public ImplicitSurface
{
public:
	RBF(const std::string& filenamePC)
	{
		// load point cloud
		m_pointcloud.ReadFromFile(filenamePC);

		// Create function samples
		double eps = 0.01f;                  
		for (unsigned int i = 0; i < m_pointcloud.GetPoints().size(); i++)
		{
			Eigen::Vector3f ptF = m_pointcloud.GetPoints()[i]; //std::vector<Eigen::Vector3f>& 是GetPoints（）的返回类型
			Eigen::Vector3f nF = m_pointcloud.GetNormals()[i]; // getnormal的返回值也是个vector
			const Vector3d pt(ptF[0], ptF[1], ptF[2]);
			const Vector3d n(nF[0], nF[1], nF[2]);

			m_funcSamp.insertSample(pt, 0);
			m_funcSamp.insertSample(pt + n*eps, eps); 
			eps *= -1;
		}

		m_numCenters = m_funcSamp.m_pos.size();
		const unsigned int dim = m_numCenters + 4;

		// build and solve the linear system of equations
		m_systemMatrix = MatrixXd(dim, dim);
		m_rhs = VectorXd(dim);
		m_coefficents = VectorXd(dim); // result of the linear system
		BuildSystem();
		SolveSystem();
	}

	double Eval(const Eigen::Vector3d& _x){
    	double result = 0.0;

   		// 遍历所有 RBF 核函数的中心点
   		for (unsigned int i = 0; i < m_numCenters; ++i)
   		{
     		   // 计算当前点与第 i 个中心点之间的欧氏距离
      		double distance = (_x - m_funcSamp.m_pos[i]).norm();

      		  // 使用 EvalBasis 计算径向基函数值
       		double basisValue = EvalBasis(distance);

       		 // 将核函数的贡献加入结果，使用系数 m_coefficents[i]
        		result += m_coefficents[i] * basisValue;
    	}

    	// 添加线性部分和常数项
    	result += m_coefficents[m_numCenters] * _x.x();
    	result += m_coefficents[m_numCenters + 1] * _x.y();
    	result += m_coefficents[m_numCenters + 2] * _x.z();
    	result += m_coefficents[m_numCenters + 3];

    	return result;
	}

private:

	double EvalBasis(double x)
	{
		return x*x*x;
	}

	void CompResidual()
	{
		double residual = 0;
		for (int x = 0; x < m_numCenters; x++)
		{
			Vector3d pt = m_funcSamp.m_pos[x];
			residual += fabs(Eval(pt) - m_funcSamp.m_val[x]);
		}
		residual /= (double)m_numCenters;

		std::cerr << "Residual: " << residual << std::endl;
	}

#define phi(i,j) EvalBasis((m_funcSamp.m_pos[i]-m_funcSamp.m_pos[j]).norm())
	//! Computes the system matrix.
	void BuildSystem()
	{
		m_systemMatrix.setZero();
		m_rhs.setZero();

		for (int i = 0; i<m_numCenters; i++) {
			for (int j = 0; j<m_numCenters; j++) 
			
			m_systemMatrix(i, j) = phi(i, j);

			m_systemMatrix.row(i).segment(m_numCenters, 3) = m_funcSamp.m_pos[i];
			m_systemMatrix.row(i)(m_numCenters + 3) = 1;

			m_systemMatrix.col(i).segment(m_numCenters, 4) = m_systemMatrix.row(i).segment(m_numCenters, 4);

			m_rhs(i) = m_funcSamp.m_val[i];
		}

		// regularizer -> smoother surface
		m_systemMatrix.diagonal() += 0.001 * VectorXd::Ones(m_numCenters + 4);
	}

	void SolveSystem()
	{
		std::cerr << "Solving RBF System" << std::endl;
		std::cerr << "Computing LU..." << std::endl;

		FullPivLU<Matrix<double, Dynamic, Dynamic>> LU(m_systemMatrix);
		m_coefficents = LU.solve(m_rhs);

		std::cerr << "Done." << std::endl;
	}

	// point cloud
	PointCloud m_pointcloud;

	//! The given function samples (at each function sample, we place a basis function).
	FunctionSamples m_funcSamp;

	//! The number of center = number of function samples.
	unsigned int m_numCenters;

	//! the right hand side of our system of linear equation. Unfortunately, float-precision is not enough, so we have to use double here.
	VectorXd m_rhs;

	//! the system matrix. Unfortunately, float-precision is not enough, so we have to use double here
	MatrixXd m_systemMatrix;

	//! store the result of the linear system here. Unfortunately, float-precision is not enough, so we have to use double here
	VectorXd m_coefficents;
};

#else

class RBF : public ImplicitSurface
{
public:
	RBF(const std::string& filenamePC)
	{
		std::cerr << "Starting RBF initialization..." << std::endl;

		// load point cloud
		std::cerr << "Reading point cloud from file: " << filenamePC << std::endl;
		 bool success = m_pointcloud.ReadFromFile(filenamePC);
		//m_pointcloud.ReadFromFile(filenamePC);
		if (!success) {
        	std::cerr << "ERROR: Failed to read point cloud file!" << std::endl;
        	return;
   		}
    	std::cerr << "Point cloud loaded successfully." << std::endl;

		// Create function samples
		double eps = 0.01f;
		std::cerr << "building first func sample" << std::endl;
		// on surface points (-> center points of the RBFs)
		for (unsigned int i = 0; i < m_pointcloud.GetPoints().size(); i++)
		{
			Eigen::Vector3f ptF = m_pointcloud.GetPoints()[i];
			Eigen::Vector3f nF = m_pointcloud.GetNormals()[i];
			const Vector3d pt(ptF[0], ptF[1], ptF[2]);
			const Vector3d n(nF[0], nF[1], nF[2]);

			m_funcSamp.insertSample(pt, 0); // on surface point => distance = 0
		}
		std::cerr << "first func sample done" << std::endl;
		// off surface points
		for (unsigned int i = 0; i < m_pointcloud.GetPoints().size(); i++)
		{
			Eigen::Vector3f ptF = m_pointcloud.GetPoints()[i];
			Eigen::Vector3f nF = m_pointcloud.GetNormals()[i];
			const Vector3d pt(ptF[0], ptF[1], ptF[2]);
			const Vector3d n(nF[0], nF[1], nF[2]);

			m_funcSamp.insertSample(pt + n*eps, eps);// off surface point => distance = eps
			eps *= -1;
		}
		 std::cerr << "Function samples created." << std::endl;
		m_numCenters = (unsigned int) m_pointcloud.GetPoints().size();
		const unsigned int dim = m_numCenters + 4;

		// build and solve the linear system of equations
		m_systemMatrix = MatrixXd(dim, dim);
		m_rhs = VectorXd(dim);
		m_coefficents = VectorXd(dim); // result of the linear system
		std::cerr << "Building system matrix..." << std::endl;
		BuildSystem();
		std::cerr << "System matrix built successfully." << std::endl;
		std::cerr << "Solving system matrix..." << std::endl;
		SolveSystem();
		std::cerr << "System matrix solved successfully." << std::endl;
		std::cerr << "RBF initialization completed." << std::endl;
	}

	double Eval(const Eigen::Vector3d& _x)
	{
		// TODO: eval the RBF function based on the coefficents stored in m_coefficents
		// the first m_numCenters entries contain the coefficients of the kernels (that can be evaluated using EvalBasis())
		// the following parameters are the coeffients for the linear and the constant part
		// the centers of the RBFs are the first m_numCenters sample points (use m_funcSamp.m_pos[i] to access them)
		// hint: Eigen provides a norm() function to compute the l2-norm of a vector (e.g. see macro phi(i,j))
		double result = 0.0;

   		 // go through all point in kern 
    	for (unsigned int i = 0; i < m_numCenters; ++i)
   		 {
        	// calculate distance usingnorm
        	double distance = (_x - m_funcSamp.m_pos[i]).norm();

        	// use EvalBasis 计算径向基函数值
       		 double basisValue = EvalBasis(distance);

        	// use m_coefficents[i] to count them in result
        	result += m_coefficents[i] * basisValue;
    	}

    	// add linear part of the RBF
    	result += m_coefficents[m_numCenters] * _x.x();
    	result += m_coefficents[m_numCenters + 1] * _x.y();

		return result;
	}

private:

	double EvalBasis(double x) // 看来RBF的核儿 就是 简单的差的三次方而已
	{
		return x*x*x;
	}

 #define phi(i,j) EvalBasis((m_funcSamp.m_pos[i]-m_funcSamp.m_pos[j]).norm()) // 这个eval basis 就是phi（ij）宏定义的实现

	//! Computes the system matrix.
	void BuildSystem()
	{
		MatrixXd A(2 * m_numCenters, m_numCenters + 4);
		VectorXd b(2 * m_numCenters);
		A.setZero();
		b.setZero();

		// TODO fill the matrix A and the vector b as described in the exercise sheet
		// note that all sample points (both on and off surface points) are stored in m_funcSamp
		// you can access matrix elements using for example A(i,j) for the i-th row and j-th column
		// similar you access the elements of the vector b, e.g. b(i) for the i-th element
		for (unsigned int i = 0; i < m_funcSamp.m_pos.size(); ++i) // 。m——pos是一个包含xyz数据的std：：vector
		{
 		   // 计算每一个样本点与每一个中心点之间的核函数值，并填入 A 的前 N 列， 也就是课件上面说过的那个AC=B那个公式描述的
  			for (unsigned int j = 0; j < m_numCenters; ++j)
   		    {
   			    A(i, j) = phi(i, j); // phi(i, j) 使用宏定义来计算径向基函数的值
   			}

    		// 他妈的 在eigen库只能生成方阵，我还得挨个往A里面添加 我真是服了
    		A(i, m_numCenters) = m_funcSamp.m_pos[i].x(); 
    		A(i, m_numCenters + 1) = m_funcSamp.m_pos[i].y();  
    		A(i, m_numCenters + 2) = m_funcSamp.m_pos[i].z(); 
    		A(i, m_numCenters + 3) = 1.0;  // 因为后面要x常数项d所以在A中对应位置需要是1.

    		// 填充 b 向量
   			 b(i) = m_funcSamp.m_val[i];  // 这个就是
		}


		// 构造一个方阵 ata 来为LU分解做铺垫
		m_systemMatrix = A.transpose() * A; 
		m_rhs = A.transpose() * b;

		// regularizer -> smoother surface
		// pushes the coefficients to zero
		double lambda = 0.0001;
		m_systemMatrix.diagonal() += lambda * lambda * VectorXd::Ones(m_numCenters + 4);
	}

	void SolveSystem()
	{
		std::cerr << "Solving RBF System" << std::endl;
		std::cerr << "Computing LU..." << std::endl;

		FullPivLU<Matrix<double, Dynamic, Dynamic>> LU(m_systemMatrix);
		m_coefficents = LU.solve(m_rhs);

		std::cerr << "Done." << std::endl;
	}

	// point cloud
	PointCloud m_pointcloud;

	//! The given function samples (at each function sample, we place a basis function).
	FunctionSamples m_funcSamp;

	//! The number of center = number of function samples.
	unsigned int m_numCenters;

	//! the right hand side of our system of linear equation. Unfortunately, float-precision is not enough, so we have to use double here.
	VectorXd m_rhs;

	//! the system matrix. Unfortunately, float-precision is not enough, so we have to use double here
	MatrixXd m_systemMatrix;

	//! store the result of the linear system here. Unfortunately, float-precision is not enough, so we have to use double here
	VectorXd m_coefficents;
};

#endif

#endif //IMPLICIT_SURFACE_H
