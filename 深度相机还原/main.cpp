#include <iostream>
#include <fstream>
#include <array>
#include "FreeImage.h"       //我觉得是这样，应该添加这个库，因为cmakelist里面特别添加了它的路径，肯定不是空穴来风
#include "Eigen.h"
#include "VirtualSensor.h"

// WANGZHENG : I don't know why we have to generate multipul meshes in while() shouldn't we generate sth like a combination of meshes?
// but I really do't know how to combine the mehes to get better view in meshlab
 
struct Vertex // 翻译应为 ”顶点“ 它包含三维数据和色彩数据，是划分网格和渲染的基本元素
{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	// 下面是4个float 元素的 变量 position 个 4个 无符号char 向量 名为 color

	// position stored as 4 floats (4th component is supposed to be 1.0)
	Eigen::Vector4f position;
	// color stored as 4 unsigned char
	Eigen::Matrix<unsigned char, 4, 1> color;
};

struct YuanShiVertex{

	Eigen::Vector4f pos; //用来储存从虚拟传感器中提取的像素的位置坐标并且用来形成一个 向量
	Eigen::Matrix<unsigned char, 4, 1> color;
};


bool WriteMesh(Vertex* vertices, unsigned int width, unsigned int height, const std::string& filename) {
    float edgeThreshold = 0.01f; // 1cm

    // 顶点总数
    unsigned int nVertices = width * height;
    // 面总数
    unsigned int nFaces = 0;

    std::ofstream outFile(filename);
    if (!outFile.is_open()) return false;

    // 写入文件头
    outFile << "COFF\n";
    outFile << nVertices << " 0 0\n";

    // 写入顶点
    for (unsigned int y = 0; y < height; ++y) {
        for (unsigned int x = 0; x < width; ++x) {
            unsigned int idx = y * width + x;

            if (vertices[idx].position.x() == MINF) {
                outFile << "0 0 0 0 0 0 0\n";
            } else {
                outFile << vertices[idx].position.x() << " "
                        << vertices[idx].position.y() << " "
                        << vertices[idx].position.z() << " "
                        << static_cast<int>(vertices[idx].color[0]) << " "
                        << static_cast<int>(vertices[idx].color[1]) << " "
                        << static_cast<int>(vertices[idx].color[2]) << " "
                        << static_cast<int>(vertices[idx].color[3]) << "\n";
            }
        }
    }

    // 写入三角形
    for (unsigned int y = 0; y < height - 1; ++y) {
        for (unsigned int x = 0; x < width - 1; ++x) {
            unsigned int idx0 = y * width + x;
            unsigned int idx1 = y * width + (x + 1);
            unsigned int idx2 = (y + 1) * width + x;
            unsigned int idx3 = (y + 1) * width + (x + 1);

            auto isValid = [&](unsigned int idx) {
                return vertices[idx].position.x() != MINF;
            };

            auto edgeLength = [&](unsigned int idxA, unsigned int idxB) {
                return (vertices[idxA].position - vertices[idxB].position).norm();
            };

            // 第一个三角形：idx0, idx1, idx2
            if (isValid(idx0) && isValid(idx1) && isValid(idx2) &&
                edgeLength(idx0, idx1) < edgeThreshold &&
                edgeLength(idx1, idx2) < edgeThreshold &&
                edgeLength(idx2, idx0) < edgeThreshold) {
                outFile << "3 " << idx0 << " " << idx1 << " " << idx2 << "\n";
                nFaces++;
            }

            // 第二个三角形：idx1, idx3, idx2
            if (isValid(idx1) && isValid(idx3) && isValid(idx2) &&
                edgeLength(idx1, idx3) < edgeThreshold &&
                edgeLength(idx3, idx2) < edgeThreshold &&
                edgeLength(idx2, idx1) < edgeThreshold) {
                outFile << "3 " << idx1 << " " << idx3 << " " << idx2 << "\n";
                nFaces++;
            }
        }
    }

    // 更新面数量
    outFile.seekp(0, std::ios::beg);
    outFile << "COFF\n";
    outFile << nVertices << " " << nFaces << " 0\n";

    outFile.close();
    return true;
}

int main()
{
	// Make sure this path points to the data folder 
	std::string filenameIn = "/config/workspace/data/rgbd_dataset_freiburg1_xyz/"; // I have to use find in VS to find the abs path it took me fucking long time to debug this
	std::string filenameBaseOut = "WangZheng_mesh";

	// load video 
	std::cout << "Initialize virtual sensor..." << std::endl;
	VirtualSensor sensor;
	if (!sensor.Init(filenameIn))
	{
		std::cout << "Failed to initialize the sensor!\nCheck file path!" << std::endl;
		return -1;
	}

	// convert video to meshes
	while (sensor.ProcessNextFrame())
	{
		// get ptr to the current depth frame
		// depth is stored in row major (get dimensions via sensor.GetDepthImageWidth() / GetDepthImageHeight())
		float* depthMap = sensor.GetDepth();

		// get ptr to the current color frame
		// color is stored as RGBX in row major (4 byte values per pixel, get dimensions via sensor.GetColorImageWidth() / GetColorImageHeight())
		BYTE* colorMap = sensor.GetColorRGBX();

		// this is the dimonsion of each frame
		unsigned int m_depthImageHeight  = sensor.GetColorImageHeight();
		unsigned int m_depthImageWidth = sensor.GetColorImageWidth();

		// get depth intrinsics
		Matrix3f depthIntrinsics = sensor.GetDepthIntrinsics(); // 从sensor 获得内参矩阵

		Matrix3f depthIntrinsicsInv = depthIntrinsics.inverse(); // 深度相机内参逆矩阵

		// 从深度相机内参矩阵中读取数据， 这个是在虚拟传感器头文件里面已经定义好的，只读取
		float fX = depthIntrinsics(0, 0);
		float fY = depthIntrinsics(1, 1);
		float cX = depthIntrinsics(0, 2);
		float cY = depthIntrinsics(1, 2);

		// compute inverse depth extrinsics
		Matrix4f depthExtrinsicsInv = sensor.GetDepthExtrinsics().inverse(); // 获得深度相机外参逆矩阵

		Matrix4f trajectory = sensor.GetTrajectory(); // 获得轨迹的矩阵

		Matrix4f trajectoryInv = sensor.GetTrajectory().inverse(); // 轨迹 的 逆矩阵

		// TODO 1: back-projection
		// write result to the vertices array below, keep pixel ordering!
		// if the depth value at idx is invalid (MINF) write the following values to the vertices array
		// vertices[idx].position = Vector4f(MINF, MINF, MINF, MINF);
		// vertices[idx].color = Vector4uc(0,0,0,0);
		// otherwise apply back-projection and transform the vertex to world space, use the corresponding color from the colormap
		Vertex* vertices = new Vertex[sensor.GetDepthImageWidth() * sensor.GetDepthImageHeight()];
		// 这个是用来储存每一个像素点的向量
		YuanShiVertex* v = new YuanShiVertex[sensor.GetDepthImageWidth() * sensor.GetDepthImageHeight()];

		//接下来时BP 世界下的vertex = 轨迹矩阵的逆 。* 外参的逆矩阵 。* 像素点的列向量（xyz1）‘
		
		for (unsigned int y = 0; y < m_depthImageHeight; ++y) {
    		for (unsigned int x = 0; x < m_depthImageWidth; ++x) {
       		 unsigned int idx = y * m_depthImageWidth + x;

       		 if (depthMap[idx] == 0) continue;  // Check for invalid depth

       		 float depth = depthMap[idx];  // 深度值 (Z)
       		 // 使用相机内参进行反向投影，计算相机坐标系中的 3D 点
       		 float X = (x - cX) * depth / fX;
      	     float Y = (y - cY) * depth / fY ;
       		 float Z = depth;
        	// store the col vec v。pos
        		v[idx].pos = Eigen::Vector4f(X, Y, Z, 1.0f);
  			}
		}
		// embeded the color into vertices
		for(int idx = 0; idx<sensor.GetDepthImageWidth() * sensor.GetDepthImageHeight(); ++idx){
			vertices[idx].position = depthExtrinsicsInv * trajectoryInv  * v[idx].pos;
			for(int i = 0; i<4; ++i){
				vertices[idx].color[i] = colorMap[4*idx + i]; 
			}	  
		}
 
		//每一个frame 在while里面都会生成一个
		std::stringstream ss;
		ss << filenameBaseOut << sensor.GetCurrentFrameCnt() << ".off";
		if (!WriteMesh(vertices, sensor.GetDepthImageWidth(), sensor.GetDepthImageHeight(), ss.str()))
		{
			std::cout << "Failed to write mesh!\nCheck file path!" << std::endl;
			return -1;
		}

		// free mem
		delete[] vertices;
		delete[] v;
	}

	return 0;
}