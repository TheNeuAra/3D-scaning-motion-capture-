#include <iostream>
#include "Eigen.h"
#include "ImplicitSurface.h"
#include "Volume.h"
#include "MarchingCubes.h"

int main()
{
	std::string filenameIn = "/config/workspace/data/pointcloud.off";
	std::string filenameOut = "hoppe.off";

	// implicit surface
	 std::cerr << "Initializing implicit surface..." << std::endl;
	ImplicitSurface* surface;
	// TODO: you have to switch between these surface types
	//surface = new Sphere(Eigen::Vector3d(0.5, 0.5, 0.5), 0.4);
	//surface = new Torus(Eigen::Vector3d(0.5, 0.5, 0.5), 0.4, 0.1);
	surface = new Hoppe(filenameIn);
	//surface = new RBF(filenameIn);
	std::cerr << "Implicit surface initialized successfully." << std::endl;
	// fill volume with signed distance values  原来SDF 是在这里生成的
	// resolution of the grid, for debugging you can reduce the resolution (-> faster)
	unsigned int mc_res = 5; // 减少了分辨率， 原来是50 的 mc_res

	Volume vol(Vector3d(-0.1,-0.1,-0.1), Vector3d(1.1,1.1,1.1), mc_res, mc_res, mc_res, 1); 
	std::cerr << "Volume initialized from (-0.1, -0.1, -0.1) to (1.1, 1.1, 1.1)" << std::endl;
	// 生成了一个从-0.1 到1.1的三维空间 
//---------------------------------------------------------------------------
for (unsigned int x = 0; x < vol.getDimX(); x++)
{
    std::cerr << "Processing x = " << x << " / " << vol.getDimX() << std::endl;
    
    for (unsigned int y = 0; y < vol.getDimY(); y++)
    {
        std::cerr << "  Processing y = " << y << " / " << vol.getDimY() << std::endl;
        
        for (unsigned int z = 0; z < vol.getDimZ(); z++)
        {
            std::cerr << "    Processing z = " << z << " / " << vol.getDimZ() << std::endl;
            
            // 转换索引到笛卡尔坐标
            Eigen::Vector3d p = vol.pos(x, y, z);
            std::cerr << "      Point in Cartesian coordinates: " << p.transpose() << std::endl;
            
            // 评估隐式函数值
            double val = surface->Eval(p);
            std::cerr << "      Eval result for point (" << x << ", " << y << ", " << z << "): " << val << std::endl;

            // 设置体素值
            vol.set(x, y, z, val);
            std::cerr << "      Set volume value at (" << x << ", " << y << ", " << z << ") to " << val << std::endl;
        }
    }
}
//---------------------------------------------------------------------------------------------
std::cerr << "Completed processing all volume cells." << std::endl;


	// extract the zero iso-surface using marching cubes
	SimpleMesh mesh;
	for (unsigned int x = 0; x < vol.getDimX() - 1; x++)
	{
		std::cerr << "Marching Cubes on slice " << x << " of " << vol.getDimX() << std::endl;

		for (unsigned int y = 0; y < vol.getDimY() - 1; y++)
		{
			for (unsigned int z = 0; z < vol.getDimZ() - 1; z++)
			{
				ProcessVolumeCell(&vol, x, y, z, 0.00f, &mesh);
			}
		}
	}

	// write mesh to file
	if (!mesh.WriteMesh(filenameOut))
	{
		std::cout << "ERROR: unable to write output file!" << std::endl;
		return -1;
	}

	delete surface;
	return 0;
}