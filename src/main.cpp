#include <RPD/RPD.hpp>
#include <IO.hpp>

#include <iostream>
#include <vector>
#include <string>

#include <Eigen/Dense>

using namespace std;
/*
int main(){
	std::string input_file_path = std::string(OUTPUT_PATH) + "test.obj";

	MeshReconstructorRPD::Reconstructor Recon;
	Recon.read_points(input_file_path);
    
    std::vector<Eigen::Vector3d> bs_points;
    std::vector<Eigen::Vector3i> bs_faces;
    Recon.build_poisson_surface(0, bs_points, bs_faces);

	return 1;
}
*/