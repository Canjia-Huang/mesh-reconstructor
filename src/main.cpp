#include <RPD/RPD.hpp>
#include <IO.hpp>

#include <iostream>
#include <vector>
#include <string>

using namespace std;
using namespace MeshReconstructorRPD;

int main(){
	std::string input_file_path = std::string(OUTPUT_PATH) + "test.obj";
	std::vector<Eigen::Vector3d> points;
	std::vector<Eigen::Vector3d> normals;
	std::vector<double> weights;

	read_points(
		input_file_path,
		points,
		normals,
		weights
	);

	return 1;
}