#include <RPD/RPD.hpp>
#include <IO.hpp>

#include <iostream>
#include <vector>
#include <string>

using namespace std;

int main(){
	std::string input_file_path = std::string(OUTPUT_PATH) + "test.obj";

	MeshReconstructorRPD::Reconstructor Recon;
	Recon.read_points(input_file_path);

	return 1;
}