#ifndef MESH_RECONSTRUCTOR_IO_HPP_
#define MESH_RECONSTRUCTOR_IO_HPP_

#include <Export.h>
#include <fstream>
#include <string>
#include <vector>
// Eigen
#include <Eigen/Dense>

namespace MeshReconstructor {
    int read_OBJ_points(
        const std::string input_file_path,
        std::vector<Eigen::Vector3d>& points,
        std::vector<Eigen::Vector3d>& normals,
        std::vector<double>& weights
    ){
        DEBUG_ONLY_COUT("");

        // init
        std::vector<Eigen::Vector3d>().swap(points);
        std::vector<Eigen::Vector3d>().swap(normals);
        std::vector<double>().swap(weights);

        // read file
        std::ifstream in(input_file_path);
	    if (!in.good()){
            DEBUG_ONLY_COUT("read input file failed!");
            return 0;
        }
        std::string sline, s0;
        std::string vertex_char = "v";
        std::string normal_char = "vt";
        while (getline(in, sline)) {
            std::istringstream ins(sline);
            ins >> s0;

            if (s0 == vertex_char) {
                Eigen::Vector3d p;
                ins >> p.x() >> p.y() >> p.z();
                points.push_back(p);
            }
            else if (s0 == normal_char){
                Eigen::Vector3d n;
                ins >> n.x() >> n.y() >> n.z();
                normals.push_back(n);
            }
        }
        in.close();

        // check input
        DEBUG_ONLY_COUT("read done");
        DEBUG_ONLY_COUT("input point size:" << " " << points.size());
        if (normals.size() > 0) {
            if (points.size() != normals.size()){
            DEBUG_ONLY_COUT("input normals size may error!");
            }
            else {
                DEBUG_ONLY_COUT("input normal size:" << " " << normals.size());
            } 
        }

        return 1;
    }

    int read_XYZ_points(
        const std::string input_file_path,
        std::vector<Eigen::Vector3d>& points,
        std::vector<Eigen::Vector3d>& normals,
        std::vector<double>& weights
    ){
        DEBUG_ONLY_COUT("");

        // init
        std::vector<Eigen::Vector3d>().swap(points);
        std::vector<Eigen::Vector3d>().swap(normals);
        std::vector<double>().swap(weights);

        // read file
        std::ifstream in(input_file_path);
	    if (!in.good()){
            DEBUG_ONLY_COUT("read input file failed!");
            return 0;
        }
        std::string sline;
        while (getline(in, sline)) {
            std::istringstream ins(sline);

            Eigen::Vector3d p;
            ins >> p.x() >> p.y() >> p.z();
            points.push_back(p);

            Eigen::Vector3d n(0, 0, 0);
            ins >> n.x() >> n.y() >> n.z();
            if (abs(n.x()) < MR_EPS && abs(n.y()) < MR_EPS && abs(n.z()) < MR_EPS){
                continue;
            }
            normals.push_back(n);

            double w = -1;
            ins >> w;
            if (w < 0){
                continue;
            }
            weights.push_back(w);
        }

        // check input
        DEBUG_ONLY_COUT("read done");
        DEBUG_ONLY_COUT("input point size:" << " " << points.size());
        if (normals.size() > 0) {
            if (points.size() != normals.size()) {
                DEBUG_ONLY_COUT("input normals size may error!");
            }
            else {
                DEBUG_ONLY_COUT("input normal size:" << " " << normals.size());
            }
        }
        if (weights.size() > 0) {
            if (points.size() != weights.size()) {
                DEBUG_ONLY_COUT("input weights size may error!");
            }
            else {
                DEBUG_ONLY_COUT("input weight size:" << " " << weights.size());
            }
        }

        return 1;
    }

    int read_OFF_points(
        const std::string input_file_path,
        std::vector<Eigen::Vector3d>& points,
        std::vector<Eigen::Vector3d>& normals,
        std::vector<double>& weights
    ){
        DEBUG_ONLY_COUT("");

        // init
        std::vector<Eigen::Vector3d>().swap(points);
        std::vector<Eigen::Vector3d>().swap(normals);
        std::vector<double>().swap(weights);

        std::string OFF_identify;
        int vertex_num;
		int face_num;
		int edge_num;
		int s;
        std::string sline;

        // read file
        std::ifstream in(input_file_path);
        if (!in.good()) {
            DEBUG_ONLY_COUT("read input file failed!");
            return 0;
        }
        {
            getline(in, sline);
            std::istringstream ins(sline);
            ins >> OFF_identify;
            if (OFF_identify != "OFF"){
                DEBUG_ONLY_COUT("off file invalid!");
                return 0;
            }
        }

        {
            getline(in, sline);
            std::istringstream ins(sline);
            ins >> vertex_num >> face_num >> edge_num;
            for (int i = 0; i < vertex_num; ++i) {
                getline(in, sline);
                std::istringstream ins_s(sline); 

                Eigen::Vector3d p;
                ins_s >> p.x() >> p.y() >> p.z();
                points.push_back(p);
            }
        }
        
        in.close();

        DEBUG_ONLY_COUT("read done");
        DEBUG_ONLY_COUT("input point size:" << " " << points.size());

        return 1;
    }

    int read_points(
        const std::string input_file_path,
        std::vector<Eigen::Vector3d>& points,
        std::vector<Eigen::Vector3d>& normals,
        std::vector<double>& weights
    ){
        DEBUG_ONLY_COUT("");

        std::string back = input_file_path.substr(input_file_path.length() - 3, input_file_path.length());
        if (back == "obj"){
            return read_OBJ_points(
                input_file_path,
                points,
                normals,
                weights
            );
        }
        else if (back == "xyz"){
            return read_XYZ_points(
                input_file_path,
                points,
                normals,
                weights
            );
        }
        else if (back == "off"){
            return read_OFF_points(
                input_file_path,
                points,
                normals,
                weights
            );
        }
        else{
            DEBUG_ONLY_COUT("input file type invalid!");
            return 0;
        }

        return 1;
    }

}

#endif