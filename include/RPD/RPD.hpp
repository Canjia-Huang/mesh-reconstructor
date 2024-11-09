#ifndef MESH_RECONSTRUCTOR_RPD_HPP_
#define MESH_RECONSTRUCTOR_RPD_HPP_

#include <Export.h>
#include <IO.hpp>
#include <RPD/RPDPoint.hpp>
#include <RPD/RPDPlane.hpp>
#include <RPD/RPDCell.hpp>
// Eigen
#include <Eigen/Dense>
// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
// CGAL Poisson Reconstruction
typedef Kernel::Point_3                     K_Point;
typedef Kernel::Vector_3				    K_Vector;
typedef std::pair<K_Point, K_Vector>	    K_Pwn;
#include <CGAL/Polyhedron_3.h>
typedef CGAL::Polyhedron_3<Kernel>		Polyhedron;
#include <CGAL/poisson_surface_reconstruction.h>

namespace MeshReconstructorRPD {
    class Reconstructor{
    public:
        int read_points(const std::string input_file_path){
            DEBUG_ONLY_COUT("");

            // get input points info
            std::vector<Eigen::Vector3d> points;
	        std::vector<Eigen::Vector3d> normals;
	        std::vector<double> weights;
            if (MeshReconstructor::read_points(
                input_file_path,
                points,
                normals,
                weights
            ) == false){
                return 0;
            }

            // process
            if (points.size() == 0){
                DEBUG_ONLY_COUT("input point size is zero!");
                return 0;
            }
            if (normals.size() != points.size()){
                DEBUG_ONLY_COUT("input normal size is invalid, use PCA to estimate");
                if (get_PCA_normals(points, normals) == false){
                    DEBUG_ONLY_COUT("PCA to get normal may error!");
                    return 0;
                }
            }
            if (weights.size() != points.size()){
                DEBUG_ONLY_COUT("input weight size is invalid, all set to zero");
                weights = std::vector<double>(points.size());
            }

            // convert to RPDPoint
            for (int i = 0, i_end = points.size(); i < i_end; ++i){
                points_.emplace_back(RPDPoint(points[i], weights[i]));
            }
            DEBUG_ONLY_COUT("read input point size:" << " " << points_.size());

            return 1;
        }
    
        int get_PCA_normals(
            std::vector<Eigen::Vector3d>& points,
            std::vector<Eigen::Vector3d>& normals
        ){
            

            return 1;
        }

        int build_poisson_base_surface(
            const double r,
            std::vector<Eigen::Vector3d> bs_points,
            std::vector<Eigen::Vector3i> bs_faces
        ){
            DEBUG_ONLY_COUT("");

            // init
            std::vector<Eigen::Vector3d>().swap(bs_points);
            std::vector<Eigen::Vector3i>().swap(bs_faces);
            
            // poisson surface reconstruction
            std::vector<K_Pwn> pwn_points;
		    for (int i = 0, i_end = points_.size(); i < i_end; ++i) {
                K_Point kp(
                    points_[i].x(), 
                    points_[i].y(), 
                    points_[i].z());
                K_Vector kv(
                    points_[i].nx(),
                    points_[i].ny(),
                    points_[i].nz());
                pwn_points.emplace_back(std::make_pair(kp, kv));
		    }

		    Polyhedron output_mesh;
		    if (CGAL::poisson_surface_reconstruction_delaunay(
                    pwn_points.begin(), pwn_points.end(),
                    CGAL::First_of_pair_property_map<K_Pwn>(),
                    CGAL::Second_of_pair_property_map<K_Pwn>(),
                    output_mesh, r) == false
		    ){
                DEBUG_ONLY_COUT("CGAL::poisson_surface_reconstruction_delaunay error!");
                return 0;
            }
            else{
                std::map<Polyhedron::Vertex_handle, int> vh_to_idx;
                for (Polyhedron::Facet_iterator fi = output_mesh.facets_begin(); fi != output_mesh.facets_end(); ++fi) {
                    std::vector<int> face_vertices_idx;
                    Polyhedron::Face_handle fh = fi;

                    Polyhedron::Halfedge_around_facet_circulator fij = fi->facet_begin();
                    do {
                        if (vh_to_idx.find(fij->vertex()) == vh_to_idx.end()) {
                            bs_points.push_back(Eigen::Vector3d(fij->vertex()->point().x(), fij->vertex()->point().y(), fij->vertex()->point().z()));
                            vh_to_idx[fij->vertex()] = bs_points.size() - 1;
                        }
                        face_vertices_idx.push_back(vh_to_idx[fij->vertex()]);
                    } while (++fij != fi->facet_begin());

                    if (face_vertices_idx.size() != 3) {
                        DEBUG_ONLY_COUT("CGAL::poisson_surface_reconstruction_delaunay output error!");
                        return 0;
                    }

                    bs_faces.push_back(Eigen::Vector3i(face_vertices_idx[0], face_vertices_idx[1], face_vertices_idx[2]));
                }
                for (auto vh : vh_to_idx) {
                    bs_points.push_back(Eigen::Vector3d(vh.first->point().x(), vh.first->point().y(), vh.first->point().z()));
                }
            }
        
            return 1;
        }
    private:
        std::vector<RPDPoint> points_;
    };
}; // namespace RPD

#endif
