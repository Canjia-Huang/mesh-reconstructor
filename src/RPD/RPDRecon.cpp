#include "RPD.h"

#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <set>
#include <cmath>

// for debug
// #define RPD_DEBUG
std::string test_file_path = "..//data//RPD_";

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Weighted_point_3 K_WPoint;
typedef K::Point_3			K_Point;
typedef K::Vector_3         K_Vector;
typedef std::pair<K_Point, K_Vector> K_Pwn;
#include <CGAL/Regular_triangulation_3.h>
typedef CGAL::Regular_triangulation_3<K> Regular_triangulation;
#include <CGAL/Weighted_point_3.h>
#include <CGAL/poisson_surface_reconstruction.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/Polyhedron_3.h>
typedef CGAL::Polyhedron_3<K> Polyhedron;
#include <CGAL/Timer.h>

// nanoflann
#include <nanoflann/nanoflann.hpp>
#include <nanoflann/utils.h>
typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud<double> >, PointCloud<double>, 3 /* dim */> my_kd_tree_t;

#define PI 3.1415926
#define EPS 1e-16
#define MAX 1e16
#define SIGN(x) ((x<0)? -1 : 1)

namespace RPD
{
    RPDRecon::RPDRecon()
    {
    }

    RPDRecon::RPDRecon(std::vector<Eigen::Vector3d> cors,
        std::vector<Eigen::Vector3d> nors,
        std::vector<bool> is_features)
    {
        if (cors.size() == nors.size() && nors.size() == is_features.size())
        {
            for (int i = 0; i < cors.size(); i++)
            {
                if (is_features[i])
                {
                    points_.push_back(RPDPoint(cors[i], nors[i]));
                    is_features_.push_back(true);
                }
                else
                {
                    points_.push_back(RPDPoint(cors[i], nors[i]));
                    is_features_.push_back(false);
                }

            }
        }
        else
        {
            RPD_DEBUG_ONLY(std::cout << "RPDrecon::RPDrecon input vector size error!" << std::endl;)
            return;
        }

#ifdef RPD_DEBUG
        if (1)
        {
            std::ofstream out(test_file_path + "input_points.obj");
            for (int i = 0; i < points_.size(); i++)
            {
                out << "v" << " " << points_[i].cor_.transpose() << std::endl;
            }
            out.close();
        }
#endif
    }

    int RPDRecon::input_points(std::vector<Eigen::Vector3d> cors,
        std::vector<Eigen::Vector3d> nors,
        std::vector<bool> is_features)
    {
        if (cors.size() == nors.size() && nors.size() == is_features.size())
        {
            for (int i = 0; i < cors.size(); i++)
            {
                if (is_features[i])
                {
                    points_.push_back(RPDPoint(cors[i], nors[i], feature_weight_));
                    is_features_.push_back(true);
                }
                else
                {
                    points_.push_back(RPDPoint(cors[i], nors[i], not_feature_weight_));
                    is_features_.push_back(false);
                }

            }
        }
        else
        {
            RPD_DEBUG_ONLY(std::cout << "RPDrecon::RPDrecon input vector size error!" << std::endl;)
                return 0;
        }
#ifdef RPD_DEBUG
        if (1)
        {
            std::ofstream out(test_file_path + "input_points.obj");
            for (int i = 0; i < points_.size(); i++)
            {
                out << "v" << " " << points_[i].cor_.transpose() << std::endl;
            }
            out.close();
        }
#endif

#ifdef RPD_DEBUG
        if (1)
        {
            RPDRecon::output_feature_points(test_file_path + "input_points_features.obj", 255, 0, 0);
        }
#endif

        return 1;
    }

    double RPDRecon::radius()
    {
        return radius_;
    }

    void RPDRecon::set_radius(double r)
    {
        RPD_DEBUG_ONLY(if (r <= 0) std::cout << "RPDRecon::set_radius radius may invalid!" << std::endl;)

        radius_ = r;
    }

    void RPDRecon::set_feature_weight(double fw)
    {
        feature_weight_ = fw;

        for (int i = 0; i < points_.size(); i++)
        {
            if (is_features_[i])
            {
                points_[i].weight_ = fw;
            }
        }

#ifdef RPD_DEBUG
        if (1)
        {
            RPDRecon::output_feature_points(test_file_path + "input_points_features.obj", 255, 0, 0);
        }
#endif
    }

    void RPDRecon::set_not_feature_weight(double nfw)
    {
        not_feature_weight_ = nfw;

        for (int i = 0; i < points_.size(); i++)
        {
            if (is_features_[i] == false)
            {
                points_[i].weight_ = nfw;
            }
        }

#ifdef RPD_DEBUG
        if (1)
        {
            RPDRecon::output_feature_points(test_file_path + "input_points_features.obj", 255, 0, 0);
        }
#endif
    }

    double RPDRecon::safetyAcos(double value)
    {
        if (value < -1.0f) value = -1.0f;
        if (value > 1.0f) value = 1.0f;
        return acos(value) * 180 / PI;
    }

    bool RPDRecon::is_feature_point(RPDPoint rp)
    {
        return (rp.weight_ == feature_weight_);
    }

    bool RPDRecon::is_feature_point(int i)
    {
        if (i < 0 || i > is_features_.size())
        {
            RPD_DEBUG_ONLY(std::cout << "RPDRecon::is_feature_point input invalid!" << std::endl;)
            return false;
        }

        return is_features_[i];
    }

    // functional
    int RPDRecon::PoissonSurfaceReconstruction(std::vector<Eigen::Vector3d>& points,
        std::vector<Eigen::Vector3i>& faces)
    {
        std::string test_func_name = "PS_";

        if (points_.size() == 0)
        {
            RPD_DEBUG_ONLY(std::cout << "PoissonSurfaceReconstruction::input point array is empty!" << std::endl;)
            return 0;
        }

        CGAL::Timer time;
        time.start();

        // init
        points.clear();
        faces.clear();

        std::vector<K_Pwn> pwn_points;
        for (int i = 0; i < points_.size(); i++)
        {
            K_Point kp(points_[i].cor_.x(), points_[i].cor_.y(), points_[i].cor_.z());
            K_Vector kv(points_[i].nor_.x(), points_[i].nor_.y(), points_[i].nor_.z());
            pwn_points.push_back(std::make_pair(kp, kv));
        }

        double average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>
            (pwn_points, 6, CGAL::parameters::point_map(CGAL::First_of_pair_property_map<K_Pwn>()));
        average_spacing *= 0.5;

        Polyhedron output_mesh;
        if (CGAL::poisson_surface_reconstruction_delaunay
        (pwn_points.begin(), pwn_points.end(),
            CGAL::First_of_pair_property_map<K_Pwn>(),
            CGAL::Second_of_pair_property_map<K_Pwn>(),
            output_mesh, average_spacing))
        {
            std::map<Polyhedron::Vertex_handle, int> vh_to_idx;
            for (Polyhedron::Facet_iterator fi = output_mesh.facets_begin(); fi != output_mesh.facets_end(); fi++)
            {
                std::vector<int> face_vertices_idx;
                Polyhedron::Face_handle fh = fi;

                Polyhedron::Halfedge_around_facet_circulator fij = fi->facet_begin();
                do {
                    if (vh_to_idx.find(fij->vertex()) == vh_to_idx.end())
                    {
                        points.push_back(Eigen::Vector3d(fij->vertex()->point().x(), fij->vertex()->point().y(), fij->vertex()->point().z()));
                        vh_to_idx[fij->vertex()] = points.size() - 1;
                    }
                    face_vertices_idx.push_back(vh_to_idx[fij->vertex()]);
                } while (++fij != fi->facet_begin());

                if (face_vertices_idx.size() != 3)
                {
                    RPD_DEBUG_ONLY(std::cout << "PoissonSurfaceReconstruction::CGAL::poisson_surface_reconstruction_delaunay output error!" << std::endl;)
                    return 0;
                }

                faces.push_back(Eigen::Vector3i(face_vertices_idx[0], face_vertices_idx[1], face_vertices_idx[2]));
            }
            for (auto vh : vh_to_idx)
            {
                points.push_back(Eigen::Vector3d(vh.first->point().x(), vh.first->point().y(), vh.first->point().z()));
            }

#ifdef RPD_DEBUG
            if (1)
            {
                std::ofstream out(test_file_path + test_func_name + "surface.obj");
                for (int i = 0; i < points.size(); i++)
                {
                    out << "v" << " " << points[i].x() << " " << points[i].y() << " " << points[i].z() << std::endl;
                }
                for (int i = 0; i < faces.size(); i++)
                {
                    out << "f" << " " << faces[i].x() + 1 << " " << faces[i].y() + 1 << " " << faces[i].z() + 1 << std::endl;
                }
                out.close();
            }
#endif
        }
        else
        {
            RPD_DEBUG_ONLY(std::cout << "PoissonSurfaceReconstruction::CGAL::poisson_surface_reconstruction_delaunay error!" << std::endl;)
            return 0;
        }

        if (points.size() == 0 || faces.size() == 0)
        {
            RPD_DEBUG_ONLY(std::cout << "PoissonSurfaceReconstruction::CGAL::poisson_surface_reconstruction_delaunay null!" << std::endl;)
            return 0;
        }

        time.stop();


        // record time count
        RecordTimeCount("Poisson_Surface_Reconstruction", time.time());

        return 1;
    }

	int RPDRecon::RegularTriangulation(std::vector<std::vector<int>>& RT_neighbors_o)
	{
        std::string test_func_name = "RT_";

        if (radius_ == 0)
        {
            RPD_DEBUG_ONLY(std::cout << "RegularTriangulation::radius is invalid!" << std::endl;)
            return 0;
        }
        if (points_.size() == 0)
        {
            RPD_DEBUG_ONLY(std::cout << "RegularTriangulation::input point array is empty!" << std::endl;)
            return 0;
        }

        CGAL::Timer time;
        time.start();

        // init
        RT_neighbors_o.clear();
        std::map<RPDPoint, int> point_to_idx;
		std::vector<Regular_triangulation::Weighted_point> w_points;
		for (int i = 0; i < points_.size(); i++)
		{
            point_to_idx[points_[i]] = i;
			w_points.push_back(K_WPoint(K_Point(points_[i].cor_.x(), points_[i].cor_.y(), points_[i].cor_.z()), pow(points_[i].weight_, 2)));
		}

        // regular triangulation
		Regular_triangulation RT(w_points.begin(), w_points.end());

        if (!RT.is_valid())
        {
            RPD_DEBUG_ONLY(std::cout << "RegularTriangulation::Regular_triangulation error!" << std::endl;)
            return 0;
        }

#ifdef RPD_DEBUG
        if (1)
        {
            std::ofstream rt_vis_out(test_file_path + test_func_name + "connection.obj");
            std::map<RPDPoint, int> point_to_idx;
            int rt_cnt = 1;
            for (const Regular_triangulation::Vertex_handle vh : RT.finite_vertex_handles())
            {
                std::vector<Regular_triangulation::Vertex_handle> f_vertices;
                RT.finite_adjacent_vertices(vh, std::back_inserter(f_vertices));

                if (point_to_idx.find(RPDPoint(vh->point().x(), vh->point().y(), vh->point().z())) == point_to_idx.end())
                {
                    point_to_idx[RPDPoint(vh->point().x(), vh->point().y(), vh->point().z())] = rt_cnt;
                    rt_vis_out << "v" << " " << vh->point().x() << " " << vh->point().y() << " " << vh->point().z() << std::endl;
                    rt_cnt++;
                }

                for (auto nb : f_vertices)
                {
                    if (point_to_idx.find(RPDPoint(nb->point().x(), nb->point().y(), nb->point().z())) == point_to_idx.end())
                    {
                        point_to_idx[RPDPoint(nb->point().x(), nb->point().y(), nb->point().z())] = rt_cnt;
                        rt_vis_out << "v" << " " << nb->point().x() << " " << nb->point().y() << " " << nb->point().z() << std::endl;
                        rt_cnt++;
                    }

                    rt_vis_out << "l" << " " << point_to_idx[RPDPoint(vh->point().x(), vh->point().y(), vh->point().z())] << " " 
                        << point_to_idx[RPDPoint(nb->point().x(), nb->point().y(), nb->point().z())] << std::endl;
                }
            }
            rt_vis_out.close();
        }
#endif

        // process regular triangulation
        std::vector<std::vector<int>> RT_neighbors;
        for (int i = 0; i < points_.size(); i++)
        {
            RT_neighbors.push_back(std::vector<int>());
        }
        for (const Regular_triangulation::Vertex_handle vh : RT.finite_vertex_handles())
        {
            int vh_idx = point_to_idx[RPDPoint(vh->point().x(), vh->point().y(), vh->point().z())];

            std::vector<Regular_triangulation::Vertex_handle> vh_neighbor;
            RT.finite_adjacent_vertices(vh, std::back_inserter(vh_neighbor));

            std::vector<int> RT_neighbors_tmp;

            for (auto nb : vh_neighbor)
            {
                int nb_idx = point_to_idx[RPDPoint(nb->point().x(), nb->point().y(), nb->point().z())];
                
                // filter out some cases
                double vh_nb_angle = safetyAcos(points_[vh_idx].nor_.dot(points_[nb_idx].nor_));
                double nb_length = (points_[vh_idx].cor_ - points_[nb_idx].cor_).norm();

                if (vh_nb_angle > connect_angle_constraint_)
                {
                    continue;
                }
                if (is_feature_point(points_[vh_idx]) || is_feature_point(points_[nb_idx]))
                {
                    if (nb_length > connect_feature_distance_ratio_ * radius_)
                    {
                        continue;
                    }
                }
                else
                {
                    if (nb_length > connect_general_distance_ratio_ * radius_)
                    {
                        continue;
                    }
                }

                RT_neighbors_tmp.push_back(nb_idx);
            }
            RT_neighbors[vh_idx] = RT_neighbors_tmp;
        }

#ifdef RPD_DEBUG
        if (1)
        {
            std::map<std::pair<int, int>, bool> edge_map;
            std::ofstream rt_vis(test_file_path + test_func_name + "connection_processed.obj");
            for (int i = 0; i < points_.size(); i++)
            {
                rt_vis << "v " << points_[i].cor_.x() << " " << points_[i].cor_.y() << " " << points_[i].cor_.z() << std::endl;
            }
            for (int i = 0; i < points_.size(); i++)
            {
                for (int j = 0; j < RT_neighbors[i].size(); j++)
                {
                    std::pair<int, int> pp = std::make_pair(std::min(i, RT_neighbors[i][j]), std::max(i, RT_neighbors[i][j]));

                    if (edge_map.find(pp) == edge_map.end())
                    {
                        edge_map[pp] = true;
                        rt_vis << "l " << i + 1 << " " << RT_neighbors[i][j] + 1 << std::endl;
                    }
                }
            }
            rt_vis.close();
        }
#endif

        RT_neighbors_o = RT_neighbors;
        time.stop();

        // record time count
        RecordTimeCount("RegularTriangulation", time.time());

        return 1;
	}

    int RPDRecon::PoissonBasedReconstruction(std::vector<Eigen::Vector3i>& output_faces)
    {
        std::string test_func_name = "PBR_";

        if (radius_ == 0)
        {
            std::vector<K_Point> KP_points;
            for (int i = 0; i < points_.size(); i++)
            {
                K_Point kp(points_[i].cor_.x(), points_[i].cor_.y(), points_[i].cor_.z());
                KP_points.push_back(kp);
            }
            radius_ = CGAL::compute_average_spacing<CGAL::Sequential_tag>(KP_points, 6);
        }

        int numProcs = omp_get_num_procs();
        omp_set_num_threads(2 * numProcs - 1);

        CGAL::Timer time;
        time.start();

        // init
        output_faces.clear();
        std::map<int, RPDCell*> PCs; // each point's RPD cell
        std::map<std::pair<int, int>, bool> RPD_edges_check; // check original RPD's edge being processed
        std::map<int, std::set<RPDPoint>> RPD_cell_connect_points; // point(i) conneted RPD cell's points
        std::map<std::pair<RPDPoint, RPDPoint>, int> RPD_cell_connects; // RPD cell's point belong to point(i)
        std::map<int, std::vector<Eigen::Vector3i>> RPD_cell_triangles; // point(i)'s cell's triangles: point(i)->RPD cell points
        

        // generate poisson model
        std::vector<Eigen::Vector3d> poisson_vertices; // poisson model's vertices
        std::vector<Eigen::Vector3i> poisson_faces; // poisson model's neighbor faces of each poisson vertex
        if (PoissonSurfaceReconstruction(poisson_vertices, poisson_faces) == 0)
        {
            RPD_DEBUG_ONLY(std::cout << "PoissonBasedReconstruction::PoissonSurfaceReconstruction error!" << std::endl;)
            return 0;
        }


        // get poisson model vertices' neighbor faces
        std::vector<std::vector<int>> poisson_neighbor_faces;
        for (int i = 0; i < poisson_vertices.size(); i++)
        {
            poisson_neighbor_faces.push_back(std::vector<int>());
        }
        for (int i = 0; i < poisson_faces.size(); i++)
        {
            poisson_neighbor_faces[poisson_faces[i].x()].push_back(i);
            poisson_neighbor_faces[poisson_faces[i].y()].push_back(i);
            poisson_neighbor_faces[poisson_faces[i].z()].push_back(i);
        }


        // build poisson model vertices' kd tree
        PointCloud<double> poisson_vertices_cloud;
        nanoflann::SearchParams params;
        const double poisson_search_radius = static_cast<double>(pow(5 * radius_, 2));
        // std::cout << poisson_search_radius << std::endl;
        poisson_vertices_cloud.pts.resize(poisson_vertices.size());
        for (int i = 0; i < poisson_vertices.size(); i++)
        {
            poisson_vertices_cloud.pts[i].x = poisson_vertices[i].x();
            poisson_vertices_cloud.pts[i].y = poisson_vertices[i].y();
            poisson_vertices_cloud.pts[i].z = poisson_vertices[i].z();
        }
        my_kd_tree_t poisson_vertices_index(3, poisson_vertices_cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
        poisson_vertices_index.buildIndex();


        // build RPD / regularity triangulation
        std::vector<std::vector<int>> RT_neighbors;
        if (RegularTriangulation(RT_neighbors) == 0)
        {
            RPD_DEBUG_ONLY(std::cout << "PoissonBasedReconstruction::RegularTriangulation error!" << std::endl;)
            return 0;
        }

        // process each point
        for (int i = 0; i < points_.size(); i++) // prepare for openmp
        {
            RPDCell* PC = new RPDCell(points_[i].cor_);
            PCs[i] = PC;
        }
//#pragma omp parallel for
        for (int i = 0; i < points_.size(); i++)
        {
            RPD_DEBUG_ONLY(if (i % 1000 == 0) std::cout << "processed point:" << " " << i << std::endl;)

            RPDCell* PC = PCs[i]; // point(i)'s cell


            // process connected RPD cell
            for (int j = 0; j < RT_neighbors[i].size(); j++)
            {
                int neighbor_idx = RT_neighbors[i][j];

                double w1 = points_[i].weight_;
                double w2 = points_[neighbor_idx].weight_;

                Eigen::Vector3d mid_point;
                if (1)
                {
                    // power diagram
                    double lambda = 0.5 + 0.5 * (w1 * w1 - w2 * w2) / (points_[i].cor_ - points_[neighbor_idx].cor_).squaredNorm();
                    mid_point = lambda * points_[i].cor_ + (1 - lambda) * points_[neighbor_idx].cor_;
                }
                else
                {
                    // voronoi diagram
                    mid_point = 0.5 * (points_[i].cor_ + points_[neighbor_idx].cor_);
                }

                Eigen::Vector3d dir = points_[neighbor_idx].cor_ - points_[i].cor_;

                RPDPlane mid_point_plane(mid_point, dir, neighbor_idx);

                PC->cutted_planes_.push_back(mid_point_plane);
            } // point(i)'s neighbors -> point(j)


            // check each contributory neighbor point
            std::vector<bool> plane_check; // check each cut plane being processed
            for (int j = 0; j < PC->cutted_planes_.size(); j++)
            {
                RPDPlane plane = PC->cutted_planes_[j];
                std::pair<int, int> cur_edge_pair = std::make_pair(std::min(i, plane.opposite_idx()), std::max(i, plane.opposite_idx()));

                if (RPD_edges_check.find(cur_edge_pair) == RPD_edges_check.end())
                {
                    RPD_edges_check[cur_edge_pair] = true;
                    plane_check.push_back(true);
                }
                else
                {
                    plane_check.push_back(false);
                }
            }


            // find poisson model's nearest face
            std::set<int> poisson_nearest_faces; // point(i)'s nearest poisson faces
            if (1)
            {
                double query_pt[3] = { points_[i].cor_.x(), points_[i].cor_.y(), points_[i].cor_.z() };
                std::vector<std::pair<uint32_t, double>> poisson_ret_matches;
                double cur_poisson_search_radius = poisson_search_radius;
                while (1)
                {
                    const size_t poisson_nMatches = poisson_vertices_index.radiusSearch(&query_pt[0], cur_poisson_search_radius, poisson_ret_matches, params);

                    if (poisson_nMatches == 0)
                    {
                        // std::cout << "cur_search_radius:" << " " << cur_poisson_search_radius << std::endl;
                        cur_poisson_search_radius *= 2;
                    }
                    else
                    {
                        for (int j = 0; j < poisson_nMatches; j++)
                        {
                            for (int jj = 0; jj < poisson_neighbor_faces[poisson_ret_matches[j].first].size(); jj++)
                            {
                                poisson_nearest_faces.insert(poisson_neighbor_faces[poisson_ret_matches[j].first][jj]);
                            }
                        }
                        break;
                    }
                }
            }
            /* // knn search poisson nearest face
            int knn_search_num = 6;
            std::vector<uint32_t> poisson_ret_index(knn_search_num);
            std::vector<double> poisson_out_dist_sqr(knn_search_num);
            const size_t poisson_nMatches = poisson_vertices_index.knnSearch(&query_pt[0], knn_search_num, &poisson_ret_index[0], &poisson_out_dist_sqr[0]);
            for (int j = 0; j < poisson_nMatches; j++)
            {
                for (int jj = 0; jj < poisson_neighbor_faces[poisson_ret_index[j]].size(); jj++)
                {
                    poisson_nearest_faces.insert(poisson_neighbor_faces[poisson_ret_index[j]][jj]);
                }
            }*/


            // process: poisson model's nearest face -> cut by RPD cell
            std::vector<std::vector<Eigen::Vector3d>> poisson_cutted_faces; // poisson_nearest_faces -> cutted by RPD planes
            std::vector<int> poisson_cutted_faces_alive; // poisson cutted faces being processed, 1:init, 0: processed
            std::map<RPDPoint, int> poisson_vertices_type; // 0: inside RPD cell, 1: outside RPD cell, -1: cutted point, -2: checked cutted point
            for (int pfi : poisson_nearest_faces)
            {
                bool pfi_cutted = false;

                poisson_vertices_type[RPDPoint(poisson_vertices[poisson_faces[pfi].x()])] = 0;
                poisson_vertices_type[RPDPoint(poisson_vertices[poisson_faces[pfi].y()])] = 0;
                poisson_vertices_type[RPDPoint(poisson_vertices[poisson_faces[pfi].z()])] = 0;

                for (int pi = 0; pi < PC->cutted_planes_.size(); pi++)
                {
                    if (plane_check[pi] == false)
                    {
                        continue;
                    }

                    RPDPlane cur_plane = PC->cutted_planes_[pi];

                    if (cur_plane.is_on_positive_side(poisson_vertices[poisson_faces[pfi].x()]))
                    {
                        poisson_vertices_type[RPDPoint(poisson_vertices[poisson_faces[pfi].x()])] = 1;
                    }
                    if (cur_plane.is_on_positive_side(poisson_vertices[poisson_faces[pfi].y()]))
                    {
                        poisson_vertices_type[RPDPoint(poisson_vertices[poisson_faces[pfi].y()])] = 1;
                    }
                    if (cur_plane.is_on_positive_side(poisson_vertices[poisson_faces[pfi].z()]))
                    {
                        poisson_vertices_type[RPDPoint(poisson_vertices[poisson_faces[pfi].z()])] = 1;
                    }

                    double f1, f2, f3;
                    f1 = cur_plane.signed_distance_to_plane(poisson_vertices[poisson_faces[pfi].x()]);
                    f2 = cur_plane.signed_distance_to_plane(poisson_vertices[poisson_faces[pfi].y()]);
                    f3 = cur_plane.signed_distance_to_plane(poisson_vertices[poisson_faces[pfi].z()]);
                    if ((f1 * f2 < 0) || (f2 * f3 < 0) || (f3 * f1 < 0))
                    {
                        pfi_cutted = true;
                    }
                }

                if (pfi_cutted)
                {
                    Eigen::Vector3i cur_p_face = poisson_faces[pfi];
                    std::vector<Eigen::Vector3d> poisson_cutted_faces_tmp;
                    poisson_cutted_faces_tmp.push_back(poisson_vertices[cur_p_face.x()]);
                    poisson_cutted_faces_tmp.push_back(poisson_vertices[cur_p_face.y()]);
                    poisson_cutted_faces_tmp.push_back(poisson_vertices[cur_p_face.z()]);

                    poisson_cutted_faces.push_back(poisson_cutted_faces_tmp);
                    poisson_cutted_faces_alive.push_back(1);
                }
            }


            // process poisson cutted faces
            for (int j = 0; j < poisson_cutted_faces.size(); j++)
            {
                RPD_DEBUG_ONLY(if (poisson_cutted_faces.size() % 1000  < 2) std::cout << "cur poisson_cutted_faces size:" << " " << poisson_cutted_faces.size() << std::endl;)

                std::vector<Eigen::Vector3d> cur_f = poisson_cutted_faces[j]; // cur poisson_cutted_face

                if (poisson_cutted_faces_alive[j] != 1)
                {
                    continue;
                }

                // cut poisson face(j)
                std::vector<Eigen::Vector3d> new_face;
                bool found_cut_face = false;
                int fdp = 0;
                for (int pi = 0; pi < PC->cutted_planes_.size(); pi++)
                {
                    RPDPlane cur_plane = PC->cutted_planes_[pi];

                    if (found_cut_face)
                    {
                        break;
                    }

                    std::vector<Eigen::Vector3d> new_face_tmp;

                    for (int vj = 0; vj < cur_f.size(); vj++)
                    {
                        Eigen::Vector3d P1, P2;
                        P1 = cur_f[vj];
                        if (vj == cur_f.size() - 1)
                        {
                            P2 = cur_f[0];
                        }
                        else
                        {
                            P2 = cur_f[vj + 1];
                        }

                        double f1, f2;
                        f1 = cur_plane.signed_distance_to_plane(P1);
                        f2 = cur_plane.signed_distance_to_plane(P2);
                        if (fabs(f1) < EPS)
                        {
                            new_face_tmp.push_back(P1);
                            continue;
                        }
                        if (fabs(f2) < EPS)
                        {
                            new_face_tmp.push_back(P1);
                            continue;
                        }
                        if (f1 * f2 < 0)
                        {
                            Eigen::Vector3d new_point;
                            f1 = fabs(f1);
                            f2 = fabs(f2);
                            new_point = P1 + (P2 - P1) * (f1 / (f1 + f2)); /// cutted point

                            new_face_tmp.push_back(P1);
                            new_face_tmp.push_back(new_point);

                            poisson_vertices_type[RPDPoint(new_point)] = -1;

                            found_cut_face = true;
                        }
                        else
                        {
                            new_face_tmp.push_back(P1);
                        }
                    }

                    if (found_cut_face)
                    {
                        new_face = new_face_tmp;
                        fdp = cur_plane.opposite_idx();
                    }
                }

                // whether is found cut face
                if (found_cut_face == false)
                {
                    bool is_point_inside = false; // if all points of this poisson face is inside the RPD cell
                    double max_distance_o = -MAX;

                    for (int vj = 0; vj < cur_f.size(); vj++)
                    {
                        if (poisson_vertices_type[RPDPoint(cur_f[vj])] == 0)
                        {
                            is_point_inside = true;
                        }

                        max_distance_o = std::max(max_distance_o, PC->get_max_distance(cur_f[vj]));
                    }
                    if (max_distance_o < EPS) // all points are on the cut planes
                    {
                        is_point_inside = true;
                    }


                    if (is_point_inside)
                    {
                        std::vector<bool> is_on_cutted_planes; // is faces' points on the RPD cell
                        std::vector<Eigen::Vector3d> on_cutted_planes_point;

                        for (int vj = 0; vj < cur_f.size(); vj++)
                        {
                            Eigen::Vector3d P1 = cur_f[vj];
                            double max_distance = PC->get_max_distance(P1);

                            if (fabs(max_distance) < EPS)
                            {
                                is_on_cutted_planes.push_back(true);
                                on_cutted_planes_point.push_back(P1);
                            }
                            else
                            {
                                is_on_cutted_planes.push_back(false);
                            }
                        }
                        if (on_cutted_planes_point.size() < 2)
                        {
                            continue;
                        }
                        // 2 or 3 points are on the RPD cell

                        std::vector<Eigen::Vector3d> drawed_points;
                        for (int vj = 0; vj < cur_f.size(); vj++)
                        {
                            Eigen::Vector3d P1, P2, Pmid;
                            bool V1, V2; // is P1, P2 on the RPD cell
                            P1 = cur_f[vj];
                            V1 = is_on_cutted_planes[vj];
                            if (vj == cur_f.size() - 1)
                            {
                                P2 = cur_f[0];
                                V2 = is_on_cutted_planes[0];
                            }
                            else
                            {
                                P2 = cur_f[vj + 1];
                                V2 = is_on_cutted_planes[vj + 1];
                            }

                            bool double_check = false; // check is all of the edge in on the RPD cell
                            Pmid = P1 + ((P2 - P1) / 2);
                            
                            double max_distance = PC->get_max_distance(Pmid);

                            if (fabs(max_distance) > EPS)
                            {
                                double_check = true;
                            }

                            if (V1 && V2 && double_check == false)
                            {
                                // all of the edge is on the RPD cell
                                drawed_points.push_back(P1);
                                drawed_points.push_back(P2);

                                RPDPoint RPDP1(P1);
                                RPDPoint RPDP2(P2);
                                RPD_cell_connect_points[i].insert(RPDP1);
                                RPD_cell_connect_points[i].insert(RPDP2);

                                if (RPDP1 < RPDP2)
                                {
                                    RPD_cell_connects[std::make_pair(RPDP1, RPDP2)] = i;
                                }
                                else
                                {
                                    RPD_cell_connects[std::make_pair(RPDP2, RPDP1)] = i;
                                }
                                /*
                                if (RPDP1 < RPDP2)
                                {
                                    if (RPD_cell_connects.find(std::make_pair(RPDP1, RPDP2)) == RPD_cell_connects.end())
                                    {
                                        RPD_cell_connects[std::make_pair(RPDP1, RPDP2)] = i;
                                    }
                                }
                                else
                                {
                                    if (RPD_cell_connects.find(std::make_pair(RPDP2, RPDP1)) == RPD_cell_connects.end())
                                    {
                                        RPD_cell_connects[std::make_pair(RPDP2, RPDP1)] = i;
                                    }
                                }*/
                            }
                        }
                    } // is_point_inside == true
                    else
                    {
                        poisson_cutted_faces_alive[j] = 0;
                    }
                } // found_cut_face == 0
                else
                {
                    std::vector<Eigen::Vector3d> new_face1, new_face2; // cut new_face to 2 face
                    bool fdd = false;

                    for (int k = 0; k < new_face.size(); k++)
                    {
                        if (poisson_vertices_type[RPDPoint(new_face[k])] == -1)
                        {
                            poisson_vertices_type[RPDPoint(new_face[k])] = -2;

                            if (fdd == false)
                            {
                                new_face1.push_back(new_face[k]);
                                new_face2.push_back(new_face[k]);
                                fdd = true;
                                continue;
                            }
                            else
                            {
                                new_face1.push_back(new_face[k]);
                                new_face2.push_back(new_face[k]);
                                fdd = false;
                                continue;
                            }
                        }

                        if (fdd == false)
                        {
                            new_face1.push_back(new_face[k]);
                        }
                        else
                        {
                            new_face2.push_back(new_face[k]);
                        }
                    }

                    poisson_cutted_faces_alive[j] = 0;
                    poisson_cutted_faces_alive.push_back(1);
                    poisson_cutted_faces_alive.push_back(1);
                    poisson_cutted_faces.push_back(new_face1);
                    poisson_cutted_faces.push_back(new_face2);

                    // std::cout << poisson_cutted_faces.size() << std::endl;
                } // found_cut_face != 0
            } // poisson cutted faces(j)


            // convert RPD cell to faces
            RPD_cell_triangles[i] = std::vector<Eigen::Vector3i>();
            for (const auto RPD_p : RPD_cell_connect_points[i])
            {
                for (int pi = 0; pi < PC->cutted_planes_.size(); pi++)
                {
                    RPDPlane plane1 = PC->cutted_planes_[pi];
                    double f1 = plane1.signed_distance_to_plane(RPD_p.cor_);

                    for (int pj = pi + 1; pj < PC->cutted_planes_.size(); pj++)
                    {
                        RPDPlane plane2 = PC->cutted_planes_[pj];
                        double f2 = plane2.signed_distance_to_plane(RPD_p.cor_);

                        if (fabs(f1) < EPS && fabs(f2) < EPS)
                        {
                            int a, b, c;
                            /*
                            a = std::min(i, std::min(plane1.opposite_idx_, plane2.opposite_idx_));
                            c = std::max(i, std::max(plane1.opposite_idx_, plane2.opposite_idx_));
                            if (i != a && i != c)
                            {
                                b = i;
                            }
                            else if (plane1.opposite_idx_ != a && plane1.opposite_idx_ != c)
                            {
                                b = plane1.opposite_idx_;
                            }
                            else
                            {
                                b = plane2.opposite_idx_;
                            }*/

                            a = i;
                            b = plane1.opposite_idx();
                            c = plane2.opposite_idx();

                            RPD_cell_triangles[i].push_back(Eigen::Vector3i(a, b, c));
                        }
                    }
                }
            }
        } // point(i)


        // check if inside other RPD
        std::map<int, bool> is_inside_other_RPD;
        for (int i = 0; i < points_.size(); i++)
        {
            is_inside_other_RPD[i] = false;
            RPDCell* PC = PCs[i];

            bool flag_RPD_p = false;
            for (const RPDPoint RPD_p : RPD_cell_connect_points[i])
            {
                bool flag_j = false;
                for (int j = 0; j < RT_neighbors[i].size(); j++)
                {
                    int n_idx = RT_neighbors[i][j];

                    if (PCs.find(n_idx) == PCs.end())
                    {
                        continue;
                    }
                    if (PCs[n_idx]->cutted_planes_.empty())
                    {
                        continue;
                    }

                    flag_j = false;
                    for (int pj = 0; pj < PCs[n_idx]->cutted_planes_.size(); pj++)
                    {
                        RPDPlane plane = PCs[n_idx]->cutted_planes_[pj];

                        double f1 = plane.signed_distance_to_plane(RPD_p.cor_);

                        if (f1 >= 0 || fabs(f1) < EPS)
                        {
                            flag_j = true;
                            break;
                        }
                    }

                    if (!flag_j)
                    {
                        break;
                    }
                }
                if (flag_j)
                {
                    flag_RPD_p = true;
                    break;
                }
            }

            if (!flag_RPD_p)
            {
                is_inside_other_RPD[i] = true;
            }
        }

#ifdef RPD_DEBUG
        if (1)
        {
            std::ofstream RPD_out(test_file_path + test_func_name + "restricted_power_diagram.obj");
            std::map<RPDPoint, int> RPDPoint_to_idx;
            int pn = 0;
            for (const auto e : RPD_cell_connects)
            {
                RPDPoint P1 = e.first.first;
                RPDPoint P2 = e.first.second;

                if (RPDPoint_to_idx.find(P1) == RPDPoint_to_idx.end())
                {
                    RPDPoint_to_idx[P1] = pn;
                    RPD_out << "v" << " " << P1.cor_.x() << " " << P1.cor_.y() << " " << P1.cor_.z() << std::endl;
                    pn++;
                }
                if (RPDPoint_to_idx.find(P2) == RPDPoint_to_idx.end())
                {
                    RPDPoint_to_idx[P2] = pn;
                    RPD_out << "v" << " " << P2.cor_.x() << " " << P2.cor_.y() << " " << P2.cor_.z() << std::endl;
                    pn++;
                }

                RPD_out << "l" << " " << RPDPoint_to_idx[P1] + 1 << " " << RPDPoint_to_idx[P2] + 1 << std::endl;
            }
            RPD_out.close();
        }
#endif


        // fix connection
        std::map<std::pair<int, int>, int> degree_of_edge;
        for (const auto cell : RPD_cell_triangles)
        {
            if (is_inside_other_RPD[cell.first] == 1)
            {
                continue;
            }

            for (int i = 0; i < cell.second.size(); i++)
            {
                int a = cell.second[i].x();
                int b = cell.second[i].y();
                int c = cell.second[i].z();

                std::pair<int, int> edge_ab(std::min(a, b), std::max(a, b));
                std::pair<int, int> edge_bc(std::min(b, c), std::max(b, c));
                std::pair<int, int> edge_ca(std::min(c, a), std::max(c, a));

                if (degree_of_edge.find(edge_ab) == degree_of_edge.end())
                {
                    degree_of_edge[edge_ab] = 0;
                }
                if (degree_of_edge.find(edge_bc) == degree_of_edge.end())
                {
                    degree_of_edge[edge_bc] = 0;
                }
                if (degree_of_edge.find(edge_ca) == degree_of_edge.end())
                {
                    degree_of_edge[edge_ca] = 0;
                }

                degree_of_edge[edge_ab]++;
                degree_of_edge[edge_bc]++;
                degree_of_edge[edge_ca]++;
            }
        }
        bool fix_done = true;
        while (fix_done)
        {
            fix_done = false;

            for (auto cell : RPD_cell_triangles)
            {
                for (int i = 0; i < cell.second.size(); i++)
                {
                    int a = cell.second[i].x();
                    int b = cell.second[i].y();
                    int c = cell.second[i].z();

                    if (a == -1 || b == -1 || c == -1)
                    {
                        continue;
                    }

                    std::pair<int, int> edge_ab(std::min(a, b), std::max(a, b));
                    std::pair<int, int> edge_bc(std::min(b, c), std::max(b, c));
                    std::pair<int, int> edge_ca(std::min(c, a), std::max(c, a));

                    if (degree_of_edge[edge_ab] < 2 || degree_of_edge[edge_bc] < 2 || degree_of_edge[edge_ca] < 2)
                    {
                        if (degree_of_edge[edge_ab] > 2 || degree_of_edge[edge_bc] > 2 || degree_of_edge[edge_ca] > 2)
                        {
                            degree_of_edge[edge_ab]--;
                            degree_of_edge[edge_bc]--;
                            degree_of_edge[edge_ca]--;

                            cell.second[i] = Eigen::Vector3i(-1, -1, -1);
                            fix_done = true;
                        }
                    }
                }
            }
        }


        // convert to connection
        std::vector<Eigen::Vector3i> output_faces_tmp;
        for (const auto cell : RPD_cell_triangles)
        {
            if (is_inside_other_RPD[cell.first] == 1)
            {
                continue;
            }
            for (int i = 0; i < cell.second.size(); i++)
            {
                int a = cell.second[i].x();
                int b = cell.second[i].y();
                int c = cell.second[i].z();

                if (a == -1 || b == -1 || c == -1)
                {
                    continue;
                }

                output_faces_tmp.push_back(Eigen::Vector3i(a, b, c));
            }
        }


        // orient face normals
#pragma omp parallel for
        for (int i = 0; i < output_faces_tmp.size(); i++)
        {
            RPDPoint P1 = points_[output_faces_tmp[i].x()];
            RPDPoint P2 = points_[output_faces_tmp[i].y()];
            RPDPoint P3 = points_[output_faces_tmp[i].z()];

            Eigen::Vector3d face_normal = (P2.cor_ - P1.cor_).cross(P3.cor_ - P2.cor_);
            
            double dot_check = P1.nor_.dot(face_normal) + P2.nor_.dot(face_normal) + P3.nor_.dot(face_normal);

            if (dot_check < 0)
            {
                output_faces_tmp[i] = Eigen::Vector3i(output_faces_tmp[i].x(), output_faces_tmp[i].z(), output_faces_tmp[i].y());
            }
        }

        output_faces = output_faces_tmp;

#ifdef RPD_DEBUG
        if (1)
        {
            std::ofstream RPDrecon_out(test_file_path + test_func_name + "reconstruction.obj");
            for (int i = 0; i < points_.size(); i++)
            {
                RPDrecon_out << "v" << " " << points_[i].cor_.x() << " " << points_[i].cor_.y() << " " << points_[i].cor_.z() << std::endl;
            }
            for (int i = 0; i < output_faces.size(); i++)
            {
                RPDrecon_out << "f" << " " << output_faces[i].x() + 1 << " " << output_faces[i].y() + 1 << " " << output_faces[i].z() + 1 << std::endl;
            }
            RPDrecon_out.close();
        }
#endif

        time.stop();

        // record time count
        RecordTimeCount("PoissonBasedReconstruction", time.time());

        return 1;
    }

    int RPDRecon::ParallelDiskReconstruction(std::vector<Eigen::Vector3i>& output_faces)
    {
        std::string test_func_name = "PDR_";

        int numProcs = omp_get_num_procs();
        omp_set_num_threads(2 * numProcs - 1);

        CGAL::Timer time;
        time.start();

        output_faces.clear();

        // init
        std::vector<RPDDisk*> disks;
        double init_disk_radius = 5 * radius_;


        // get points' weight range
        double max_weight = -MAX;
        double min_weight = MAX;
        for (int i = 0; i < points_.size(); i++)
        {
            max_weight = std::max(max_weight, points_[i].weight_);
            min_weight = std::min(min_weight, points_[i].weight_);
        }


        // build points' kd tree
        PointCloud<double> point_cloud;
        point_cloud.pts.resize(points_.size());
        for (int i = 0; i < points_.size(); i++)
        {
            point_cloud.pts[i].x = points_[i].cor_.x();
            point_cloud.pts[i].y = points_[i].cor_.y();
            point_cloud.pts[i].z = points_[i].cor_.z();
        }
        my_kd_tree_t point_cloud_index(3, point_cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
        point_cloud_index.buildIndex();
        nanoflann::SearchParams params;
        // const double search_radius = static_cast<double>(std::pow(sqrt(init_disk_radius * init_disk_radius + (max_weight * max_weight - min_weight * min_weight) * radius_ * radius_) + init_disk_radius, 2));
        const double search_radius = static_cast<double>(std::pow(2 * init_disk_radius, 2));


        // process each point's disk
        std::map<int, std::map<int, bool>> RPD_neighbors;
        for (int i = 0; i < points_.size(); i++) // prepare for openmp
        {
            RPDDisk* disk = new RPDDisk(points_[i].cor_, points_[i].nor_, 8, init_disk_radius);
            disks.push_back(disk);
        }
#pragma omp parallel for
        for (int i = 0; i < points_.size(); i++)
        {
            RPDDisk* disk = disks[i];


            // find nearest points
            double query_pt[3] = { points_[i].cor_.x(), points_[i].cor_.y(), points_[i].cor_.z() };
            std::vector<std::pair<uint32_t, double>> ret_matches;
            const size_t nMatches = point_cloud_index.radiusSearch(&query_pt[0], search_radius, ret_matches, params);


            // cut the disk
            for (int j = 0; j < nMatches; j++)
            {
                if (ret_matches[j].second < EPS)
                {
                    continue;
                }

                Eigen::Vector3d mid_point = 0.5 * (points_[i].cor_ + points_[ret_matches[j].first].cor_);
                Eigen::Vector3d dir = points_[i].cor_ - points_[ret_matches[j].first].cor_;
                RPDPlane plane(mid_point, dir, ret_matches[j].first);

                disk->cut_by_plane(plane);
            }
            

            // convert to map connection
            std::map<int, bool> disk_map;
            for (int j = 0; j < disk->border_planes_.size(); j++)
            {
                RPDPlane cur_plane = disk->border_planes_[j];

                if (cur_plane.opposite_idx() == -1)
                {
                    continue;
                }
                disk_map[cur_plane.opposite_idx()] = true;
            }
            RPD_neighbors[i] = disk_map;
        }


        // select triangles
        std::vector<Eigen::Vector3i> good_triangles;
        std::vector<Eigen::Vector3i> not_bad_triangles;
        std::vector<Eigen::Vector3i> bad_triangles;
        std::map<std::tuple<int, int, int>, bool> triangle_map;
        for (int i = 0; i < points_.size(); i++)
        {
            RPDDisk* disk = disks[i];

            std::vector<bool> is_good_plane;
            for (int j = 0; j < disk->border_planes_.size(); j++)
            {
                int opposite_idx = disk->border_planes_[j].opposite_idx();

                if (opposite_idx == -1)
                {
                    is_good_plane.push_back(false);
                }
                else
                {
                    if (RPD_neighbors[opposite_idx].find(i) == RPD_neighbors[opposite_idx].end())
                    {
                        is_good_plane.push_back(false);
                    }
                    else
                    {
                        is_good_plane.push_back(true);
                    }
                }
            }

            for (int j = 0; j < disk->border_planes_.size(); j++)
            {
                int next_j = j + 1;
                if (next_j == disk->border_planes_.size())
                {
                    next_j = 0;
                }

                int a = i;
                int b = disk->border_planes_[j].opposite_idx();
                int c = disk->border_planes_[next_j].opposite_idx();

                if (b == -1 || c == -1)
                {
                    continue;
                }

                int aa = std::min({ a, b, c });
                int bb;
                int cc = std::max({ a, b, c });
                if (a != aa && a != cc)
                {
                    bb = a;
                }
                if (b != aa && b != cc)
                {
                    bb = b;
                }
                if (c != aa && c != cc)
                {
                    bb = c;
                }
                std::tuple<int, int, int> tri_tuple(aa, bb, cc);
                if (triangle_map.find(tri_tuple) == triangle_map.end())
                {
                    triangle_map[tri_tuple] = true;
                }
                else
                {
                    continue;
                }
                

                if (is_good_plane[j] && is_good_plane[next_j])
                {
                    good_triangles.push_back(Eigen::Vector3i(a, b, c));
                }
                else if (is_good_plane[j] || is_good_plane[next_j])
                {
                    not_bad_triangles.push_back(Eigen::Vector3i(a, b, c));
                }
                else
                {
                    bad_triangles.push_back(Eigen::Vector3i(a, b, c));
                }
            }
        }


        // merge triangles
        std::vector<Eigen::Vector3i> output_faces_tmp;
        std::map<std::pair<int, int>, int> edge_cnt;
        std::vector<int> good_triangles_features;
        for (int i = 0; i < good_triangles.size(); i++)
        {
            int p1 = good_triangles[i].x();
            int p2 = good_triangles[i].y();
            int p3 = good_triangles[i].z();

            int feature_cnt = 0;
            if (is_feature_point(points_[p1]))
            {
                feature_cnt++;
            }
            if (is_feature_point(points_[p2]))
            {
                feature_cnt++;
            }
            if (is_feature_point(points_[p3]))
            {
                feature_cnt++;
            }

            good_triangles_features.push_back(feature_cnt);
        }
        for (int fn = 3; fn >= 0; fn--)
        {
            for (int i = 0; i < good_triangles.size(); i++)
            {
                if (good_triangles_features[i] == fn)
                {
                    int p1 = good_triangles[i].x();
                    int p2 = good_triangles[i].y();
                    int p3 = good_triangles[i].z();

                    std::pair<int, int> e12(std::min(p1, p2), std::max(p1, p2));
                    std::pair<int, int> e23(std::min(p2, p3), std::max(p2, p3));
                    std::pair<int, int> e31(std::min(p3, p1), std::max(p3, p1));

                    if (edge_cnt.find(e12) == edge_cnt.end())
                    {
                        edge_cnt[e12] = 0;
                    }
                    if (edge_cnt.find(e23) == edge_cnt.end())
                    {
                        edge_cnt[e23] = 0;
                    }
                    if (edge_cnt.find(e31) == edge_cnt.end())
                    {
                        edge_cnt[e31] = 0;
                    }

                    if (edge_cnt[e12] > 1 || edge_cnt[e23] > 1 || edge_cnt[e31] > 1)
                    {
                        continue;
                    }

                    output_faces_tmp.push_back(good_triangles[i]);

                    edge_cnt[e12]++;
                    edge_cnt[e23]++;
                    edge_cnt[e31]++;
                }
            }
        }
        for (int i = 0; i < not_bad_triangles.size(); i++)
        {
            int p1 = not_bad_triangles[i].x();
            int p2 = not_bad_triangles[i].y();
            int p3 = not_bad_triangles[i].z();

            std::pair<int, int> e12(std::min(p1, p2), std::max(p1, p2));
            std::pair<int, int> e23(std::min(p2, p3), std::max(p2, p3));
            std::pair<int, int> e31(std::min(p3, p1), std::max(p3, p1));

            if (edge_cnt.find(e12) == edge_cnt.end())
            {
                edge_cnt[e12] = 0;
            }
            if (edge_cnt.find(e23) == edge_cnt.end())
            {
                edge_cnt[e23] = 0;
            }
            if (edge_cnt.find(e31) == edge_cnt.end())
            {
                edge_cnt[e31] = 0;
            }

            if (edge_cnt[e12] > 1 || edge_cnt[e23] > 1 || edge_cnt[e31] > 1)
            {
                continue;
            }

            output_faces_tmp.push_back(not_bad_triangles[i]);

            edge_cnt[e12]++;
            edge_cnt[e23]++;
            edge_cnt[e31]++;
        }
        

#ifdef RPD_DEBUG
        if (1)
        {
            std::ofstream disk_out(test_file_path + test_func_name + "disks.obj");
            std::map<RPDPoint, int> rp_to_idx;
            int rp_point_cnt = points_.size();
            for (int i = 0; i < points_.size(); i++)
            {
                disk_out << "v" << " " << points_[i].cor_.transpose() << " " << "255 0 0 " << std::endl;
            }
            for (int i = 0; i < points_.size(); i++)
            {
                for (int j = 0; j < disks[i]->border_points_.size(); j++)
                {
                    int next_j = j + 1;
                    if (next_j == disks[i]->border_points_.size())
                    {
                        next_j = 0;
                    }

                    RPDPoint P1 = disks[i]->border_points_[j];
                    RPDPoint P2 = disks[i]->border_points_[next_j];
                    if (rp_to_idx.find(P1) == rp_to_idx.end())
                    {
                        disk_out << "v" << " " << P1.cor_.transpose() << " " << "255 255 255 " << std::endl;
                        rp_to_idx[P1] = rp_point_cnt;
                        rp_point_cnt++;
                    }
                    if (rp_to_idx.find(P2) == rp_to_idx.end())
                    {
                        disk_out << "v" << " " << P2.cor_.transpose() << " " << "255 255 255 " << std::endl;
                        rp_to_idx[P2] = rp_point_cnt;
                        rp_point_cnt++;
                    }

                    disk_out << "f" << " " << i + 1 << " " << rp_to_idx[P1] + 1 << " " << rp_to_idx[P2] + 1 << std::endl;
                }
            }
            disk_out.close();
        }
#endif

#ifdef RPD_DEBUG
        if (1)
        {
            std::ofstream disk_out(test_file_path + test_func_name + "disks_connection.obj");
            for (int i = 0; i < points_.size(); i++)
            {
                disk_out << "v" << " " << points_[i].cor_.transpose() << std::endl;
            }
            std::map<std::pair<int, int>, bool> pair_check;
            for (int i = 0; i < points_.size(); i++)
            {
                for (const auto mp : RPD_neighbors[i])
                {
                    if (pair_check.find(std::make_pair(i, mp.first)) == pair_check.end())
                    {
                        pair_check[std::make_pair(i, mp.first)] = true;
                        disk_out << "l" << " " << i + 1 << " " << mp.first + 1 << std::endl;
                    }
                }
            }
            disk_out.close();
        }
#endif

#ifdef RPD_DEBUG
        if (1)
        {
            std::ofstream good_out(test_file_path + test_func_name + "triangle_good.obj");
            for (int i = 0; i < points_.size(); i++)
            {
                good_out << "v" << " " << points_[i].cor_.transpose() << std::endl;
            }
            for (int i = 0; i < good_triangles.size(); i++)
            {
                good_out << "f" << " " << good_triangles[i].x() + 1 << " " << good_triangles[i].y() + 1 << " " << good_triangles[i].z() + 1 << std::endl;
            }
            good_out.close();

            std::ofstream not_bad_out(test_file_path + test_func_name + "triangle_not_bad.obj");
            for (int i = 0; i < points_.size(); i++)
            {
                not_bad_out << "v" << " " << points_[i].cor_.transpose() << std::endl;
            }
            for (int i = 0; i < not_bad_triangles.size(); i++)
            {
                not_bad_out << "f" << " " << not_bad_triangles[i].x() + 1 << " " << not_bad_triangles[i].y() + 1 << " " << not_bad_triangles[i].z() + 1 << std::endl;
            }
            not_bad_out.close();

            std::ofstream bad_out(test_file_path + test_func_name + "triangle_bad.obj");
            for (int i = 0; i < points_.size(); i++)
            {
                bad_out << "v" << " " << points_[i].cor_.transpose() << std::endl;
            }
            for (int i = 0; i < bad_triangles.size(); i++)
            {
                bad_out << "f" << " " << bad_triangles[i].x() + 1 << " " << bad_triangles[i].y() + 1 << " " << bad_triangles[i].z() + 1 << std::endl;
            }
            bad_out.close();
        }
#endif

#ifdef RPD_DEBUG
        if (1)
        {
            std::ofstream out(test_file_path + test_func_name + "reconstruction.obj");
            for (int i = 0; i < points_.size(); i++)
            {
                out << "v" << " " << points_[i].cor_.transpose() << std::endl;
            }
            for (int i = 0; i < output_faces_tmp.size(); i++)
            {
                out << "f" << " " << output_faces_tmp[i].x() + 1 << " " << output_faces_tmp[i].y() + 1 << " " << output_faces_tmp[i].z() + 1 << std::endl;
            }
            out.close();
        }
#endif

        time.stop();

        RecordTimeCount("ParallelDiskReconstruction", time.time());

        return 1;
    }

    //others
    void RPDRecon::RecordTimeCount(std::string func_name,
        double time)
    {
        if (time_count_.find(func_name) == time_count_.end())
        {
            time_count_[func_name] = time;
        }
        else
        {
            time_count_[func_name] += time;
        }
    }

    void RPDRecon::PrintTimeCount()
    {
        std::cout << "---TIME COUNT---" << std::endl;
        for (auto tc : time_count_)
        {
            std::cout << tc.first << ":" << "\t" << tc.second << ".sec" << std::endl;
        }
    }

    void RPDRecon::output_feature_points(std::string filename,
        double r,
        double g,
        double b)
    {
        r = std::min(std::max(r, double(0)), double(255));
        g = std::min(std::max(g, double(0)), double(255));
        b = std::min(std::max(b, double(0)), double(255));

        Eigen::Vector3d color(r, g, b);
        double min_weight = MAX;
        double max_weight = -MAX;
        for (int i = 0; i < points_.size(); i++)
        {
            min_weight = std::min(min_weight, points_[i].weight_);
            max_weight = std::max(max_weight, points_[i].weight_);
        }
        if (abs(max_weight - min_weight) > EPS)
        {
            std::ofstream out(filename);
            for (int i = 0; i < points_.size(); i++)
            {
                out << "v" << " " << points_[i].cor_.transpose() << " " << (color * (points_[i].weight_ - min_weight) / (max_weight - min_weight)).transpose() << std::endl;
            }
            out.close();
        }
    }
} //namespace RPD
