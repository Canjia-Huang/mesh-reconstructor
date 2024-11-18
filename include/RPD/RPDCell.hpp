#ifndef MESH_RECONSTRUCTOR_RPD_RPDCELL_HPP_
#define MESH_RECONSTRUCTOR_RPD_RPDCELL_HPP_

#include <Export.h>
#include <Point.hpp>
#include <RPD/RPDPlane.hpp>
#include <vector>
#include <list>
#include <unordered_set>
// Eigen
#include <Eigen/Dense>

using namespace MeshReconstructor;

namespace MeshReconstructorRPD {
	struct pair_hash {
		template <class T1, class T2>
		size_t operator () (std::pair<T1, T2> const& pair) const {
			size_t h1 = std::hash<T1>()(pair.first);
			size_t h2 = std::hash<T2>()(pair.second);
			return h1 ^ h2;
		}
	};

	class RPDCell {
	public:
		RPDCell() {
			cor_ = Eigen::Vector3d(0, 0, 0);
		}

		~RPDCell() {
			cor_ = Eigen::Vector3d(0, 0, 0);
			(std::vector<RPDPlane>()).swap(cutted_planes_);
			(std::vector<std::pair<Eigen::Vector3i, Point> >()).swap(cutted_vertices_);
		}

		Eigen::Vector3d cor() const { return cor_; }
		std::vector<RPDPlane>& cutted_planes() { return cutted_planes_; }
		std::vector<std::pair<Eigen::Vector3i, Point> >& cutted_vertices() { return cutted_vertices_; }

		void init_cube(
			Eigen::Vector3d& cor, 
			double r
		) {
			cutted_planes_.emplace_back(RPDPlane(cor.x() + r, cor.y(), cor.z(), 1, 0, 0));
			cutted_planes_.emplace_back(RPDPlane(cor.x() - r, cor.y(), cor.z(), -1, 0, 0));
			cutted_planes_.emplace_back(RPDPlane(cor.x(), cor.y() + r, cor.z(), 0, 1, 0));
			cutted_planes_.emplace_back(RPDPlane(cor.x(), cor.y() - r, cor.z(), 0, -1, 0));
			cutted_planes_.emplace_back(RPDPlane(cor.x(), cor.y(), cor.z() + r, 0, 0, 1));
			cutted_planes_.emplace_back(RPDPlane(cor.x(), cor.y(), cor.z() - r, 0, 0, -1));

			cutted_vertices_.emplace_back(std::make_pair(Eigen::Vector3i(0, 2, 4), Point(cor.x() + r, cor.y() + r, cor.z() + r)));
			cutted_vertices_.emplace_back(std::make_pair(Eigen::Vector3i(2, 1, 4), Point(cor.x() - r, cor.y() + r, cor.z() + r)));
			cutted_vertices_.emplace_back(std::make_pair(Eigen::Vector3i(1, 3, 4), Point(cor.x() - r, cor.y() - r, cor.z() + r)));
			cutted_vertices_.emplace_back(std::make_pair(Eigen::Vector3i(3, 0, 4), Point(cor.x() + r, cor.y() - r, cor.z() + r)));
			cutted_vertices_.emplace_back(std::make_pair(Eigen::Vector3i(2, 0, 5), Point(cor.x() + r, cor.y() + r, cor.z() - r)));
			cutted_vertices_.emplace_back(std::make_pair(Eigen::Vector3i(1, 2, 5), Point(cor.x() - r, cor.y() + r, cor.z() - r)));
			cutted_vertices_.emplace_back(std::make_pair(Eigen::Vector3i(3, 1, 5), Point(cor.x() - r, cor.y() - r, cor.z() - r)));
			cutted_vertices_.emplace_back(std::make_pair(Eigen::Vector3i(0, 3, 5), Point(cor.x() + r, cor.y() - r, cor.z() - r)));
		}

		bool cut_by_plane(const RPDPlane& plane) {
			std::list<Eigen::Vector3i> cutted_vertices;
			// find cutted vertices
			for (auto v_it = cutted_vertices_.begin(); v_it != cutted_vertices_.end();) {
				if (plane.is_on_positive_side((*v_it).second.cor())) {
					cutted_vertices.push_back((*v_it).first);
					v_it = cutted_vertices_.erase(v_it);
				}
				else {
					v_it++;
				}
			}
			if (cutted_vertices.size() == 0) {
				return false;
			}

			// find cutted edge -> real cutted vertices
			std::unordered_set<int> remain_planes_set;
			for (auto v_it : cutted_vertices_) {
				remain_planes_set.insert((v_it.first)[0]);
				remain_planes_set.insert((v_it.first)[1]);
				remain_planes_set.insert((v_it.first)[2]);
			}
			std::unordered_set<std::pair<int, int>, pair_hash> cutted_edges;
			for (auto it = cutted_vertices.begin(); it != cutted_vertices.end(); it++) {
				std::pair<int, int> edge1 = std::make_pair((*it)[0], (*it)[1]);
				std::pair<int, int> edge2 = std::make_pair((*it)[1], (*it)[2]);
				std::pair<int, int> edge3 = std::make_pair((*it)[2], (*it)[0]);

				auto f_e1 = cutted_edges.find(std::make_pair((*it)[1], (*it)[0]));
				if (f_e1 == cutted_edges.end()) {
					cutted_edges.insert(edge1);
				}
				else {
					cutted_edges.erase(f_e1);
				}

				auto f_e2 = cutted_edges.find(std::make_pair((*it)[2], (*it)[1]));
				if (f_e2 == cutted_edges.end()) {
					cutted_edges.insert(edge2);
				}
				else {
					cutted_edges.erase(f_e2);
				}

				auto f_e3 = cutted_edges.find(std::make_pair((*it)[0], (*it)[2]));
				if (f_e3 == cutted_edges.end()) {
					cutted_edges.insert(edge3);
				}
				else {
					cutted_edges.erase(f_e3);
				}
			}

			// process
			std::list<int> dangling_planes;
			while (cutted_edges.size() > 0) {
				int cur_cutted_edges_nb = cutted_edges.size();

				for (auto it = cutted_edges.begin(); it != cutted_edges.end();) {
					if (dangling_planes.size() == 0) {
						dangling_planes.push_back((*it).first);
						dangling_planes.push_back((*it).second);
						it = cutted_edges.erase(it);
						continue;
					}

					if ((*it).first != dangling_planes.back() &&
						(*it).second != dangling_planes.front()) {
						it++;
					}
					else {
						if ((*it).first == dangling_planes.back() && (*it).second == dangling_planes.front()) {

						}
						else {
							if ((*it).first == dangling_planes.back()) {
								dangling_planes.push_back((*it).second);
							}
							if ((*it).second == dangling_planes.front()) {
								dangling_planes.push_front((*it).first);
							}
						}

						it = cutted_edges.erase(it);
					}
				}

				if (cutted_edges.size() == cur_cutted_edges_nb) break;
			}

			// renew
			if (dangling_planes.size() > 2) {
				cutted_planes_.emplace_back(plane);

				int planes_cnt = 0;
				std::list<int>::iterator prev_it = dangling_planes.end();
				prev_it--;
				for (std::list<int>::iterator it = dangling_planes.begin(); it != dangling_planes.end(); ++it) {
					RPDPlane* plane1 = &(cutted_planes_[(*prev_it)]);
					RPDPlane* plane2 = &(cutted_planes_[(*it)]);

					Eigen::Vector3d cutted_point = get_tri_planes_cutted_point((*plane1), (*plane2), plane);
					if (std::isnan(cutted_point.x()) || std::isnan(cutted_point.y()) || std::isnan(cutted_point.z()) ||
						std::isinf(cutted_point.x()) || std::isinf(cutted_point.y()) || std::isinf(cutted_point.z())) {

					}
					else {
						cutted_vertices_.emplace_back(
							std::make_pair(
								Eigen::Vector3i((*prev_it), (*it), cutted_planes_.size() - 1),
								Point(cutted_point)));
					}
					prev_it = it;
				}
			}

			return true;
		}

	protected:
		Eigen::Vector3d cor_;
		std::vector<RPDPlane> cutted_planes_;
		std::vector<std::pair<Eigen::Vector3i, Point> > cutted_vertices_;
	};
}

#endif
