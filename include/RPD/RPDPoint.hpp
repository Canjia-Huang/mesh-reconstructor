#ifndef MESH_RECONSTRUCTOR_RPD_RPDPOINT_HPP_
#define MESH_RECONSTRUCTOR_RPD_RPDPOINT_HPP_

#include <Export.h>
#include <Point.hpp>

namespace MeshReconstructorRPD {
	class MESH_RECONSTRUCTOR_EXPORT RPDPoint: public MeshReconstructor::Point {
	public:
		RPDPoint() {
			cor_ = Eigen::Vector3d(0, 0, 0);
			weight_ = 0;
		}

		~RPDPoint() {
			cor_ = Eigen::Vector3d(0, 0, 0);
			weight_ = 0;
		}

		RPDPoint(double x, double y, double z, double w) {
			cor_.x() = x;
			cor_.y() = y;
			cor_.z() = z;
			weight_ = w;
		}

		RPDPoint(Eigen::Vector3d& cor, double w) {
			cor_ = cor;
			weight_ = w;
		}

		double& w() { return weight_; }

		void operator =(const RPDPoint& p) {
			this->cor_ = p.cor_;
			this->weight_ = p.weight_;
		}

		bool operator ==(const RPDPoint& p) const {
			if ((cor_.x() - p.cor_.x()) < MR_EPS && (cor_.x() - p.cor_.x()) > -MR_EPS) {
				if ((cor_.y() - p.cor_.y()) < MR_EPS && (cor_.y() - p.cor_.y()) > -MR_EPS) {
					if ((cor_.z() - p.cor_.z()) < MR_EPS && (cor_.z() - p.cor_.z()) > -MR_EPS) {
						if ((weight_ - p.weight_) < MR_EPS && (weight_ - p.weight_) > -MR_EPS) {
							return true;
						}
					}
				}
			}
			return false;
		}

	protected:
		double weight_;
	};
}

#endif
