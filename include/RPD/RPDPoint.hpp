#ifndef MESH_RECONSTRUCTOR_RPD_RPDPOINT_HPP_
#define MESH_RECONSTRUCTOR_RPD_RPDPOINT_HPP_

#include <Export.h>
#include <Point.hpp>

namespace MeshReconstructorRPD {
	class RPDPoint: public MeshReconstructor::Point_with_normal {
	public:
		RPDPoint() {
			cor_ = Eigen::Vector3d(0, 0, 0);
			nor_ = Eigen::Vector3d(0, 0, 0);
			weight_ = 0;
		}

		~RPDPoint() {
			cor_ = Eigen::Vector3d(0, 0, 0);
			nor_ = Eigen::Vector3d(0, 0, 0);
			weight_ = 0;
		}

		RPDPoint(double x, double y, double z, double w) {
			cor_.x() = x;
			cor_.y() = y;
			cor_.z() = z;
			weight_ = w;
		}

		RPDPoint(double x, double y, double z, double nx, double ny, double nz, double w) {
			cor_.x() = x;
			cor_.y() = y;
			cor_.z() = z;
			nor_.x() = nx;
			nor_.y() = ny;
			nor_.z() = nz;
			weight_ = w;
		}

		RPDPoint(Eigen::Vector3d& cor, double w) {
			cor_ = cor;
			weight_ = w;
		}

		RPDPoint(Eigen::Vector3d& cor, Eigen::Vector3d& nor, double w) {
			cor_ = cor;
			nor_ = nor;
			weight_ = w;
		}

		double& w() { return weight_; }

		void operator =(const RPDPoint& p) {
			this->cor_ = p.cor_;
			this->nor_ = p.nor_;
			this->weight_ = p.weight_;
		}

	protected:
		double weight_;
	};
}

#endif
