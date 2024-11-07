#ifndef MESH_RECONSTRUCTOR_POINT_HPP_
#define MESH_RECONSTRUCTOR_POINT_HPP_

#include <Export.h>
// Eigen
#include <Eigen/Dense>

namespace MeshReconstructor {
	class Point {
	public:
		Point() {
			cor_ = Eigen::Vector3d(0, 0, 0);
		};

		~Point() {
			cor_ = Eigen::Vector3d(0, 0, 0);
		}

		Point(double x, double y, double z) {
			cor_.x() = x;
			cor_.y() = y;
			cor_.z() = z;
		}

		Point(const Eigen::Vector3d& cor) {
			cor_ = cor;
		}

		double& x() { return cor_.x();}

		double& y() { return cor_.y();}

		double& z() { return cor_.z();}

		Eigen::Vector3d& cor() { return cor_; }

		Point operator *(double s) {
			Point rp;
			rp.cor_.x() = s * this->cor_.x();
			rp.cor_.y() = s * this->cor_.y();
			rp.cor_.z() = s * this->cor_.z();
			return rp;
		}

		Point operator +(const Point& p) {
			Point rp;
			rp.cor_.x() = this->cor_.x() + p.cor_.x();
			rp.cor_.y() = this->cor_.y() + p.cor_.y();
			rp.cor_.z() = this->cor_.z() + p.cor_.z();
			return rp;
		}

		Point operator -(const Point& p) {
			Point rp;
			rp.cor_.x() = this->cor_.x() - p.cor_.x();
			rp.cor_.y() = this->cor_.y() - p.cor_.y();
			rp.cor_.z() = this->cor_.z() - p.cor_.z();
			return rp;
		}

		void operator =(const Point& p) {
			this->cor_ = p.cor_;
		}

		bool operator <(const Point& p) const {
			if ((cor_.x() - p.cor_.x()) < MR_EPS && (cor_.x() - p.cor_.x()) > -MR_EPS) {
				if ((cor_.y() - p.cor_.y()) < MR_EPS && (cor_.y() - p.cor_.y()) > -MR_EPS) {
					return (cor_.z() < p.cor_.z());
				}
				return (cor_.y() < p.cor_.y());
			}
			return (cor_.x() < p.cor_.x());
		}

		bool operator ==(const Point& p) const {
			if ((cor_.x() - p.cor_.x()) < MR_EPS && (cor_.x() - p.cor_.x()) > -MR_EPS) {
				if ((cor_.y() - p.cor_.y()) < MR_EPS && (cor_.y() - p.cor_.y()) > -MR_EPS) {
					if ((cor_.z() - p.cor_.z()) < MR_EPS && (cor_.z() - p.cor_.z()) > -MR_EPS) {
						return true;
					}
				}
			}
			return false;
		}

	protected:
		Eigen::Vector3d cor_;
	};
}

#endif
