#ifndef MESH_RECONSTRUCTOR_RPD_RPDPLANE_HPP_
#define MESH_RECONSTRUCTOR_RPD_RPDPLANE_HPP_

#include <Export.h>
// Eigen
#include <Eigen/Dense>

namespace MeshReconstructorRPD {
	class MESH_RECONSTRUCTOR_EXPORT RPDPlane {
	public:
		RPDPlane() {
			a_ = 0;
			b_ = 0;
			c_ = 0;
			d_ = 0;
			opposite_id_ = -1;
		}

		~RPDPlane() {
			a_ = 0;
			b_ = 0;
			c_ = 0;
			d_ = 0;
			opposite_id_ = -1;
		}

		RPDPlane(
			double cx, double cy, double cz,
			double nx, double ny, double nz) {
			a_ = nx;
			b_ = ny;
			c_ = nz;
			d_ = -cx * nx - cy * ny - cz * nz;
		}

		RPDPlane(const Eigen::Vector3d& cor, const Eigen::Vector3d& nor) {
			a_ = nor.x();
			b_ = nor.y();
			c_ = nor.z();
			d_ = -cor.dot(nor);
		}

		double a() const { return a_; }
		double b() const { return b_; }
		double c() const { return c_; }
		double d() const { return d_; }
		int& opposite_id() { return opposite_id_; }

		double sign_dis(
			double x, double y, double z) const {
			return a_ * x + b_ * y + c_ * z + d_;
		}

		double sign_dis(const Eigen::Vector3d& cor) const {
			return a_ * cor.x() + b_ * cor.y() + c_ * cor.z() + d_;
		}

		bool is_on_positive_side(
			double x, double y, double z) const {
			return (sign_dis(x, y, z) > MR_EPS);
		}

		bool is_on_positive_side(const Eigen::Vector3d& cor) const {
			return (sign_dis(cor) > MR_EPS);
		}

	protected:
		double a_;
		double b_;
		double c_;
		double d_;
		int opposite_id_;
	};

	Eigen::Vector3d get_tri_planes_cutted_point(
		const RPDPlane& P1,
		const RPDPlane& P2,
		const RPDPlane& P3
	) {
		double d = P1.a() * (P2.b() * P3.c() - P2.c() * P3.b()) - P1.b() * (P2.a() * P3.c() - P2.c() * P3.a()) + P1.c() * (P2.a() * P3.b() - P2.b() * P3.a());
		double x = -P1.d() * (P2.b() * P3.c() - P2.c() * P3.b()) + P2.d() * (P1.b() * P3.c() - P1.c() * P3.b()) - P3.d() * (P1.b() * P2.c() - P1.c() * P2.b());
		double y = P1.d() * (P2.a() * P3.c() - P2.c() * P3.a()) - P2.d() * (P1.a() * P3.c() - P1.c() * P3.a()) + P3.d() * (P1.a() * P2.c() - P1.c() * P2.a());
		double z = P1.d() * (P2.b() * P3.a() - P2.a() * P3.b()) - P2.d() * (P1.b() * P3.a() - P1.a() * P3.b()) + P3.d() * (P1.b() * P2.a() - P1.a() * P2.b());
		return Eigen::Vector3d(x / d, y / d, z / d);
	}
}

#endif
