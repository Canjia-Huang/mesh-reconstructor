#include "RPD.h"

// Eigen
#include <Eigen/Dense>

namespace RPD
{
    // RPD Plane
    RPDPlane::RPDPlane()
    {
        cor_ = Eigen::Vector3d(0, 0, 0);
        nor_ = Eigen::Vector3d(0, 0, 0);
        a_ = 0;
        b_ = 0;
        c_ = 0;
        d_ = 0;
    }

    RPDPlane::RPDPlane(Eigen::Vector3d cor,
        Eigen::Vector3d nor)
    {
        cor_ = cor;
        nor_ = nor.normalized();
        a_ = nor.x();
        b_ = nor.y();
        c_ = nor.z();
        d_ = -cor.dot(nor);
    }

    RPDPlane::RPDPlane(Eigen::Vector3d cor,
        Eigen::Vector3d nor,
        int op_id)
    {
        cor_ = cor;
        nor_ = nor.normalized();
        a_ = nor.x();
        b_ = nor.y();
        c_ = nor.z();
        d_ = -cor.dot(nor);
        opposite_idx_ = op_id;
    }

    bool RPDPlane::is_on_positive_side(Eigen::Vector3d p)
    {
        return (signed_distance_to_plane(p) > 0);
    }

    double RPDPlane::signed_distance_to_plane(Eigen::Vector3d p)
    {
        // return (p - cor_).dot(nor_);
        return a_ * p.x() + b_ * p.y() + c_ * p.z() + d_;
    }

    Eigen::Vector3d RPDPlane::cor()
    {
        return cor_;
    }

    Eigen::Vector3d RPDPlane::nor()
    {
        return nor_;
    }

    int RPDPlane::opposite_idx()
    {
        return opposite_idx_;
    }

    void RPDPlane::set_opposite_idx(int id)
    {
        opposite_idx_ = id;
    }

} //namespace RPD
