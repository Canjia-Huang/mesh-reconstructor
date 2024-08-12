#include "RPD.h"

// Eigen
#include <Eigen/Dense>

namespace RPD
{
    // RPD Plane
    RPDPlane::RPDPlane()
    {
        cor_ = Eigen::Vector3d(0, 0, 0);
        nor_ = Eigen::Vector3d(1, 0, 0);
    }

    RPDPlane::RPDPlane(Eigen::Vector3d cor,
        Eigen::Vector3d nor)
    {
        cor_ = cor;
        nor_ = nor.normalized();
    }

    RPDPlane::RPDPlane(Eigen::Vector3d cor,
        Eigen::Vector3d nor,
        int op_id)
    {
        cor_ = cor;
        nor_ = nor.normalized();
        opposite_idx_ = op_id;
    }

    bool RPDPlane::is_on_positive_side(Eigen::Vector3d p)
    {
        return ((p - cor_).dot(nor_) > 0);
    }

    double RPDPlane::signed_distance_to_plane(Eigen::Vector3d p)
    {
        return (p - cor_).dot(nor_);
    }
} //namespace RPD
