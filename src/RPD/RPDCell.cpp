#include "RPD.h"

#include <vector>
#include <cmath>

// Eigen
#include <Eigen/Dense>

#define MAX 1e16

namespace RPD
{
    RPDCell::RPDCell()
    {
        cutted_planes_ = std::vector<RPDPlane>();
    }

    RPDCell::RPDCell(Eigen::Vector3d center)
    {
        cutted_planes_ = std::vector<RPDPlane>();
    }

    bool RPDCell::cut_by_plane(RPDPlane plane)
    {
        cutted_planes_.push_back(plane);

        return true;
    }

    double RPDCell::get_max_distance(Eigen::Vector3d p)
    {
        if (cutted_planes_.size() == 0)
        {
            return 0;
        }

        double max_f = -MAX;
        for (int i = 0; i < cutted_planes_.size(); i++)
        {
            double f = cutted_planes_[i].signed_distance_to_plane(p);
            max_f = std::max(max_f, f);
        }
        return max_f;
    }
} //namespace RPD
