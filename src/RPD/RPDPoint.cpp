#include "RPD.h"

// Eigen
#include <Eigen/Dense>

namespace RPD
{
    RPDPoint::RPDPoint()
    {
        cor_ = Eigen::Vector3d(0, 0, 0);
        nor_ = Eigen::Vector3d(0, 0, 0);
    }

    RPDPoint::RPDPoint(double x,
        double y,
        double z)
    {
        cor_ = Eigen::Vector3d(x, y, z);
        nor_ = Eigen::Vector3d(0, 0, 0);
    }

    RPDPoint::RPDPoint(Eigen::Vector3d cor)
    {
        cor_ = cor;
        nor_ = Eigen::Vector3d(0, 0, 0);
    }

    RPDPoint::RPDPoint(Eigen::Vector3d cor,
        Eigen::Vector3d nor)
    {
        cor_ = cor;
        nor_ = nor;
    };

    RPDPoint::RPDPoint(Eigen::Vector3d cor,
        Eigen::Vector3d nor,
        double weight)
    {
        cor_ = cor;
        nor_ = nor;
        weight_ = weight;
    };
} //namespace RPD
