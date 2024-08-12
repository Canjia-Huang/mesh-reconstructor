#include "RPD.h"

#include <iostream>
#include <fstream>

// Eigen
#include <Eigen/Dense>

#define PI 3.1415926
#define EPS 1e-16
#define SIGN(x) ((x<0)? -1 : 1)

namespace RPD
{
    RPDDisk::RPDDisk()
    {
        cor_ = Eigen::Vector3d(0, 0, 0);
        nor_ = Eigen::Vector3d(0, 0, 0);
    }

    RPDDisk::RPDDisk(Eigen::Vector3d cor,
        Eigen::Vector3d nor)
    {
        cor_ = cor;
        nor_ = nor;

        initialization_border(segment_, radius_);
    }

    RPDDisk::RPDDisk(Eigen::Vector3d cor,
        Eigen::Vector3d nor,
        double radius)
    {
        cor_ = cor;
        nor_ = nor;
        radius_ = radius;

        initialization_border(segment_, radius_);
    }

    RPDDisk::RPDDisk(Eigen::Vector3d cor,
        Eigen::Vector3d nor,
        int seg)
    {
        cor_ = cor;
        nor_ = nor;
        segment_ = seg;

        initialization_border(segment_, radius_);
    }

    RPDDisk::RPDDisk(Eigen::Vector3d cor,
        Eigen::Vector3d nor,
        int seg,
        double radius)
    {
        cor_ = cor;
        nor_ = nor;
        segment_ = seg;
        radius_ = radius;

        initialization_border(segment_, radius_);
    }

    int RPDDisk::initialization_border(int seg,
        double radius)
    {
        if (seg < 3)
        {
            std::cout << "RPDDisk::initialization_border input seg is invalid!" << std::endl;
            return 0;
        }

        Eigen::Vector3d init_vec(0, 0, 0);
        double dot_value = cor_.dot(nor_);
        if (cor_.norm() < EPS)
        {
            if (fabs(nor_.x()) < EPS)
            {
                init_vec.x() = 1;
            }
            else
            {
                if (fabs(nor_.y()) < EPS)
                {
                    init_vec.y() = 1;
                }
                else
                {
                    if (fabs(nor_.z()) < EPS)
                    {
                        init_vec.z() = 1;
                    }
                    else
                    {
                        init_vec.x() = 1;
                        init_vec.z() = -nor_.x() / nor_.z();
                    }
                }
            }
        }
        else if (fabs(dot_value) < EPS)
        {
            init_vec = cor_;
        }
        else
        {
            if (fabs(nor_.x()) < EPS)
            {
                if (fabs(nor_.y()) < EPS)
                {
                    if (fabs(nor_.z()) < EPS)
                    {
                        std::cout << "RPDDisk::initialization_border normal is invalid!" << std::endl;
                        return 0;
                    }
                    else
                    {
                        init_vec.z() = dot_value / nor_.z();
                    }
                }
                else
                {
                    init_vec.y() = dot_value / nor_.y();
                }
            }
            else
            {
                init_vec.x() = dot_value / nor_.x();
            }

            init_vec = init_vec - cor_;
        }
        // init_vec.normalize();

        for (int i = 0; i < seg; i++)
        {
            double ang = 2 * PI / double(seg) * double(i);
            Eigen::Vector3d res = (cos(ang) * init_vec + (1 - cos(ang)) * (init_vec.dot(nor_)) * nor_ + sin(ang) * (nor_.cross(init_vec))).normalized();

            border_points_.push_back(RPDPoint(cor_ + radius * res));
        }
        for (int i = 0; i < seg; i++)
        {
            int next_i = i + 1;
            if (next_i == seg)
            {
                next_i = 0;
            }

            Eigen::Vector3d mid_point = 0.5 * (border_points_[i].cor_ + border_points_[next_i].cor_);
            border_planes_.push_back(RPDPlane(mid_point, cor_ - mid_point, -1));
        }

        return 1;
    }

    bool RPDDisk::cut_by_plane(RPDPlane plane)
    {
        bool is_cutted = false;

        int center_sign = SIGN(plane.signed_distance_to_plane(cor_));

        std::vector<double> distance_to_plane;
        std::vector<int>    sign_to_plane;
        for (int i = 0; i < border_points_.size(); i++)
        {
            double dis = plane.signed_distance_to_plane(border_points_[i].cor_);
            distance_to_plane.push_back(fabs(dis));
            sign_to_plane.push_back(SIGN(dis));

            if (SIGN(dis) != center_sign)
            {
                is_cutted = true;
            }
        }

        // process
        if (is_cutted)
        {
            int e1p1, e1p2, e2p1, e2p2;
            int cutted_edge_cnt = 0;
            for (int i = 0; i < border_points_.size(); i++)
            {
                int next_i = i + 1;
                if (next_i == border_points_.size())
                {
                    next_i = 0;
                }

                int sign1 = sign_to_plane[i];
                int sign2 = sign_to_plane[next_i];

                if (sign1 != sign2)
                {
                    cutted_edge_cnt++;
                    if (sign1 == center_sign)
                    {
                        e1p1 = i;
                        e1p2 = next_i;
                    }
                    else if (sign2 == center_sign)
                    {
                        e2p1 = i;
                        e2p2 = next_i;
                    }
                }
            }
            if (cutted_edge_cnt != 2)
            {
                std::cout << "RPDDisk::cut_by_plane error!" << std::endl;
                return false;
            }


            RPDPoint new_point1(border_points_[e1p1].cor_ + (border_points_[e1p2].cor_ - border_points_[e1p1].cor_) * (distance_to_plane[e1p1] / (distance_to_plane[e1p1] + distance_to_plane[e1p2])));
            RPDPoint new_point2(border_points_[e2p1].cor_ + (border_points_[e2p2].cor_ - border_points_[e2p1].cor_) * (distance_to_plane[e2p1] / (distance_to_plane[e2p1] + distance_to_plane[e2p2])));

            // renew
            std::vector<RPDPoint> border_points_tmp;
            std::vector<RPDPlane> border_planes_tmp;
            for (int i = 0; i < border_points_.size(); i++)
            {
                if (sign_to_plane[i] != center_sign)
                {
                    continue;
                }

                if (i == e1p1)
                {
                    border_points_tmp.push_back(border_points_[e1p1]);
                    border_planes_tmp.push_back(border_planes_[e1p1]);

                    border_points_tmp.push_back(new_point1);
                    border_planes_tmp.push_back(plane);

                    border_points_tmp.push_back(new_point2);
                    border_planes_tmp.push_back(border_planes_[e2p1]);
                }
                else
                {
                    border_points_tmp.push_back(border_points_[i]);
                    border_planes_tmp.push_back(border_planes_[i]);
                }
            }
            border_points_ = border_points_tmp;
            border_planes_ = border_planes_tmp;
        }

        return is_cutted;
    }

    void RPDDisk::output_disk_OBJ(std::string filepath)
    {
        std::ofstream out(filepath);
        out << "v " << cor_.transpose() << std::endl;
        for (int i = 0; i < border_points_.size(); i++)
        {
            out << "v " << border_points_[i].cor_.transpose() << std::endl;

            if (i == border_points_.size() - 1)
            {
                out << "f " << 1 << " " << i + 2 << " " << 2 << std::endl;
            }
            else
            {
                out << "f " << 1 << " " << i + 2 << " " << i + 3 << std::endl;
            }
        }
        out.close();
    }
} //namespace RPD
