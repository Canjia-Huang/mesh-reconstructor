#ifndef RPDRECONSTRUCTION_H
#define RPDRECONSTRUCTION_H

#include <Export.h>
#include <vector>
#include <string>
#include <map>

// Eigen
#include <Eigen/Dense>

namespace RPD
{
	class MESH_RECONSTRUCTOR_EXPORT RPDPoint
	{
	public:
		RPDPoint();
		RPDPoint(double x,
			double y,
			double z);
		RPDPoint(Eigen::Vector3d cor);
		RPDPoint(Eigen::Vector3d cor,
			Eigen::Vector3d nor);
		RPDPoint(Eigen::Vector3d cor,
			Eigen::Vector3d nor,
			double weight);
		RPDPoint& operator =(const RPDPoint& p)
		{
			if (this != &p)
			{
				this->cor_ = p.cor_;
				this->nor_ = p.nor_;
				this->weight_ = p.weight_;
			}
			return *this;
		}
		bool operator <(const RPDPoint& p) const
		{
			if ((cor_.x() - p.cor_.x()) < 1e-12 && (cor_.x() - p.cor_.x()) > -1e-12)
			{
				if ((cor_.y() - p.cor_.y()) < 1e-12 && (cor_.y() - p.cor_.y()) > -1e-12)
				{
					return (cor_.z() < p.cor_.z());
				}
				return (cor_.y() < p.cor_.y());
			}
			return (cor_.x() < p.cor_.x());
		}
		bool operator ==(const RPDPoint& p) const
		{
			if ((cor_.x() - p.cor_.x()) < 1e-12 && (cor_.x() - p.cor_.x()) > -1e-12)
			{
				if ((cor_.y() - p.cor_.y()) < 1e-12 && (cor_.y() - p.cor_.y()) > -1e-12)
				{
					if ((cor_.z() - p.cor_.z()) < 1e-12 && (cor_.z() - p.cor_.z()) > -1e-12)
					{
						return true;
					}
				}
			}
			return false;
		}
	public:
		Eigen::Vector3d cor_;
		Eigen::Vector3d nor_;
		double			weight_ = 0;
		bool			is_feature_ = false;
	}; // class RPDPoint


	class MESH_RECONSTRUCTOR_EXPORT RPDPlane
	{
	public:
		RPDPlane();
		RPDPlane(Eigen::Vector3d cor,
			Eigen::Vector3d nor);
		RPDPlane(Eigen::Vector3d cor,
			Eigen::Vector3d nor,
			int op_id);

		bool is_on_positive_side(Eigen::Vector3d p);
		double signed_distance_to_plane(Eigen::Vector3d p);

		Eigen::Vector3d cor();
		Eigen::Vector3d nor();
		int opposite_idx();
		void set_opposite_idx(int id);
	private:
		Eigen::Vector3d cor_;
		Eigen::Vector3d nor_;
		double			a_;
		double			b_;
		double			c_;
		double			d_;
		int opposite_idx_ = -1;
	}; // class RPDPlane


	class MESH_RECONSTRUCTOR_EXPORT RPDDisk
	{
	private:
		int initialization_border(int seg,
			double radius);
	public:
		RPDDisk();
		RPDDisk(Eigen::Vector3d cor,
			Eigen::Vector3d nor);
		RPDDisk(Eigen::Vector3d cor,
			Eigen::Vector3d nor,
			double radius);
		RPDDisk(Eigen::Vector3d cor,
			Eigen::Vector3d nor,
			int seg);
		RPDDisk(Eigen::Vector3d cor,
			Eigen::Vector3d nor,
			int seg,
			double radius);

		bool cut_by_plane(RPDPlane plane);

		void output_disk_OBJ(std::string filepath);
	public:
		Eigen::Vector3d cor_;
		Eigen::Vector3d nor_;
		std::vector<RPDPoint> border_points_;
		std::vector<RPDPlane> border_planes_;
		int		segment_	= 8;
		double	radius_		= 0.1;
	}; // class RPDDisk


	class MESH_RECONSTRUCTOR_EXPORT RPDCell
	{
	private:
	public:
		RPDCell();
		RPDCell(Eigen::Vector3d center);
		bool cut_by_plane(RPDPlane plane);
		double get_max_distance(Eigen::Vector3d p);
	private:
		double				  size_ = 0.1;
	public:
		std::vector<RPDPlane> cutted_planes_;
	}; // class RPDCell


	class MESH_RECONSTRUCTOR_EXPORT RPDRecon
	{
	private:
		double safetyAcos(double value);
		bool is_feature_point(RPDPoint rp);
		bool is_feature_point(int i);
		void output_feature_points(std::string filename,
			double r, 
			double g, 
			double b);

		int PoissonSurfaceReconstruction(std::vector<Eigen::Vector3d>& points,
			std::vector<Eigen::Vector3i>& faces);
		int RegularTriangulation(std::vector<std::vector<int>>& RT_neighbors);

		void RecordTimeCount(std::string func_name,
			double time);
	public:
		RPDRecon();
		RPDRecon(std::vector<Eigen::Vector3d> cors,
			std::vector<Eigen::Vector3d> nors,
			std::vector<bool> is_features);
		int input_points(std::vector<Eigen::Vector3d> cors,
			std::vector<Eigen::Vector3d> nors,
			std::vector<bool> is_features);
		double radius();
		void set_radius(double r);
		void set_feature_weight(double fw);
		void set_not_feature_weight(double nfw);

		int PoissonBasedReconstruction(std::vector<Eigen::Vector3i>& output_faces);
		int ParallelDiskReconstruction(std::vector<Eigen::Vector3i>& output_faces);

		void PrintTimeCount();

	private:
		std::vector<RPDPoint> points_;
		std::vector<bool>	  is_features_;

		double				  radius_			  = 0;

		double				  feature_weight_	  = 0;
		double				  not_feature_weight_ = 0;

		double				  connect_angle_constraint_		  = 100;
		double				  connect_general_distance_ratio_ = 3;
		double				  connect_feature_distance_ratio_ = 5;

		std::map<std::string, double> time_count_;
	public:
	}; // class RPDRecon

}; // namespace RPD

#endif
