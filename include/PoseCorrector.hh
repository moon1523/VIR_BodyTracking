#ifndef POSECORRECTOR_HH
#define POSECORRECTOR_HH

#include <iostream>
#include <vector>

#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/SparseCore>
#include <Eigen/StdVector>

// Author: Sungho Moon

#define rad2deg 180.0 / M_PI
#define deg2rad M_PI / 180.0

using namespace std;

typedef vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>> 
    RotationList;

static Eigen::RowVector3d projectOntoPlane(const Eigen::RowVector3d& v, const Eigen::RowVector3d& n) {
    Eigen::RowVector3d normal = n.normalized();
    double distance = v.dot(normal);
    return v - distance * normal;
}

class PoseCorrector
{
public:
    PoseCorrector(const Eigen::MatrixXd& _C,
                  const Eigen::MatrixXi& _BE,
                  const Eigen::VectorXi& _P,
                  const Eigen::RowVector3d& _isocenter);
    ~PoseCorrector();

    bool ROOT_JointsCorrection();
    void CHEST_JointsCorrection();
    void LEGS_JointsCorrection();
    void ARMS_JointsCorrection();
    void LIMBS_JointsCorrection();
    void HEAD_xyCorrection();

    void use_default_Larm();
    void use_default_Rarm();

    void set_C1(const Eigen::MatrixXd& _C1) { C1 = _C1; }
    void set_C1_prev(const Eigen::MatrixXd& _C1_prev) { C1 = _C1_prev; }
    Eigen::MatrixXd get_C1() { return C1; }
    void set_head(const Eigen::RowVector3d& xyz) { C1.row(27) = xyz; }
    RotationList get_vQ() { return vQ; }
    RotationList get_vQ_corr() { return vQ_corr; }
    
private:
    void set_rotation(const Eigen::RowVector3d& _v0, const Eigen::RowVector3d& v1, 
                      const Eigen::Quaterniond& PQ, Eigen::Quaterniond &Q) {
        Eigen::Vector3d v0 = _v0.normalized();
        Eigen::Vector3d Rv0 = PQ.matrix().transpose() * v1.normalized().transpose();
        Eigen::Quaterniond R( Eigen::AngleAxisd( acos(v0.dot(Rv0)), v0.cross(Rv0).normalized() ) );
        R.normalize();
        if (std::isnan(R.w())) {
            R = Eigen::Quaterniond::Identity();
        }
        Q = PQ * R;
        Q.normalize();
    }

private:
    Eigen::MatrixXd C;
    Eigen::MatrixXi BE;
    Eigen::VectorXi P;
    Eigen::RowVector3d isocenter;
    Eigen::MatrixXd C1, C1_prev;
    Eigen::VectorXd bL;
    Eigen::MatrixXd nbV, bV;

    Eigen::RowVector3d root_trans, clav_trans;

    Eigen::Matrix3d root0, clav0, chest0, head0;
    Eigen::Matrix3d root1, clav1, chest1, head1;
    Eigen::Quaterniond rootQ, clavQ, chestQ, headQ;

    double angle;
    Eigen::RowVector3d proj_zPlane, proj_yPlane, proj_xPlane;
    bool isFront, isUp, isOut;
    
    RotationList vQ;
    RotationList vQ_corr;

    Eigen::Quaterniond q;
};

#endif
