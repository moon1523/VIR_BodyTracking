#ifndef MultiPhantomDeformer_HH_
#define MultiPhantomDeformer_HH_

#include <iostream>
#include <fstream>
#include <sstream>
#include <functional>
#include <algorithm>
#include <vector>

#include <igl/readPLY.h>
#include <igl/readTGF.h>
#include <igl/directed_edge_parents.h>
#include <igl/LinSpaced.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include <igl/writePLY.h>
#include <igl/bbw.h>
#include <igl/boundary_conditions.h>
#include <igl/forward_kinematics.h>
#include <igl/dqs.h>

#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/SparseCore>
#include <Eigen/StdVector>

#include "functions.hh"
#include "PoseCorrector.hh"


using namespace std;

class MultiPhantomDeformer
{
  public:
    MultiPhantomDeformer(const vector<string>& member_phantom);
    ~MultiPhantomDeformer();

    void RunPhantom(const map<int, Eigen::MatrixXf>& phantom_pose);
    void RunDummy(const map<int, Eigen::MatrixXf>& dummy_pose);
    void LoadPhantom();
    void LoadDummy();

    void GetU(map<int, Eigen::MatrixXd>& U) { U = this->U; }
    void GetVo(map<int, Eigen::MatrixXd>& Vo) { Vo = this->Vo; }
    void GetFo(map<int, Eigen::MatrixXi>& Fo) { Fo = this->Fo; }
    void GetC1(map<int, Eigen::MatrixXd>& C1) { C1 = this->C1; }
    void GetBE(map<int, Eigen::MatrixXi>& BE) { BE = this->BE; }
    void GetC_NEW(map<int, Eigen::MatrixXd>& C_NEW) { C_NEW = this->C_NEW; }
    void GetViewIds(set<int>& view_ids) { view_ids = this->view_ids; }

    void GetDummyU(map<int, Eigen::MatrixXd>& U) { U = this->U_dummy; }
    void GetDummyVo(map<int, Eigen::MatrixXd>& Vo) { Vo = this->Vo_dummy; }
    void GetDummyFo(map<int, Eigen::MatrixXi>& Fo) { Fo = this->Fo_dummy; }
    void GetDummyC1(map<int, Eigen::MatrixXd>& C1) { C1 = this->C1_dummy; }
    void GetDummyBE(map<int, Eigen::MatrixXi>& BE) { BE = this->BE_dummy; }
    void GetDummyC_NEW(map<int, Eigen::MatrixXd>& C_NEW) { C_NEW = this->C_NEW_dummy; }
    void GetDummyViewIds(set<int>& view_ids) { view_ids = this->view_ids_dummy; }
    
    
  private:
    vector<int> tgf2zed;
    vector<string> member_phantom;

    map<int, Eigen::MatrixXd> C;
    map<int, Eigen::MatrixXi> BE;
    map<int, Eigen::VectorXi> P;
    map<int, Eigen::MatrixXd> Vo;
    map<int, Eigen::MatrixXi> Fo;
    map<int, Eigen::MatrixXd> Wo;
    map<int, Eigen::MatrixXi> FT;
    map<int, Eigen::MatrixXd> VT;
    map<int, Eigen::MatrixXi> TT;
    map<int, vector<map<int, double>>> weightMap;
    map<int, PoseCorrector*> pc;
    map<int, Eigen::MatrixXd> C1, C1_prev;
    map<int, Eigen::MatrixXd> U;
    map<int, Eigen::MatrixXd> C_NEW;
    set<int> view_ids;

    map<int, Eigen::MatrixXd> C_dummy;
    map<int, Eigen::MatrixXi> BE_dummy;
    map<int, Eigen::VectorXi> P_dummy;
    map<int, Eigen::MatrixXd> Vo_dummy;
    map<int, Eigen::MatrixXi> Fo_dummy;
    map<int, Eigen::MatrixXd> Wo_dummy;
    map<int, Eigen::MatrixXi> FT_dummy;
    map<int, Eigen::MatrixXd> VT_dummy;
    map<int, Eigen::MatrixXi> TT_dummy;
    map<int, vector<map<int, double>>> weightMap_dummy;
    map<int, PoseCorrector*> pc_dummy;
    map<int, Eigen::MatrixXd> C1_dummy, C1_prev_dummy;
    map<int, Eigen::MatrixXd> U_dummy;
    map<int, Eigen::MatrixXd> C_NEW_dummy;
    set<int> view_ids_dummy;
};


#endif