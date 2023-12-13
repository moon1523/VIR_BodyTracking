#ifndef MultiPoseTracker_HH_
#define MultiPoseTracker_HH_

#include "ClientPublisher.hh"

#include <sl/Camera.hpp>
#include <sl/Fusion.hpp>

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>

using namespace std;

typedef vector< tuple<string, float, Eigen::Vector3f> > FACEINFO;

class MultiPoseTracker
{
  public:
    MultiPoseTracker(string record_zeds);
    ~MultiPoseTracker();

    void CollectData(int cut, map<int, Eigen::Affine3f>& frameData_monitor_pose);
    void Run(const FACEINFO& _face_info, const Eigen::Affine3f& wcam_a);
    bool End() { return frame_count >= collected_frame ? true : false; }
    map<int, Eigen::MatrixXf> Get_PhantomPose() { return phantom_pose; }
    map<int, Eigen::MatrixXf> Get_DummyPose() { return dummy_pose; }
    bool IsFirstFaceDetected() { return isFirstFaceDetected; }

    // Access Functions
    void Set_TopView(string sn) { topView = stoi(sn); }
    void Set_CornerView(string sn) { cornerView = stoi(sn); }
    void Set_Monitor_TopView(string sn) { monitor_topView = stoi(sn); }
    void Set_Monitor_BotView(string sn) { monitor_botView = stoi(sn); }
    void Set_Cam_Json(string json) { json_file = json; }
    void Set_MemberIDs(vector<string> ids) { member_id = ids; } 
    vector<sl::FusionConfiguration> Get_ZEDconf() { return configurations; }
    vector<sl::CameraParameters> Get_ZEDparam() { return camera_parameters; }
    map<int, map<int, cv::Mat>> Get_FrameData_Images() { return camera_images; }
    map<int, map<int, sl::Objects>> Get_FrameData_Objects() { return camera_objects; }
    map<int, map<int, sl::Bodies>> Get_FrameData_Bodies() { return camera_bodies; }
    map<int, map<int, sl::Mat>> Get_FrameData_PointClouds() { return camera_pointClouds; }
    map<int, map<int, sl::Transform>> Get_CamPoses() { return camera_poses; }
    
  private:
    void convertAffineEigen2SL(const Eigen::Affine3f &inp, sl::Transform &out);
    void transformObjects(const Eigen::Affine3f &cam_pose, sl::Objects &camera_objects);
    void transformBodies(const Eigen::Affine3f &cam_pose, sl::Bodies &camera_bodies);
    float Point2LineDistance2(const Eigen::Vector3f& zed_head,
                              const Eigen::Vector3f& face_center,
                              const Eigen::Affine3f& wcam_a);

  private:
    string svo_prefix, json_file;
    std::vector<ClientPublisher> zeds;

    int topView;
    int cornerView;
    int monitor_topView;
    int monitor_botView;
    ofstream ofs_BV;
    int frame_count;
    map<int, map<int, cv::Mat>> camera_images;
    map<int, map<int, sl::Objects>> camera_objects;
    map<int, map<int, sl::Bodies>> camera_bodies;
    map<int, map<int, sl::Transform>> camera_poses;
    map<int, map<int, sl::Mat>> camera_pointClouds;
    int collected_frame;
    FACEINFO face_info;
    map<int, set<int>> fID2hIDs;
    set<set<int>> hIDs_sets;

    map<int, Eigen::MatrixXf> phantom_pose;
    map<int, Eigen::MatrixXf> dummy_pose;

    bool isFirstFaceDetected;
    vector<sl::FusionConfiguration> configurations;
    vector<sl::CameraParameters> camera_parameters;
    vector<string> member_id;
    
};


#endif