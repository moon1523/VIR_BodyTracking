#ifndef FaceDetector_HH_
#define FaceDetector_HH_

#include <Python.h>
#include <numpy/arrayobject.h>
#include <iostream>
#include <vector>
#include <tuple>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>

using namespace std;

typedef vector< tuple<string, float, Eigen::Vector3f> > FACEINFO;

class FaceDetector
{
  public:
    FaceDetector(string record_file);
    ~FaceDetector();

    void Initialize();
    void Run(bool isView=false);
    void CollectData(int frameNo);
    
    void Set_CameraPose(const Eigen::Affine3f& a) { A_wcam = a; }
    void Set_TrainModel(string _train_model)  { train_model = _train_model; }
    void Set_MemberIDs(vector<string> _ids) {
      member_id = _ids;
      num_workers = to_string(member_id.size());
    }
    FACEINFO Get_FaceInfo()      { return face_info; }
    string   Get_Name(int id)      { return get<0>(face_info[id]); }
    float    Get_Similarity(int id) { return get<1>(face_info[id]); }
    Eigen::Vector3f Get_3D_center(int id) { return get<2>(face_info[id]); }
    map<int, FACEINFO> Get_FrameData_FaceInfo() { return frameData_faceInfo; }
  
  private:
    void Import_Modules();
    void Load_Models();
    
  private:
    string train_model;
    string num_workers;
    vector<string> member_id;

    PyObject* main_module;
    PyObject* global_dict;
    
    FACEINFO face_info;
    float focal_length;
    int frameNo;

    Eigen::Affine3f A_wcam;
    string record_file;
    map<int, FACEINFO> frameData_faceInfo;

};


#endif