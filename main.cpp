///////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2023, STEREOLABS.
//
// All rights reserved.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
///////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>

using namespace std;

#include "MonitorTracker.hh"
#include "FaceDetector.hh"
#include "MultiPoseTracker.hh"
#include "MultiPhantomDeformer.hh"

#define ENABLE_GL_GUI 0
#if ENABLE_GL_GUI
#include "GLViewer.hpp"
#else
#include <igl/opengl/glfw/Viewer.h>
#endif

// Handle the CTRL-C keyboard signal
#include <signal.h>
bool exit_app = false;
void nix_exit_handler(int s) {
    exit_app = true;
}
// Set the function to handle the CTRL-C
void SetCtrlHandler() {
    struct sigaction sigIntHandler;
    sigIntHandler.sa_handler = nix_exit_handler;
    sigemptyset(&sigIntHandler.sa_mask);
    sigIntHandler.sa_flags = 0;
    sigaction(SIGINT, &sigIntHandler, NULL);
}

void PrintLogo();
void VIR_camera_pose(string VIR_cams_json, map<string, Eigen::Affine3f>& camera_init_pose);

int main(int argc, char **argv)
{
    // 1030_take1
    // string sn_zed_head             = "24380057";
    // string sn_zed_body_monitor_bot = "20332044";
    // string sn_zed_body_monitor_top = "29134678";
    // string sn_zed_body_corner      = "21929705";
    // string sn_kinect_glass         = "000440922112";
    // string sn_kinect_monitor       = "000913194512";
    // string sn_webcam_face          = "0";
    // string VIR_Data_PATH = string(getenv("VIR_Data"));
    // string VIR_cams_json = VIR_Data_PATH + "record/1030/VIR_cams.json";
    // map<string, Eigen::Affine3f> camera_init_pose;
    // VIR_camera_pose(VIR_cams_json, camera_init_pose);
    // string monitor_camParam  = VIR_Data_PATH + "camera/kinect_QHD_intrinsic.yml";
    // string monitor_detParam  = VIR_Data_PATH + "camera/detector_params.yml";
    // string monitor_init_pose = VIR_Data_PATH + "record/1030/take1/monitor_a0.txt";
    // string record_monitor    = VIR_Data_PATH + "record/1030/take1/CAM_1.avi";
    // string record_face       = VIR_Data_PATH + "record/1030/take1/CAM_0.avi";
    // string record_zeds       = VIR_Data_PATH + "record/1030/take1/ZED_";
    // string train_model       = VIR_Data_PATH + "face_model/1030_Epoch_17.pt";
    // vector<string> member_id = { "0_SunghoMoon", "1_GaheeSon", "2_JaehyoKim", 
    //                              "3_BanghoShin", "4_HyeonilKim", "5_SuhyeonKim" };
    // vector<string> member_phantom = { "M_H170W65", "F_H160W55", "M_H170W65", 
    //                                   "M_H170W75", "M_H170W85", "F_H160W55" };

    // 1031_take18
    string sn_zed_head             = "24380057";
    string sn_zed_body_monitor_bot = "20332044";
    string sn_zed_body_monitor_top = "29134678";
    string sn_zed_body_corner      = "21929705";
    string sn_kinect_glass         = "000440922112";
    string sn_kinect_monitor       = "000913194512";
    string sn_webcam_face          = "0";
    string VIR_Data_PATH = string(getenv("VIR_Data"));
    string VIR_cams_json = VIR_Data_PATH + "record/1031/VIR_cams.json";
    map<string, Eigen::Affine3f> camera_init_pose;
    VIR_camera_pose(VIR_cams_json, camera_init_pose);
    string monitor_camParam  = VIR_Data_PATH + "camera/kinect_QHD_intrinsic.yml";
    string monitor_detParam  = VIR_Data_PATH + "camera/detector_params.yml";
    string monitor_init_pose = VIR_Data_PATH + "record/1031/take18_picc_88F/monitor_a0.txt";
    string record_monitor    = VIR_Data_PATH + "record/1031/take18_picc_88F/CAM_1.avi";
    string record_face       = VIR_Data_PATH + "record/1031/take18_picc_88F/CAM_0.avi";
    string record_zeds       = VIR_Data_PATH + "record/1031/take18_picc_88F/ZED_";
    string train_model       = VIR_Data_PATH + "face_model/1031_Epoch_17.pt";
    vector<string> member_id = { "0_Doctor", "1_Assist", "2_Radiographer_M0", 
                                 "3_Nurse_F0", "4_Nurse_F1", "5_Radiographer_M1", "6_Radiographer_M2" };
    vector<string> member_phantom = { "M_H170W65", "M_H170W65", "M_H170W65", 
                                      "F_H160W55", "F_H160W55", "M_H170W65", "M_H170W65" };
    
    PrintLogo();
    int cut = 200;
    if (argc > 1) cut = atoi(argv[1]);
    
    MonitorTracker* MT = new MonitorTracker(record_monitor);
    MT->Set_MonitorPose0(monitor_init_pose);
    MT->Set_CameraPose(camera_init_pose[sn_kinect_monitor]);
    MT->Set_Parameters(monitor_camParam, monitor_detParam);
    MT->CollectData(cut);
    map<int, Eigen::Affine3f> frameData_monitor_pose = MT->Get_FrameData_Monitor_Pose();
    map<int, cv::Mat> frameData_monitor_image = MT->Get_FrameData_Monitor_Image();

    FaceDetector* FD = new FaceDetector(record_face);
    // FD->Set_CameraPose(camera_init_pose[sn_webcam_face]);
    FD->Set_TrainModel(train_model);
    FD->Set_MemberIDs(member_id);
    FD->Initialize();
    // FD->CollectData(cut);
    // map<int, FACEINFO> frameData_faceInfo = FD->Get_FrameData_FaceInfo();
    
    MultiPoseTracker* MPT = new MultiPoseTracker(record_zeds);
    MPT->Set_TopView(sn_zed_head);
    MPT->Set_CornerView(sn_zed_body_corner);
    MPT->Set_Monitor_TopView(sn_zed_body_monitor_top);
    MPT->Set_Monitor_BotView(sn_zed_body_monitor_bot);
    MPT->Set_Cam_Json(VIR_cams_json);
    MPT->Set_MemberIDs(member_id);
    MPT->CollectData(cut, frameData_monitor_pose);
    map<int, map<int, cv::Mat>> frameData_camera_images = MPT->Get_FrameData_Images();


#if ENABLE_GL_GUI
    // GL VIEWER
    char key = ' ';
    bool gl_viewer_available = true;
    bool quit = false;
    int frame_count(0);
    GLViewer viewer;
    vector<sl::FusionConfiguration> configurations = MPT->Get_ZEDconf();
    vector<sl::CameraParameters> camera_parameters = MPT->Get_ZEDparam();
    viewer.init(argc, argv, camera_parameters, configurations);
    map<int, map<int, sl::Objects>> camera_objects = MPT->Get_FrameData_Objects();
    map<int, map<int, sl::Bodies>> camera_bodies = MPT->Get_FrameData_Bodies();
    map<int, map<int, sl::Transform>> camera_poses = MPT->Get_CamPoses();
    map<int, map<int, sl::Mat>> camera_pointClouds = MPT->Get_FrameData_PointClouds();
    for (int f=0; f<cut; f++)
    {
        cout << f << endl;
        viewer.updateData(camera_poses[f], camera_objects[f], camera_bodies[f], camera_pointClouds[f]);
        cv::waitKey(50);
        gl_viewer_available = viewer.isAvailable();
    }
    viewer.exit();
    return 0;
#else
    MultiPhantomDeformer* MPD = new MultiPhantomDeformer(member_phantom);

    igl::opengl::glfw::Viewer vr;
    vr.data().point_size = 5;
    vr.data().line_width = 2;
    vr.data().show_lines = false;
    vr.data().show_overlay_depth = false;
    vr.core().is_animating = false;
    vr.core().animation_max_fps = 5;
    vr.core().background_color = Eigen::Vector4f(0.f, 0.f, 0.f, 1.0);
    Eigen::RowVector3d white(1., 1., 1.);

    map<int, Eigen::MatrixXd> U;
    map<int, Eigen::MatrixXd> C_NEW;
    map<int, Eigen::MatrixXd> Vo; MPD->GetVo(Vo);
    map<int, Eigen::MatrixXi> Fo; MPD->GetFo(Fo);
    map<int, Eigen::MatrixXi> BE; MPD->GetBE(BE);

    map<int, Eigen::MatrixXd> U_dummy;
    map<int, Eigen::MatrixXd> C_NEW_dummy;
    map<int, Eigen::MatrixXd> Vo_dummy; MPD->GetDummyVo(Vo_dummy);
    map<int, Eigen::MatrixXi> Fo_dummy; MPD->GetDummyFo(Fo_dummy);
    map<int, Eigen::MatrixXi> BE_dummy; MPD->GetDummyBE(BE_dummy);

    Eigen::MatrixXd colours(7,3);
    colours <<  255,  0,  0, // red
                255, 165,  0, // organge
                0, 255,  0, // green
                0,  100, 255, // blue
                75,  50, 130, // indigo
                238, 130, 238, // violet
                255, 192, 203; // pink
                // 0, 255, 255; // cyan
    colours /= 255.0;
    for (int id=0; id<member_phantom.size(); id++) {
        vr.append_mesh();
        vr.data(id).set_mesh(Vo[id], Fo[id]);
        vr.data(id).set_colors(colours.row(id%6-1));
        vr.data(id).point_size = 5;
        vr.data(id).line_width = 2;
        vr.data(id).show_lines = false;
        vr.data(id).show_overlay_depth = false;
    }
    for (int id=member_phantom.size(); id<member_phantom.size()*2; id++) {
        vr.append_mesh();
        vr.data(id).set_mesh(Vo_dummy[id], Fo_dummy[id]);
        vr.data(id).set_colors(white);
        vr.data(id).point_size = 5;
        vr.data(id).line_width = 2;
        vr.data(id).show_lines = false;
        vr.data(id).show_overlay_depth = false;
    }

    Eigen::RowVector3d isocenter(151.3, 48.9, 90.1);
    Eigen::MatrixXd V_carm;
    Eigen::MatrixXi F_carm;
    igl::readPLY("c-arm.ply", V_carm, F_carm);
    Eigen::Affine3d a = Eigen::Translation3d(-isocenter) * Eigen::AngleAxisd(0.5*M_PI, -Eigen::Vector3d::UnitX()) * Eigen::Translation3d(isocenter);
    V_carm = (a.matrix().topLeftCorner(3,3) * V_carm.transpose()).transpose();
    V_carm.rowwise() += isocenter;
    vr.data().set_mesh(V_carm, F_carm);
    vr.data().set_colors(white);

    set<int> view_ids, view_ids_dummy;
    int frameNo(0);
    vr.callback_key_down = [&](igl::opengl::glfw::Viewer& _vr, unsigned char key, int mod)->bool {
        switch( key ) {
            case ' ':
                vr.core().is_animating = !vr.core().is_animating;
            case ']':
                FD->Run(true);
                MPT->Run(FD->Get_FaceInfo(), 
                         frameData_monitor_pose[frameNo] * camera_init_pose[sn_webcam_face]);
                // cv::Mat head_image = frameData_camera_images[frameNo][stoi(sn_zed_head)];
                // cv::imshow(sn_zed_head, head_image);
                // cv::Mat monitor_image = frameData_monitor_image[frameNo];
                // cv::resize(monitor_image, monitor_image, cv::Size(1280,720));
                // cv::imwshow(sn_kinect_monitor, monitor_image);
                // cv::waitKey(1);
                // monitor_image.release();
                MPD->RunPhantom(MPT->Get_PhantomPose());
                MPD->GetU(U);
                MPD->GetC_NEW(C_NEW);
                MPD->GetViewIds(view_ids);
                for (int id=0; id<member_phantom.size(); id++) {
                    if (view_ids.find(id) != view_ids.end())    vr.data(id).is_visible = true;
                    else                                        vr.data(id).is_visible = false;
                }
                MPD->RunDummy(MPT->Get_DummyPose());
                MPD->GetDummyU(U_dummy);
                MPD->GetDummyC_NEW(C_NEW_dummy);
                MPD->GetDummyViewIds(view_ids_dummy);
                for (int id=member_phantom.size(); id<member_phantom.size()*2; id++) {
                    if (view_ids_dummy.find(id) != view_ids_dummy.end())    vr.data(id).is_visible = true;
                    else                                                    vr.data(id).is_visible = false;
                }
                frameNo++;
                cv::waitKey(1);
                break;
        }
        vr.data().clear_edges();
        vr.data().clear_points();
        vr.data().add_edges(RowVector3d::Zero(), RowVector3d(100,0,0), RowVector3d(255,0,0));
        vr.data().add_edges(RowVector3d::Zero(), RowVector3d(0,100,0), RowVector3d(0,255,0));
        vr.data().add_edges(RowVector3d::Zero(), RowVector3d(0,0,100), RowVector3d(0,0,255));

        for (auto u: U) {
            vr.data(u.first).set_vertices(u.second);
            vr.data(u.first).compute_normals();
        }
        for (auto c: C_NEW) {
            vr.data(c.first).clear_points();
            vr.data(c.first).clear_edges();
            vr.data(c.first).add_points(c.second, white);
            for (int i=0; i<BE[c.first].rows(); i++) {
                vr.data(c.first).add_edges(c.second.row(BE[c.first](i,0)), c.second.row(BE[c.first](i,1)), colours.row(c.first%6-1));
            }
        }
        
        for (auto u: U_dummy) {
            vr.data(u.first).set_vertices(u.second);
            vr.data(u.first).compute_normals();
        }
        for (auto c: C_NEW_dummy) {
            vr.data(c.first).clear_points();
            vr.data(c.first).clear_edges();
            vr.data(c.first).add_points(c.second, white);
            for (int i=0; i<BE_dummy[c.first].rows(); i++) {
                vr.data(c.first).add_edges(c.second.row(BE_dummy[c.first](i,0)), c.second.row(BE_dummy[c.first](i,1)), white);
            }
        }
    };

    vr.callback_pre_draw = [&](igl::opengl::glfw::Viewer _vr)->bool {
        if (!vr.core().is_animating)
            return false;
    };
    vr.launch();
    delete MPD;
#endif

    delete MT;
    // delete FD;
    delete MPT;
    
    

    cout << "EXIT_SUCCESS" << endl;
    return EXIT_SUCCESS;
}


void PrintLogo()
{
    std::cout << "===================================" << std::endl;
    std::cout << " __  __     ______      ____       " << std::endl;
    std::cout << "/\\ \\/\\ \\   /\\__  _\\    /\\  _`\\     " << std::endl;
    std::cout << "\\ \\ \\ \\ \\  \\/_/\\ \\/    \\ \\ \\_\\ \\   " << std::endl;
    std::cout << " \\ \\ \\ \\ \\    \\ \\ \\     \\ \\ ,  /   " << std::endl;
    std::cout << "  \\ \\ \\_/ \\    \\_\\ \\__   \\ \\ \\\\ \\  " << std::endl;
    std::cout << "   \\ `\\___/    /\\_____\\   \\ \\_\\ \\_\\" << std::endl;
    std::cout << "    `\\/__/     \\/_____/    \\/_/\\/ /" << std::endl << std::endl;
    std::cout << "               Body Tracking Module" << std::endl;
    std::cout << "===================================" << std::endl;
}

void VIR_camera_pose(string VIR_cams_json, map<string, Eigen::Affine3f>& camera_init_pose) {
    ifstream ifs(VIR_cams_json);
    if (!ifs.is_open()) {
        cerr << VIR_cams_json << " is not opened" << endl;
        exit(1);
    }
    string dump, dump_prev;
    while(getline(ifs,dump)) {
        stringstream ss(dump);
        ss >> dump;
        if (dump == "\"input\":") {
            string sn = dump_prev.substr(1, dump_prev.size()-3);
            Eigen::Vector3f rvec, tvec;
            while (getline(ifs, dump)) {
                stringstream sa(dump);
                sa >> dump;
                if (dump == "\"rotation\":") {
                    getline(ifs, dump); float rx = stof(dump.substr(0,dump.size()-1));
                    getline(ifs, dump); float ry = stof(dump.substr(0,dump.size()-1));
                    getline(ifs, dump); float rz = stof(dump);
                    rvec = Eigen::Vector3f(rx,ry,rz);
                }
                if (dump == "\"translation\":") {
                    getline(ifs, dump); float tx = stof(dump.substr(0,dump.size()-1));
                    getline(ifs, dump); float ty = stof(dump.substr(0,dump.size()-1));
                    getline(ifs, dump); float tz = stof(dump);
                    tvec = Eigen::Vector3f(tx,ty,tz);
                    break;
                }
            }
            Eigen::Affine3f a = Eigen::Affine3f::Identity();
            a.translate(tvec);
            a.rotate( (Eigen::AngleAxisf(rvec.norm(), rvec.normalized())).toRotationMatrix() );;
            camera_init_pose[sn] = a;
        }
        dump_prev = dump;
    }
    ifs.close();
}