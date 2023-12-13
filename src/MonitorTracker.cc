#include "MonitorTracker.hh"


MonitorTracker::MonitorTracker(string record_file)
: nMarkers(4), cumulCount(0)
{
    params = cv::aruco::DetectorParameters::create();
    dictionary = cv::aruco::generateCustomDictionary(nMarkers,4,2);
    markerLength = 18;
    
    if (!record_file.empty()) {
        cap.open(record_file);
    }
    else {
        cap.open(findCameraIndex(), cv::CAP_V4L2);
        cap.set(cv::CAP_PROP_FOURCC, cv::VideoWriter::fourcc('M','J','P','G'));
        cap.set(cv::CAP_PROP_FRAME_WIDTH,  2560);
        cap.set(cv::CAP_PROP_FRAME_HEIGHT, 1440);
        cap.set(cv::CAP_PROP_FPS, 15);
        cap.set(cv::CAP_PROP_AUTOFOCUS, 0);
    }
    
    if (!cap.isOpened()) {
        cerr << "Error! Unabale to open camera" << endl;
        exit(1);
    }
}

void MonitorTracker::CollectData(int frameNo)
{
    for (int f=0; f<frameNo; f++)
    {
        cap >> image;
        Run(image);
        frameData_monitor_pose[f] = Get_Monitor_Pose();
        frameData_monitor_image[f] = image;
        Render();
        char key = (char)cv::waitKey(1);
        if (key == 'q')
            break;
        cout << "\rcollected " << f << " frames" << flush;
    }
    cap.release();
    cv::destroyAllWindows();
}

void MonitorTracker::Set_Parameters(string camParam, string detParam)
{
    cv::FileStorage fs(camParam, cv::FileStorage::READ);
    if (!fs.isOpened())
        return;
    fs["camera_matrix"] >> camMatrix;
    fs["distortion_coefficients"] >> distCoeffs;
    cv::FileStorage fs2(detParam, cv::FileStorage::READ);
    if (!fs2.isOpened())
        return;
    fs2["adaptiveThreshWinSizeMin"] >> params->adaptiveThreshWinSizeMin;
    fs2["adaptiveThreshWinSizeMax"] >> params->adaptiveThreshWinSizeMax;
    fs2["adaptiveThreshWinSizeStep"] >> params->adaptiveThreshWinSizeStep;
    fs2["adaptiveThreshConstant"] >> params->adaptiveThreshConstant;
    fs2["minMarkerPerimeterRate"] >> params->minMarkerPerimeterRate;
    fs2["maxMarkerPerimeterRate"] >> params->maxMarkerPerimeterRate;
    fs2["polygonalApproxAccuracyRate"] >> params->polygonalApproxAccuracyRate;
    fs2["minCornerDistanceRate"] >> params->minCornerDistanceRate;
    fs2["minDistanceToBorder"] >> params->minDistanceToBorder;
    fs2["minMarkerDistanceRate"] >> params->minMarkerDistanceRate;
    fs2["cornerRefinementMethod"] >> params->cornerRefinementMethod;
    fs2["cornerRefinementWinSize"] >> params->cornerRefinementWinSize;
    fs2["cornerRefinementMaxIterations"] >> params->cornerRefinementMaxIterations;
    fs2["cornerRefinementMinAccuracy"] >> params->cornerRefinementMinAccuracy;
    fs2["markerBorderBits"] >> params->markerBorderBits;
    fs2["perspectiveRemovePixelPerCell"] >> params->perspectiveRemovePixelPerCell;
    fs2["perspectiveRemoveIgnoredMarginPerCell"] >> params->perspectiveRemoveIgnoredMarginPerCell;
    fs2["maxErroneousBitsInBorderRate"] >> params->maxErroneousBitsInBorderRate;
    fs2["minOtsuStdDev"] >> params->minOtsuStdDev;
    fs2["errorCorrectionRate"] >> params->errorCorrectionRate;
}

MonitorTracker::~MonitorTracker()
{
    cv::destroyAllWindows();
}

void MonitorTracker::Run(cv::Mat &image)
{
    image.copyTo(display);

    // detect markers
    vector<int> ids;
    vector<vector<cv::Point2f>> corners, rejected, cornersWhole;
    cv::aruco::detectMarkers(image, dictionary, corners, ids, params, rejected);
    
    if (ids.size() == 0)
        return;
    
    cv::aruco::drawDetectedMarkers(display, corners, ids);

    // check the status
    bool isNewPose(false);
    for (int i=0; i<ids.size(); i++)
    {
        vector<cv::Point2f> points;
        for (auto p: corners[i]) {
            points.push_back( cv::Point2f(p.x, p.y) );
        }
        cornersWhole.push_back(points);
        if (isNewPose)
            continue;
        if (corner_cumul.find(ids[i]) != corner_cumul.end())
        {
            cv::Point2f oldCen(0,0), newCen(0,0);
            for (int n=0; n<4; n++)
            {
                newCen += points[n];
                oldCen += corner_cumul[ids[i]][n];
            }
            cv::Point2f vec = oldCen - newCen;
            if (vec.dot(vec) > 10) // distance check
            {
                corner_cumul.clear();
                isNewPose = true;
                As.clear();
                rveci_cumuls.clear();
                tveci_cumuls.clear();
                ai_cumuls.clear();
                cumulCount = 0;
            }
        }
    }
    for (int i=0; i<ids.size(); i++)
        corner_cumul[ids[i]] = cornersWhole[i];

    vector<cv::Vec3d> rvecs, tvecs; // must double type
    cv::aruco::estimatePoseSingleMarkers(corners, markerLength, camMatrix, distCoeffs, rvecs, tvecs);

    vector<tuple<int, cv::Vec3f, cv::Vec3f>> poses;
    for (int i=0; i<ids.size(); i++)
        poses.push_back(make_tuple(ids[i], rvecs[i], tvecs[i]));
    // sort by id
    sort(poses.begin(), poses.end(), [](const tuple<int, cv::Vec3f, cv::Vec3f> &a, const tuple<int, cv::Vec3f, cv::Vec3f> &b) {
        return get<0>(a) < get<0>(b);
    });
    cumulCount++;

    map<int, cv::Vec3f> _tvecs;
    for (auto itr: poses) {
        int i = get<0>(itr);
        cv::Vec3f rr = get<1>(itr);
        cv::Vec3f tt = get<2>(itr);

        // rvec
        float angle = norm(rr);
        Eigen::Vector3f axis( rr(0) / angle, rr(1) / angle, rr(2) / angle);
        Eigen::Quaternionf q(Eigen::AngleAxisf(angle, axis));
        q.normalize();

        if (cumulCount > 1) {
            qi_cumuls[i] = qi_cumuls[i].slerp(1.f/ (cumulCount + 1.f), q);                
        }
        else {
            qi_cumuls[i] = q;
        }

        Eigen::AngleAxisf avg(qi_cumuls[i]);
        cv::Vec3f rvec;
        cv::eigen2cv(avg.axis(), rvec);
        rvec *= avg.angle();
        rveci_cumuls[i] = rvec;

        // tvec
        _tvecs[i] += tt;
        if (cumulCount > 1) {
            _tvecs[i] += tveci_cumuls[i] * (cumulCount - 1);
            _tvecs[i] /= (float)cumulCount;
        }
        tveci_cumuls[i] = _tvecs[i];
        
        cv::cv2eigen(tveci_cumuls[i], ti_cumuls[i]);
        
        ai_cumuls[i].linear() = qi_cumuls[i].normalized().toRotationMatrix();
        ai_cumuls[i].translation() = ti_cumuls[i] * 0.01;


        if (ai_cumuls.find(i) != ai_cumuls.end()) {
            Eigen::Affine3f A = A_viewCam * ai_cumuls[i] * a0_cumuls[i].inverse() * A_viewCam_inverse;
            As.push_back( A );
            cv::aruco::drawAxis(display, camMatrix, distCoeffs, rveci_cumuls[i], tveci_cumuls[i], markerLength * 0.5f);
        }
    }
    if (As.size() > 0) {
        A_monitor = averaging_Affine3f(As);
    }
        

    cv::putText(display, "Tracking monitor ...", cv::Point(10, 20), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(255.f, 0.f, 0.f), 1.2);
    cv::putText(display, "cummulated data #: " + to_string(cumulCount), cv::Point(10, 40), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(255.f, 0.f, 0.f), 1.2);
    

}

void MonitorTracker::Render()
{
    cv::resize(display, display, cv::Size(1280, 720));
    cv::imshow("Render", display);
}

void MonitorTracker::Render(cv::Mat& _image)
{
    cv::resize(_image, _image, cv::Size(1280, 720));
    cv::imshow("Render", _image);
}