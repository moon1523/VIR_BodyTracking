#include "MultiPoseTracker.hh"

MultiPoseTracker::MultiPoseTracker(string record_zeds)
: frame_count(-1), svo_prefix(record_zeds), isFirstFaceDetected(false)
{
    cout << "MultiPoseTracker ==========================" << endl;
    ofs_BV.open("BV.txt");
}

MultiPoseTracker::~MultiPoseTracker()
{
    ofs_BV.close();
}

void MultiPoseTracker::CollectData(int cut, map<int, Eigen::Affine3f>& frameData_monitor_pose)
{
    cout << "svo_prefix: " << svo_prefix << endl;
    cout << "json_file: " << json_file << endl;
    cout << "topView: " << topView << endl;
    cout << "cut: " << cut << endl;

    // Read json file containing the configuration of your multicamera setup.
    configurations  = sl::readFusionConfigurationFile(json_file, sl::COORDINATE_SYSTEM::IMAGE, sl::UNIT::METER);
    if (configurations.empty()) {
        std::cout << "Empty configuration File." << std::endl;
        exit(1);
    }

    // Check if the ZED camera should run within the same process or if they are running on the edge.
    Eigen::Affine3f camPose_Eigen; 
    sl::Transform camPose_SL;
    for (sl::FusionConfiguration &conf : configurations)
    {
        ClientPublisher zed;
        if (conf.serial_number != topView && conf.serial_number != cornerView &&
            conf.serial_number != monitor_topView && conf.serial_number != monitor_botView)
            continue;
        // if the ZED camera should run locally, then start a thread to handle it
        std::cout << "Try to open ZED " << conf.serial_number << ".." << std::flush;
        auto state = zed.open(conf.input_type, conf.serial_number, svo_prefix);
        if (state)
        {
            std::cout << ". ready !" << std::endl;
            sl::float3 t_sl = conf.pose.getTranslation();
            sl::float4 q_sl = conf.pose.getOrientation();
            Eigen::Translation3f t(t_sl.x, t_sl.y, t_sl.z);
            Eigen::Quaternionf q(q_sl.w, q_sl.x, q_sl.y, q_sl.z); q.normalize();
            Eigen::Affine3f a = t * q;

            camera_parameters.push_back(zed.getCameraParameters());
            zed.setCameraPose_Eigen(a);
            zed.setCameraPose_SL(conf.pose);

            // main loop
            
            for (int f=0; f<cut; f++)
            {
                if (!zed.next())
                    break;
                
                if (conf.serial_number == monitor_topView || conf.serial_number == monitor_botView)
                    camPose_Eigen = frameData_monitor_pose[f] * zed.getCameraPose_Eigen();
                else
                    camPose_Eigen = zed.getCameraPose_Eigen();
                
                convertAffineEigen2SL(camPose_Eigen, camPose_SL);
                camera_poses[f][conf.serial_number] = camPose_SL;

                if (conf.serial_number == topView)
                {
                    auto obj = zed.getObjects();
                    if (!obj.is_tracked)
                        continue;
                    transformObjects(camPose_Eigen, obj);
                    obj.object_list.erase(std::remove_if(obj.object_list.begin(), obj.object_list.end(),[](sl::ObjectData iter)
                                        { return iter.head_position.z < 1.3; }),
                                        obj.object_list.end());
                    obj.object_list.erase(std::remove_if(obj.object_list.begin(), obj.object_list.end(),[](sl::ObjectData iter)
                                        { return iter.head_position.z > 2.f; }),
                                        obj.object_list.end());
                    camera_objects[f][conf.serial_number] = obj;
                }
                
                auto body = zed.getBodies();
                if (!body.is_tracked)
                    continue;
                transformBodies(camPose_Eigen, body);
                camera_bodies[f][conf.serial_number] = body;

                camera_images[f][conf.serial_number] = zed.getImage_CV();
                
                cout << "\rcollected " << f << " frames" << flush;
            }
            zed.close();
            cout << endl;
        }
    }
    collected_frame = camera_objects.size();
}

void MultiPoseTracker::Run(const FACEINFO& _face_info, const Eigen::Affine3f& wcam_a)
{
    frame_count++;
    cout << "Frame count: " << frame_count << endl;
    phantom_pose.clear();
    if (frame_count >= collected_frame)
        return;

    sl::Bodies fusedBody, fusedBody_BV;
    fusedBody.body_format = sl::BODY_FORMAT::BODY_34;
    fusedBody_BV.body_format = sl::BODY_FORMAT::BODY_34;
    ofs_BV << "f " << frame_count << endl;
    if (!camera_objects[frame_count][topView].is_tracked)
        return;

    // object tracking - head
    int det_n = camera_objects[frame_count][topView].object_list.size(); // detected head # in topView
    for (auto &obj : camera_objects[frame_count][topView].object_list) {
        if (obj.tracking_state != sl::OBJECT_TRACKING_STATE::OK)
            det_n--;
    }
    if (det_n == 0) return;
    Eigen::MatrixXf ref_Heads(det_n, 2);
    Eigen::MatrixXf ref_Heads_3D(det_n, 3);
    Eigen::VectorXi ref_IDs(det_n);
    int cnt_n = 0;
    for (auto &obj : camera_objects[frame_count][topView].object_list) {
        if (obj.tracking_state != sl::OBJECT_TRACKING_STATE::OK)
            continue;
        ref_Heads.row(cnt_n) = Eigen::RowVector2f(obj.head_position.x, obj.head_position.y);
        ref_Heads_3D.row(cnt_n) = Eigen::RowVector3f(obj.head_position.x, obj.head_position.y, obj.head_position.z);
        ref_IDs(cnt_n++) = obj.id;
    }
    
    // body tracking - head
    map<int, Eigen::MatrixXf> camera_heads_oneframe;
    for (auto it : camera_bodies[frame_count]) {
        // it.first: sn, it.second: bodies
        if (it.second.body_list.size() == 0)
            continue;
        Eigen::MatrixXf dat(it.second.body_list.size(), 2);
        for (int i = 0; i < dat.rows(); i++) {
            auto head = it.second.body_list[i].keypoint[26];
            dat.row(i) = Eigen::RowVector2f(head.x, head.y); // dat.row(body id) = head position
        }
        camera_heads_oneframe[it.first] = dat; // heads[sn] = dat
    }

    // matching object/body tracking head IDs
    for (auto &obj: camera_objects[frame_count][topView].object_list) {
        if (obj.tracking_state != sl::OBJECT_TRACKING_STATE::OK)
            continue;
        Eigen::RowVector2f center(obj.head_position.x, obj.head_position.y);
        for (auto dat: camera_heads_oneframe) {
            Eigen::Index idx, idx0;
            // if the head is too far from the center, skip it
            if ((dat.second.rowwise() - center).rowwise().squaredNorm().minCoeff(&idx) > 0.25)
                continue;
            (ref_Heads.rowwise() - dat.second.row(idx)).rowwise().squaredNorm().minCoeff(&idx0);
            // matching IDs
            if (ref_IDs(idx0) != obj.id)
                continue;
        }
    }

    cout << "head ids: " << ref_IDs.transpose() << endl;
    // matching IDs
    face_info = _face_info;
    for (int i=0; i<face_info.size(); i++) {
        if (!isFirstFaceDetected) isFirstFaceDetected = true;
        string fid_name = get<0>(face_info[i]);
        int fid = stoi(fid_name.substr(0, fid_name.find("_")));
        string name = fid_name.substr(fid_name.find("_")+1);
        double sim = get<1>(face_info[i]);
        // Eigen::Vector3f face_center = get<2>(face_info[i]);
        Eigen::Vector3f face_center = get<2>(face_info[i]);

        float min_dist2 = FLT_MAX;
        int idx(-1);
        for (int j=0; j<ref_Heads_3D.rows(); j++) {
            Eigen::Vector3f zed_head = ref_Heads_3D.row(j);
            float dist2 = Point2LineDistance2(zed_head, face_center, wcam_a);
            if (dist2 > 4000) continue; // 20 cm
            if (min_dist2 > dist2) {
                min_dist2 = dist2;
                idx = ref_IDs(j);
            }
        }
        cout << "head/face/name/sim : " << idx << "/" << fid << "/" << name << "/" << sim << endl;
        if (idx < 0)
            continue;

        // Matching face_id and head_id ------------------------------
        bool isMatch(false);
        // check the head_id is already matched
        for (auto it : fID2hIDs) { 
            if (it.second.find(idx) != it.second.end()) {
                isMatch = true;
                break;
            }
        }
        if (isMatch) continue;
            fID2hIDs[fid].insert(idx);
        // ----------------------------------------------------------
    }

    // just for check /////////////////////
    cout << "Lists (faceID -> headIDs)----" << endl;
    for (auto it : fID2hIDs) {
        cout << "[" << member_id[it.first] << "]: " << it.first << " -> "; 
        for (auto it2 : it.second) {
            cout << it2 << " ";
        } cout << endl;
    }
    cout << endl;
    //////////////////////////////////////

    for (auto &obj : camera_objects[frame_count][topView].object_list)
    {
        if (obj.tracking_state != sl::OBJECT_TRACKING_STATE::OK)
            continue;
        Eigen::RowVector2f center(obj.head_position.x, obj.head_position.y);
        Eigen::MatrixXf fusedPosture, fusedBV;
        Eigen::VectorXf fusedConf, fusedConf_BV;

        for (auto it : camera_heads_oneframe)
        {
            int sn = it.first;
            Eigen::MatrixXf heads = it.second;
            Eigen::Index idx, idx0;
            // if the head is too far from the center, skip it
            if ((it.second.rowwise() - center).rowwise().squaredNorm().minCoeff(&idx) > 0.25)
                continue;
            (ref_Heads.rowwise() - heads.row(idx)).rowwise().squaredNorm().minCoeff(&idx0);
            // matching IDs
            if (ref_IDs(idx0) != obj.id)
                continue;
            
            sl::BodyData *body = &camera_bodies[frame_count][sn].body_list[idx];
            if (fusedPosture.rows() == 0)
            {
                fusedPosture = Eigen::MatrixXf::Zero(34, 3);
                fusedConf = Eigen::VectorXf::Zero(34);
                fusedBV = Eigen::MatrixXf::Zero(sl::BODY_34_BONES.size(), 4);
                fusedConf_BV = Eigen::VectorXf::Zero(sl::BODY_34_BONES.size());
            }
            // avg. keypoint
            for (int i = 0; i < 34; i++)
            {
                auto p = body->keypoint[i];
                float conf = body->keypoint_confidence[i];
                fusedPosture.row(i) += Eigen::RowVector3f(p.x, p.y, p.z) * conf;
                fusedConf(i) += conf;
            }
            // avg. BV
            for (size_t i = 0; i < sl::BODY_34_BONES.size(); i++)
            {
                auto bv = body->keypoint[(int)sl::BODY_34_BONES[i].second] - body->keypoint[(int)sl::BODY_34_BONES[i].first];
                float norm = bv.norm();
                bv = bv / norm;
                double conf = (body->keypoint_confidence[(int)sl::BODY_34_BONES[i].first] + body->keypoint_confidence[(int)sl::BODY_34_BONES[i].second]) * 0.5;
                fusedBV.row(i) += Eigen::RowVector4f(bv.x, bv.y, bv.z, norm) * conf * conf;
                fusedConf_BV(i) += conf * conf;
            }
        }

        if (fusedPosture.rows() > 0)
        {
            // avg. kpt
            fusedBody.is_tracked = true;
            sl::BodyData body;
            body.id = obj.id;
            body.tracking_state = sl::OBJECT_TRACKING_STATE::OK;
            std::vector<sl::float3> keypoints;
            
            bool isTrack(false);
            for (auto it: fID2hIDs) {
                if (it.second.find(obj.id) != it.second.end()) {
                    body.id = it.first;
                    isTrack = true;
                    break;
                }
                else {
                    isTrack = false;
                }
            }
            // if (!isTrack) {
            //     continue;
            // }
            Eigen::MatrixXf pose = Eigen::MatrixXf::Zero(35, 3);
            pose.row(34) = Eigen::RowVector3f(center(0)*100, center(1)*100, 0.);
            if (isTrack) ofs_BV << "p " << body.id << " " << center(0) * 100 << " " << center(1) * 100 << endl;
            else         ofs_BV << "d " << body.id << " " << center(0) * 100 << " " << center(1) * 100 << endl;

            for (int i = 0; i < 34; i++)
            {
                Eigen::RowVector3f p = fusedPosture.row(i) / fusedConf(i);
                keypoints.push_back(sl::float3(p(0), p(1), p(2)));
                p *= 100; // meter->cm
            }
            body.keypoint = keypoints;
            fusedBody.body_list.push_back(body);

            // avg. bv
            fusedBody_BV.is_tracked = true;
            Eigen::MatrixXf kpt = Eigen::MatrixXf::Zero(34, 3);
            kpt.row(0) = Eigen::RowVector3f(keypoints[0].x, keypoints[0].y, keypoints[0].z);
            for (int i = 0; i < sl::BODY_34_BONES.size(); i++)
            {
                Eigen::RowVector4f b = fusedBV.row(i) / fusedConf_BV(i);
                kpt.row((int)sl::BODY_34_BONES[i].second) = kpt.row((int)sl::BODY_34_BONES[i].first) + b.leftCols(3).normalized() * b(3);
            }
            if (isnan(kpt(30, 0) + kpt(28, 0)))
            {
                // left ear/eye
                kpt.row(28) = Eigen::RowVector3f(keypoints[28].x, keypoints[28].y, keypoints[28].z);
                kpt.row(29) = Eigen::RowVector3f(keypoints[29].x, keypoints[29].y, keypoints[29].z);
                // right ear/eye
                kpt.row(30) = Eigen::RowVector3f(keypoints[30].x, keypoints[30].y, keypoints[30].z);
                kpt.row(31) = Eigen::RowVector3f(keypoints[31].x, keypoints[31].y, keypoints[31].z);
            }
            keypoints.clear();
            for (int i = 0; i < 34; i++)
            {
                keypoints.push_back(sl::float3(kpt(i, 0), kpt(i, 1), kpt(i, 2)));
                pose.row(i) = kpt.row(i) * 100;
                ofs_BV << kpt.row(i) * 100 << endl;
            }
            body.keypoint = keypoints;
            fusedBody_BV.body_list.push_back(body);

            if (isTrack) phantom_pose[body.id] = pose;
            else         dummy_pose[body.id+member_id.size()] = pose;
        }
        else
        {
            bool isTrack(false);
            for (auto it: fID2hIDs) {
                if (it.second.find(obj.id) != it.second.end()) {
                    obj.id = it.first;
                    isTrack = true;
                    break;
                }
                else {
                    isTrack = false;
                }
            }
            if(isTrack) {
                ofs_BV << "p " << -obj.id << " " << center(0) * 100 << " " << center(1) * 100 << endl;
                phantom_pose[obj.id].row(0) = Eigen::RowVector2f(center(0)*100, center(1)*100);
            }
            else {
                ofs_BV << "d " << -obj.id << " " << center(0) * 100 << " " << center(1) * 100 << endl;
                dummy_pose[obj.id+member_id.size()].row(0) = Eigen::RowVector2f(center(0)*100, center(1)*100);
            }        
        }
    }
}


void MultiPoseTracker::convertAffineEigen2SL(const Eigen::Affine3f &inp, sl::Transform &out)
{
    Eigen::Vector3f T = inp.translation();
    Eigen::Quaternionf Q = Eigen::Quaternionf(inp.rotation());
    Q.normalize();
    sl::float3 T_sl(T.x(), T.y(), T.z());
    sl::float4 Q_sl(Q.x(), Q.y(), Q.z(), Q.w());
    out.setTranslation(T_sl);
    out.setOrientation(Q_sl);
}

void MultiPoseTracker::transformObjects(const Eigen::Affine3f &cam_pose, sl::Objects &camera_objects)
{
    for (auto &head : camera_objects.object_list)
    {
        for (auto &bb : head.head_bounding_box)
        {
            Eigen::Vector3f bb_tf = cam_pose * Eigen::Vector3f(bb.x, bb.y, bb.z + 0.2);
            bb.x = bb_tf.x();
            bb.y = bb_tf.y();
            bb.z = bb_tf.z();
        }
        Eigen::Vector3f head_tf = cam_pose * Eigen::Vector3f(head.head_position.x, head.head_position.y, head.head_position.z);
        head.head_position.x = head_tf.x();
        head.head_position.y = head_tf.y();
        head.head_position.z = head_tf.z();
    }
}
void MultiPoseTracker::transformBodies(const Eigen::Affine3f &cam_pose, sl::Bodies &camera_bodies)
{
    for (auto &body : camera_bodies.body_list)
    {
        for (auto &kpt : body.keypoint)
        {
            Eigen::Vector3f kpt_tf = cam_pose * Eigen::Vector3f(kpt.x, kpt.y, kpt.z);
            kpt.x = kpt_tf.x();
            kpt.y = kpt_tf.y();
            kpt.z = kpt_tf.z();
        }
        for (auto &pt : body.bounding_box)
        {
            Eigen::Vector3f pt_tf = cam_pose * Eigen::Vector3f(pt.x, pt.y, pt.z);
            pt.x = pt_tf.x();
            pt.y = pt_tf.y();
            pt.z = pt_tf.z();
        }
    }
}

float MultiPoseTracker::Point2LineDistance2(const Eigen::Vector3f& zed_head,
                                            const Eigen::Vector3f& face_center,
                                            const Eigen::Affine3f& wcam_a)
{
    // 1. zed head data -> cm
    Eigen::Vector3f zed_head_cm = zed_head * 100;
    // 2. face center data (100 cm dist image, optional)
    // Eigen::Vector3f face_center_s = face_center * (100/face_center.z());
    // 3. transform face_center to world coordinate
    Eigen::Affine3f wcam_a_cm = wcam_a;
    wcam_a_cm.translation() *= 100;
    Eigen::Vector3f face_center_tf = wcam_a_cm * face_center;
    // 4. wcam_position
    Eigen::Vector3f wcam_pos = wcam_a_cm.translation();
    // 5. line vector (l,m,n)
    Eigen::Vector3f line(face_center_tf.x() - wcam_pos.x(), 
                         face_center_tf.y() - wcam_pos.y(),
                         face_center_tf.z() - wcam_pos.z());
    line.normalize();
    float l = line.x();
    float m = line.y();
    float n = line.z();
    // 6. calculate parameter t
    float t = l * (zed_head_cm.x() - wcam_pos.x()) + m * (zed_head_cm.y() - wcam_pos.y()) + n * (zed_head_cm.z() - wcam_pos.z())
            / (l * l + m * m + n * n);
    // 7. calculate point to line squared distance
    Eigen::Vector3f p(wcam_pos.x() + l * t, wcam_pos.y() + m * t, wcam_pos.z() + n * t);
    float dist2 = (p - zed_head_cm).squaredNorm();

    // cout check
    // cout << "--- Point2LineDistance2 ---" << endl;
    // cout << "zed_head (mm): " << zed_head_cm.transpose() * 10 << endl;
    // cout << "face_center (px): " << face_center.transpose() << endl;
    // // cout << "face_center_s (cm): " << face_center_s.transpose() << endl;
    // cout << "face_center_tf (px): " << face_center_tf.transpose() << endl;
    // cout << "wcam_pos (mm): " << wcam_pos.transpose() * 10 << endl;
    // cout << "wcam_xyzw: " << Eigen::Quaternionf(wcam_a.rotation()).coeffs().transpose() << endl;
    // cout << "l,m,n: " << l << ", " << m << ", " << n << endl;
    // cout << "t: " << t << endl;
    // cout << "p: " << p.transpose() << endl;
    // cout << "dist2: " << dist2 << endl;
    // cout << endl;


    return dist2;
}
