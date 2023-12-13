#include "MultiPhantomDeformer.hh"

MultiPhantomDeformer::MultiPhantomDeformer(const vector<string>& _member_phantom)
{
    cout << "MultiPhantomDeformer ==========================" << endl;
    member_phantom = _member_phantom;
    LoadPhantom();
    LoadDummy();
    tgf2zed = {0, 1, 2, 4, 5, 6, 7, 3, 11, 12, 13, 14, 18, 19, 20, 22, 23, 24, 9, 16, 21, 25, 26};
    for (int id=0; id<member_phantom.size(); id++) {
        C1[id] = Eigen::MatrixXd::Constant(28, 3, -1001);
        C1_prev[id] = Eigen::MatrixXd::Constant(28, 3, -1001);
    }
    for (int id=member_phantom.size(); id<member_phantom.size()*2; id++) {
        C1_dummy[id] = Eigen::MatrixXd::Constant(28, 3, -1001);
        C1_prev_dummy[id] = Eigen::MatrixXd::Constant(28, 3, -1001);
    }
}

MultiPhantomDeformer::~MultiPhantomDeformer()
{
    for (auto& p: pc) {
        delete p.second;
    }
}

void MultiPhantomDeformer::RunPhantom(const map<int, Eigen::MatrixXf>& phantom_pose)
{
    view_ids.clear();
    for (auto p: phantom_pose) {
        int id = p.first;
        Eigen::MatrixXd pose = p.second.cast<double>();
            
        if (pose.rows() == 35) {
            for (int i=0; i < 23; i++) {
                C1[id].row(i) = pose.row(tgf2zed[i]);
            }
            C1[id].row(23) = pose.row(28); // eyeL
            C1[id].row(24) = pose.row(30); // eyeR
            C1[id].row(25) = pose.row(29); // earL
            C1[id].row(26) = pose.row(31); // earR
            C1[id].row(27) = pose.row(34); // head XY
            pc[id]->set_C1(C1[id]);
            view_ids.insert(id);
        }
        else if ( pose.rows() == 1 && C1_prev[id](0,2) > -1000 )
        {
            C1[id] = C1_prev[id];
            C1[id].row(27) = pose.row(0);
            pc[id]->set_C1(C1[id]);                
            view_ids.insert(id);
        }
        else {
            continue;
        }

        if ( !pc[id]->ROOT_JointsCorrection() )
                continue;
        pc[id]->CHEST_JointsCorrection();
        pc[id]->LEGS_JointsCorrection();
        pc[id]->ARMS_JointsCorrection();
        pc[id]->LIMBS_JointsCorrection();
        pc[id]->HEAD_xyCorrection();
        C1[id] = pc[id]->get_C1();
    }
        
    // prevent arm overlapping
    for (auto i: C1) {
        int id = i.first;
        if (C1[id](0,2) < -1000) continue;
        else if (view_ids.find(id) == view_ids.end()) continue;
        for (auto j: C1) {
            if (C1[j.first](0,2) < -1000) continue;
            else if (view_ids.find(j.first) == view_ids.end()) continue;
            else if (id == j.first) continue;
            double L0 = ( C1[id].row(6) - C1[j.first].row(0) ).squaredNorm();
            double L1 = ( C1[id].row(6) - C1[j.first].row(1) ).squaredNorm();
            double L2 = ( C1[id].row(6) - C1[j.first].row(2) ).squaredNorm();
            double L3 = ( C1[id].row(5) - C1[j.first].row(0) ).squaredNorm();
            double L4 = ( C1[id].row(5) - C1[j.first].row(1) ).squaredNorm();
            double L5 = ( C1[id].row(5) - C1[j.first].row(2) ).squaredNorm();
            if (L0 < 900 || L1 < 900 || L2 < 900 || L3 < 900 || L4 < 900 || L5 < 900) { // 30 cm
                pc[id]->use_default_Larm();
            }
            double R0 = ( C1[id].row(11) - C1[j.first].row(0) ).squaredNorm();
            double R1 = ( C1[id].row(11) - C1[j.first].row(1) ).squaredNorm();
            double R2 = ( C1[id].row(11) - C1[j.first].row(2) ).squaredNorm();
            double R3 = ( C1[id].row(10) - C1[j.first].row(0) ).squaredNorm();
            double R4 = ( C1[id].row(10) - C1[j.first].row(1) ).squaredNorm();
            double R5 = ( C1[id].row(10) - C1[j.first].row(2) ).squaredNorm();
            if (R0 < 900 || R1 < 900 || R2 < 900 || R3 < 900 || R4 < 900 || R5 < 900) { // 30 cm
                pc[id]->use_default_Rarm();
            }
        }
        C1_prev[id] = C1[id];
        pc[id]->set_C1_prev(C1_prev[id]);
    }
    
    for (auto p: phantom_pose) {
        int id = p.first;
        if (view_ids.find(id) == view_ids.end()) continue;            
        RotationList vQ = pc[id]->get_vQ();
        vector<Eigen::Vector3d> vT(BE[id].rows(), Eigen::Vector3d::Zero());
        Eigen::MatrixXd C_new(23, 3);
        C_new.row(0) = C1[id].row(0);
        
        for (int i=1; i<BE[id].rows(); i++) {
            C_new.row( BE[id](i,1) ) = C_new.row(BE[id](i,0)) + ( vQ[i].matrix() * ( C[id].row(BE[id](i,1)) - C[id].row(BE[id](i,0)) ).transpose() ).transpose();
            Eigen::Affine3d a = Eigen::Translation3d(C_new.row(BE[id](i,0))) * vQ[i].matrix() * Eigen::Translation3d(-C[id].row(BE[id](i,0)));
            vT[i] = a.translation();
        }
        C_new.row( BE[id](0,1) ) = C_new.row(BE[id](0,0)) + ( vQ[0].matrix() * (C[id].row(BE[id](0,1)) - C[id].row(BE[id](0,0))).transpose() ).transpose();
        C_NEW[id] = C_new;
        Eigen::Affine3d a = Eigen::Translation3d(C_new.row(BE[id](0,0))) * vQ[0].matrix() * Eigen::Translation3d(-C[id].row(BE[id](0,0)));
        vT[0] = a.translation();
        myDqs(Vo[id], weightMap[id], vQ, vT, U[id]);
    }

}


void MultiPhantomDeformer::LoadPhantom()
{
    for (int id=0; id<member_phantom.size(); id++) {
        string VIR_Data_PATH = string(getenv("VIR_Data"));
        VIR_Data_PATH += "phantoms/" + member_phantom[id] + "/";
        string tgf_file = VIR_Data_PATH + member_phantom[id] + "_DS.tgf";
        string ply_file = VIR_Data_PATH + member_phantom[id] + "_DS.ply";
        string Wo_file = VIR_Data_PATH + member_phantom[id] + "_DS_Wo.dmat";
        string ELE_file = VIR_Data_PATH + member_phantom[id] + "_DS_ELE.dmat";
        string NODE_file = VIR_Data_PATH + member_phantom[id] + "_DS_NODE.dmat";
        cout << "VIR_Data_PATH: " << VIR_Data_PATH << endl;
        cout << "Load " << tgf_file << endl;
        cout << "Load " << ply_file << endl;
        cout << "Load " << Wo_file << endl;
        cout << "Load " << ELE_file << endl;
        cout << "Load " << NODE_file << endl;
        if ( !igl::readTGF( tgf_file, C[id], BE[id]) ) {
            cout << "There is no " << tgf_file << endl;
            return;
        }
        igl::directed_edge_parents(BE[id], P[id]);
            
        // read SHELL
        igl::readPLY(ply_file, Vo[id], Fo[id]);
        int numVo = Vo[id].rows();
        
        if ( !igl::readDMAT(Wo_file, Wo[id]) ) {
            // TETRAHEDRALIZATION
            Vo[id].conservativeResize(Vo[id].rows() + C[id].rows(), 3);
            Vo[id].bottomRows( C[id].rows() ) = C[id];
            for (int i=0; i<BE[id].rows();i ++)
            {
                int num = ( C[id].row(BE[id](i,0)) - C[id].row(BE[id](i,1)) ).norm();
                Eigen::RowVector3d itvl = ( C[id].row(BE[id](i,0)) - C[id].row(BE[id](i,1)) ) / (double)num;
                Eigen::MatrixXd boneP(num-1, 3);
                for (int n=1; n<num; n++)
                    boneP.row(n-1) = C[id].row(BE[id](i,1)) + itvl * n;
                Vo[id].conservativeResize(Vo[id].rows() + boneP.rows(), 3);
                Vo[id].bottomRows( boneP.rows() ) = boneP;
            }
            igl::copyleft::tetgen::tetrahedralize(Vo[id], Fo[id], "qp/0.0001YT0.000000001", VT[id], TT[id], FT[id]);
            FT[id].resize(0, 0);
            igl::writeDMAT(ELE_file, TT[id], false);
            igl::writeDMAT(NODE_file, VT[id], false);
            Vo[id] = Vo[id].topRows(numVo);

            // BBW - BONE
            Eigen::MatrixXd bc;
            Eigen::VectorXi b;
            if ( igl::boundary_conditions(VT[id], TT[id], C[id], Eigen::VectorXi(), BE[id], Eigen::MatrixXi(), b, bc) )
                cout << "boundary condition was generated for " << b.rows() << " vertices." << endl;
            else {
                cout << "boundary condition was not generated." << endl;
                return;
            }

            igl::BBWData bbw_data;
            bbw_data.active_set_params.max_iter = 10;
            bbw_data.verbosity = 2;
            igl::normalize_row_sums(bc, bc);
            igl::bbw(VT[id], TT[id], b, bc, bbw_data, Wo[id]);
            Wo[id] = Wo[id].topRows(numVo);
            igl::writeDMAT(Wo_file, Wo[id], false);
        }

        for (int i=0; i<Wo[id].rows(); i++) {
            map<int, double> w;
            for (int j=0; j<Wo[id].cols(); j++)
                w[j] = Wo[id](i, j);
            weightMap[id].push_back(w);
        }
        Eigen::RowVector3d isocenter(151.3, 28.9, 90.1);
        pc.insert( make_pair(id, new PoseCorrector(C[id], BE[id], P[id], isocenter)) );
    }    
}

void MultiPhantomDeformer::LoadDummy()
{
    string VIR_Data_PATH = string(getenv("VIR_Data"));
    VIR_Data_PATH += "phantoms/amMRCP_DS/amMRCP_DS";
    string tgf_file = VIR_Data_PATH + ".tgf";
    string ply_file = VIR_Data_PATH + ".ply";
    string Wo_file = VIR_Data_PATH + "_Wo.dmat";
    string ELE_file = VIR_Data_PATH + "_ELE.dmat";
    string NODE_file = VIR_Data_PATH + "_NODE.dmat";
    cout << "VIR_Data_PATH: " << VIR_Data_PATH << endl;
    cout << "Load " << tgf_file << endl;
    cout << "Load " << ply_file << endl;
    cout << "Load " << Wo_file << endl;
    cout << "Load " << ELE_file << endl;
    cout << "Load " << NODE_file << endl;

    for (int id=member_phantom.size(); id<member_phantom.size()*2; id++) {
        if ( !igl::readTGF( tgf_file, C_dummy[id], BE_dummy[id]) ) {
            cout << "There is no " << tgf_file << endl;
            return;
        }
        igl::directed_edge_parents(BE_dummy[id], P_dummy[id]);
            
        // read SHELL
        igl::readPLY(ply_file, Vo_dummy[id], Fo_dummy[id]);
        int numVo = Vo_dummy[id].rows();
        
        if ( !igl::readDMAT(Wo_file, Wo_dummy[id]) ) {
            // TETRAHEDRALIZATION
            Vo_dummy[id].conservativeResize(Vo_dummy[id].rows() + C_dummy[id].rows(), 3);
            Vo_dummy[id].bottomRows( C_dummy[id].rows() ) = C_dummy[id];
            for (int i=0; i<BE_dummy[id].rows();i ++)
            {
                int num = ( C_dummy[id].row(BE_dummy[id](i,0)) - C_dummy[id].row(BE_dummy[id](i,1)) ).norm();
                Eigen::RowVector3d itvl = ( C_dummy[id].row(BE_dummy[id](i,0)) - C_dummy[id].row(BE_dummy[id](i,1)) ) / (double)num;
                Eigen::MatrixXd boneP(num-1, 3);
                for (int n=1; n<num; n++)
                    boneP.row(n-1) = C_dummy[id].row(BE_dummy[id](i,1)) + itvl * n;
                Vo_dummy[id].conservativeResize(Vo_dummy[id].rows() + boneP.rows(), 3);
                Vo_dummy[id].bottomRows( boneP.rows() ) = boneP;
            }
            igl::copyleft::tetgen::tetrahedralize(Vo_dummy[id], Fo_dummy[id], "qp/0.0001YT0.000000001", VT_dummy[id], TT_dummy[id], FT_dummy[id]);
            FT_dummy[id].resize(0, 0);
            igl::writeDMAT(ELE_file, TT_dummy[id], false);
            igl::writeDMAT(NODE_file, VT_dummy[id], false);
            Vo_dummy[id] = Vo_dummy[id].topRows(numVo);

            // BBW - BONE
            Eigen::MatrixXd bc;
            Eigen::VectorXi b;
            if ( igl::boundary_conditions(VT_dummy[id], TT_dummy[id], C_dummy[id], Eigen::VectorXi(), BE_dummy[id], Eigen::MatrixXi(), b, bc) )
                cout << "boundary condition was generated for " << b.rows() << " vertices." << endl;
            else {
                cout << "boundary condition was not generated." << endl;
                return;
            }

            igl::BBWData bbw_data;
            bbw_data.active_set_params.max_iter = 10;
            bbw_data.verbosity = 2;
            igl::normalize_row_sums(bc, bc);
            igl::bbw(VT_dummy[id], TT_dummy[id], b, bc, bbw_data, Wo_dummy[id]);
            Wo_dummy[id] = Wo_dummy[id].topRows(numVo);
            igl::writeDMAT(Wo_file, Wo_dummy[id], false);
        }

        for (int i=0; i<Wo_dummy[id].rows(); i++) {
            map<int, double> w;
            for (int j=0; j<Wo_dummy[id].cols(); j++)
                w[j] = Wo_dummy[id](i, j);
            weightMap_dummy[id].push_back(w);
        }
        Eigen::RowVector3d isocenter(151.3, 28.9, 90.1);
        pc_dummy.insert( make_pair(id, new PoseCorrector(C_dummy[id], BE_dummy[id], P_dummy[id], isocenter)) );
    } 
}

void MultiPhantomDeformer::RunDummy(const map<int, Eigen::MatrixXf>& dummy_pose)
{
    view_ids_dummy.clear();
    for (auto d: dummy_pose) {
        int id = d.first;
        if (id >= member_phantom.size()*2) {
            id = id%member_phantom.size() + member_phantom.size();
        }
        Eigen::MatrixXd pose = d.second.cast<double>();
        if (pose.rows() == 35) {
            for (int i=0; i < 23; i++) {
                C1_dummy[id].row(i) = pose.row(tgf2zed[i]);
            }
            C1_dummy[id].row(23) = pose.row(28); // eyeL
            C1_dummy[id].row(24) = pose.row(30); // eyeR
            C1_dummy[id].row(25) = pose.row(29); // earL
            C1_dummy[id].row(26) = pose.row(31); // earR
            C1_dummy[id].row(27) = pose.row(34); // head XY
            pc_dummy[id]->set_C1(C1_dummy[id]);
            view_ids_dummy.insert(id);
        }
        else if ( pose.rows() == 1 && C1_prev_dummy[id](0,2) > -1000 )
        {
            C1_dummy[id] = C1_prev_dummy[id];
            C1_dummy[id].row(27) = pose.row(0);
            pc_dummy[id]->set_C1(C1_dummy[id]);
            view_ids_dummy.insert(id);
        }
        else {
            continue;
        }

        if ( !pc_dummy[id]->ROOT_JointsCorrection() )
                continue;
        pc_dummy[id]->CHEST_JointsCorrection();
        pc_dummy[id]->LEGS_JointsCorrection();
        pc_dummy[id]->ARMS_JointsCorrection();
        pc_dummy[id]->LIMBS_JointsCorrection();
        pc_dummy[id]->HEAD_xyCorrection();
        C1_dummy[id] = pc_dummy[id]->get_C1();
    }
    // prevent arm overlapping
    for (auto i: C1_dummy) {
        int id = i.first;
        // if (id >= member_phantom.size()*2) {
        //     id = id%member_phantom.size() + member_phantom.size();
        // }
        if (C1_dummy[id](0,2) < -1000) continue;
        else if (view_ids_dummy.find(id) == view_ids_dummy.end()) continue;
        for (auto j: C1_dummy) {
            if (C1_dummy[j.first](0,2) < -1000) continue;
            else if (view_ids_dummy.find(j.first) == view_ids_dummy.end()) continue;
            else if (id == j.first) continue;
            double L0 = ( C1_dummy[id].row(6) - C1_dummy[j.first].row(0) ).squaredNorm();
            double L1 = ( C1_dummy[id].row(6) - C1_dummy[j.first].row(1) ).squaredNorm();
            double L2 = ( C1_dummy[id].row(6) - C1_dummy[j.first].row(2) ).squaredNorm();
            double L3 = ( C1_dummy[id].row(5) - C1_dummy[j.first].row(0) ).squaredNorm();
            double L4 = ( C1_dummy[id].row(5) - C1_dummy[j.first].row(1) ).squaredNorm();
            double L5 = ( C1_dummy[id].row(5) - C1_dummy[j.first].row(2) ).squaredNorm();
            if (L0 < 900 || L1 < 900 || L2 < 900 || L3 < 900 || L4 < 900 || L5 < 900) { // 30 cm
                pc_dummy[id]->use_default_Larm();
            }
            double R0 = ( C1_dummy[id].row(11) - C1_dummy[j.first].row(0) ).squaredNorm();
            double R1 = ( C1_dummy[id].row(11) - C1_dummy[j.first].row(1) ).squaredNorm();
            double R2 = ( C1_dummy[id].row(11) - C1_dummy[j.first].row(2) ).squaredNorm();
            double R3 = ( C1_dummy[id].row(10) - C1_dummy[j.first].row(0) ).squaredNorm();
            double R4 = ( C1_dummy[id].row(10) - C1_dummy[j.first].row(1) ).squaredNorm();
            double R5 = ( C1_dummy[id].row(10) - C1_dummy[j.first].row(2) ).squaredNorm();
            if (R0 < 900 || R1 < 900 || R2 < 900 || R3 < 900 || R4 < 900 || R5 < 900) { // 30 cm
                pc_dummy[id]->use_default_Rarm();
            }
        }

        if (C1_dummy[id].row(0) == C1_prev_dummy[id].row(0)) {
            view_ids_dummy.erase(id);
        }

        C1_prev_dummy[id] = C1_dummy[id];
        pc_dummy[id]->set_C1_prev(C1_prev_dummy[id]);
    }
    
    for (auto d: dummy_pose) {
        int id = d.first;
        if (id >= member_phantom.size()*2) {
            id = id%member_phantom.size() + member_phantom.size();
        }
        if (view_ids_dummy.find(id) == view_ids_dummy.end()) continue;            
        RotationList vQ = pc_dummy[id]->get_vQ();
        vector<Eigen::Vector3d> vT(BE_dummy[id].rows(), Eigen::Vector3d::Zero());
        Eigen::MatrixXd C_new(23, 3);
        C_new.row(0) = C1_dummy[id].row(0);
        
        for (int i=1; i<BE_dummy[id].rows(); i++) {
            C_new.row( BE_dummy[id](i,1) ) = C_new.row(BE_dummy[id](i,0)) + ( vQ[i].matrix() * ( C_dummy[id].row(BE_dummy[id](i,1)) - C_dummy[id].row(BE_dummy[id](i,0)) ).transpose() ).transpose();
            Eigen::Affine3d a = Eigen::Translation3d(C_new.row(BE_dummy[id](i,0))) * vQ[i].matrix() * Eigen::Translation3d(-C_dummy[id].row(BE_dummy[id](i,0)));
            vT[i] = a.translation();
        }
        C_new.row( BE_dummy[id](0,1) ) = C_new.row(BE_dummy[id](0,0)) + ( vQ[0].matrix() * (C_dummy[id].row(BE_dummy[id](0,1)) - C_dummy[id].row(BE_dummy[id](0,0))).transpose() ).transpose();
        C_NEW_dummy[id] = C_new;
        Eigen::Affine3d a = Eigen::Translation3d(C_new.row(BE_dummy[id](0,0))) * vQ[0].matrix() * Eigen::Translation3d(-C_dummy[id].row(BE_dummy[id](0,0)));
        vT[0] = a.translation();
        myDqs(Vo_dummy[id], weightMap_dummy[id], vQ, vT, U_dummy[id]);
    }
}