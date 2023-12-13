#include "PoseCorrector.hh"

PoseCorrector::PoseCorrector(const Eigen::MatrixXd& _C, 
                             const Eigen::MatrixXi& _BE,
                             const Eigen::VectorXi& _P,
                             const Eigen::RowVector3d& _isocenter)
{    
    C = _C;
    BE = _BE;
    P = _P;
    isocenter = _isocenter;

    // save bone length and vector
    bL.resize(BE.rows());
    nbV.resize(BE.rows(), 3);
    bV.resize(BE.rows(), 3);
    for (int i=0; i<BE.rows(); i++) {
        bL(i) = ( C.row(BE(i,0)) - C.row(BE(i,1)) ).norm();
        nbV.row(i) = ( C.row(BE(i,1)) - C.row(BE(i,0)) ).normalized();
        bV.row(i) = C.row(BE(i,1)) - C.row(BE(i,0));
    }

    // root
    root0.col(0) = ( C.row(12) - C.row(15) ).normalized();
    root0.col(2) = nbV.row(1);
    root0.col(1) = root0.col(2).cross(root0.col(0)).normalized();
    root0.col(0) = root0.col(1).cross(root0.col(2)).normalized();
    // clav
    clav0.col(0) = ( C.row(3) - C.row(8) ).normalized();
    clav0.col(2) = nbV.row(8);
    clav0.col(1) = clav0.col(2).cross(clav0.col(0)).normalized();
    clav0.col(0) = clav0.col(1).cross(clav0.col(2)).normalized();
    // chest
    chest0.col(0) = (root0.col(0) + clav0.col(0)).normalized();
    chest0.col(2) = nbV.row(2);
    chest0.col(1) = chest0.col(2).cross(chest0.col(0)).normalized();
    chest0.col(0) = chest0.col(1).cross(chest0.col(2)).normalized();
    // head
    head0.col(0) = root0.col(0);
    head0.col(2) = (C.row(22)-C.row(7)-Eigen::RowVector3d::UnitZ()*8.).normalized(); //adjusted
    // head0.col(2) = (C.row(22)-C.row(7)).normalized(); //adjusted
    head0.col(1) = head0.col(2).cross(head0.col(0)).normalized();
    head0.col(2) = head0.col(0).cross(head0.col(1)).normalized();

    root_trans = Eigen::RowVector3d::Zero();
    clav_trans = Eigen::RowVector3d::Zero();    
}

PoseCorrector::~PoseCorrector()
{

}

bool PoseCorrector::ROOT_JointsCorrection()
{
    // Check root
    if ( C1(0,2) < -1000 ) {
        if ( C1_prev.rows() == 28 && C1_prev.row(0)(2) > -1000) { // use previous frame data
            C1.row(0) = C1_prev.row(0);
        } 
        else {
            // cout << id << "'s pelvis joint is not detected" << endl;
        }
        return false;
    }

    // calculate root translation (present - previous frames)
    if ( C1_prev.rows() == 28 && C1_prev.row(0)(2) > -1000) {
        root_trans = C1.row(0) - C1_prev.row(0);
    }

    // spineN
    if ( C1(1,2) < -1000) {
        if ( C1_prev.rows() == 28 && C1_prev.row(1)(2) > -1000) {
            C1.row(1) = C1_prev.row(1) + root_trans;
        } else {
            C1.row(1) = C1.row(0) + bL(1) * Eigen::RowVector3d::UnitZ();
        }
    }

    // hipL, hipR
    if ( C1(12,2) < -1000 && C1(15,2) > -1000 ) { // hipL not detected, hipR detected
        C1.row(12) = C1.row(0) + (C1.row(0) - C1.row(15));
    }
    else if ( C1(12,2) > -1000 && C1(15,2) < -1000 ) { // hipL detected, hipR not detected
        C1.row(15) = C1.row(0) + (C1.row(0) - C1.row(12));
    }
    else if ( C1(12,2) < -1000 && C1(15,2) < -1000) { // both not detected
        if ( C1_prev.rows() == 28 && C1_prev.row(12)(2) > -1000 && C1_prev.row(15)(2) > -1000) { // use previous frame data
            C1.row(12) = C1_prev.row(12) + root_trans;
            C1.row(15) = C1_prev.row(15) + root_trans;
        }
        else {
            if (C1(3,2) > -1000 && C1(8,2) > -1000) {
                Eigen::RowVector3d clavR2clavL = C1.row(3) - C1.row(8);
                C1.row(12) = C1.row(0) + clavR2clavL.normalized() * bL(14);
                C1.row(15) = C1.row(0) - clavR2clavL.normalized() * bL(18);
            }
            else {
                // cout << id << "'s hip joints are not detected" << endl;
                return false;
            }
        }
    }
    
    // Adjust hip joints
    Eigen::RowVector3d adjust = C1.row(0) - (C1.row(12) + C1.row(15)) * 0.5;
    C1.row(12) += adjust;
    C1.row(15) += adjust;
    C1(12,2) = C1(0,2);
    C1(15,2) = C1(0,2);

    // Get Root Rotation
    root1.col(0) = ( C1.row(12) - C1.row(15) ).normalized(); // pelvis line
    root1.col(2) = ( C1.row(1) - C1.row(0) ).normalized();
    root1.col(1) = root1.col(2).cross(root1.col(0)).normalized();
    root1.col(0) = root1.col(1).cross(root1.col(2)).normalized();
    rootQ = Eigen::Quaterniond(root1 * root0.transpose()); rootQ.normalize();

    // Set vQ
    vQ = RotationList(BE.rows());
    for (int i=0;i<BE.rows();i++) {
        vQ[i] = rootQ;
    }
    vQ_corr = vQ;

    return true;
}

void PoseCorrector::CHEST_JointsCorrection()
{
    // spineC
    if ( C1(2,2) < -1000) {
        if ( C1_prev.rows() == 28 && C1_prev.row(2)(2) > -1000) {
            C1.row(2) = C1_prev.row(2) + root_trans;
        } else {
            C1.row(2) = C1.row(1) + ( vQ[1].matrix() * bV.row(2).transpose() ).transpose();
        }
    }
    
    // @@ Clavicle joints
    // Matching //////////////////////////
    Eigen::RowVector3d hipR2hipL = C1.row(12) - C1.row(15);
    if ( C1(3,2) < -1000 && C1(8,2) > -1000 ) { // clavL not detected, clavR detected
        Eigen::RowVector3d spineC2clavR = C1.row(8) - C1.row(2);
        double projLen = -hipR2hipL.dot(spineC2clavR) / hipR2hipL.norm();
        C1.row(3) = C1.row(8) + hipR2hipL.normalized() * projLen * 2;
    }
    else if ( C1(3,2) > -1000 && C1(8,2) < -1000 ) { // clavL detected, clavR not detected
        Eigen::RowVector3d spineC2clavL = C1.row(3) - C1.row(2);
        double projLen = hipR2hipL.dot(spineC2clavL) / hipR2hipL.norm();
        C1.row(8) = C1.row(3) - hipR2hipL.normalized() * projLen * 2;
    }
    else if ( C1(3,2) < -1000 && C1(8,2) < -1000) { // both not detected
        if ( C1_prev.rows() == 28 && C1_prev.row(3)(2) > -1000 && C1_prev.row(8)(2)) { // use previous frame data
            C1.row(3) = C1_prev.row(3) + root_trans;
            C1.row(8) = C1_prev.row(8) + root_trans;
        }
        else {
            C1.row(3) = C1.row(2) + ( vQ[2].matrix() * bV.row(3).transpose() ).transpose();
            C1.row(8) = C1.row(2) + ( vQ[2].matrix() * bV.row(9).transpose() ).transpose();
        }
    }

    // Neck
    C1.row(7) = ( C1.row(3) + C1.row(8) ) * 0.5;

    // calculate clav/chest translation (present - previous frames)
    if ( C1_prev.rows() == 28 && C1_prev.row(7)(2) > -1000 )
        clav_trans = C1.row(7) - C1_prev.row(7);

    // Get Clav
    clav1.col(0) = ( C1.row(3) - C1.row(8) ).normalized();
    clav1.col(2) = ( C1.row(7) - C1.row(2) ).normalized();
    clav1.col(1) = clav1.col(2).cross(clav1.col(0)).normalized();
    clav1.col(0) = clav1.col(1).cross(clav1.col(2)).normalized();
    clavQ = Eigen::Quaterniond(clav1 * clav0.transpose()); clavQ.normalize();

    // Get Chest
    chest1.col(0) = (root1.col(0) + clav1.col(0)).normalized();
    chest1.col(2) = (C1.row(2) - C1.row(1)).normalized();
    chest1.col(1) = chest1.col(2).cross(chest1.col(0)).normalized();
    chest1.col(0) = chest1.col(1).cross(chest1.col(2)).normalized();
    chestQ = Eigen::Quaterniond(chest1 * chest0.transpose()); chestQ.normalize();
    
    // Head
    if ( C1(22,2) < -1000 ) {
        if ( C1_prev.rows() == 28 && C1_prev(22,2) > -1000 ) {
            C1.row(22) = C1_prev.row(22) + clav_trans;
        } else {
            C1.row(22) = C1.row(7) + ( clavQ.matrix() * bV.row(0).transpose() ).transpose();
        }
    }
    // Eigen::RowVector3d neck2head = (C1.row(22) - C1.row(7)).normalized();
    // if (projectOntoPlane(neck2head, clav1.col(2)).dot(clav1.col(1)) > 0) { // Back
    //     // C1.row(22) = C1.row(7) + ( clavQ.matrix() * bV.row(0).transpose() ).transpose();
    //     // cout << "back" << endl;
    //     C1.row(22) = C1.row(7) + neck2head * bL(0);
    // }
    // else { // Front
    //     proj_zPlane = projectOntoPlane(neck2head, clav1.col(2)).normalized();
    //     angle = acos( proj_zPlane.dot(clav1.col(0)) ) * rad2deg;
    //     if (angle < 30) {
    //         // cout << "frontL" << endl;
    //         Eigen::Matrix3d rot = Eigen::AngleAxisd( (angle-30) * deg2rad, clav1.col(2) ).toRotationMatrix();
    //         C1.row(22) = C1.row(7) + ( rot * neck2head.transpose() ).transpose();
    //     } else if (angle > 150) {
    //         // cout << "frontR" << endl;
    //         Eigen::Matrix3d rot = Eigen::AngleAxisd( (angle-150) * deg2rad, clav1.col(2) ).toRotationMatrix();
    //         C1.row(22) = C1.row(7) + ( rot * neck2head.transpose() ).transpose();
    //     }
    // }

    // Get Head rotation
    bool isEyes = false;
    bool isEars = false;
    Eigen::RowVector3d eye_line, ear_line;
    
    if      ( C1(23,2) > -1000 && C1(24,2) > -1000 ) {
        isEyes = true;
        eye_line = (C1.row(23)-C1.row(24)).normalized();
    }  
    else if ( C1(25,2) > -1000 && C1(26,2) > -1000 )  {
        isEars = true;
        ear_line = (C1.row(25)-C1.row(26)).normalized();
    }

    if (isEyes && isEars) {
        head1.col(0) = (eye_line + ear_line).normalized();
    } else if (isEyes && !isEars) {
        head1.col(0) = eye_line;
    } else if (!isEyes && isEars) {
        head1.col(0) = ear_line;
    } else {
        head1.col(0) = clav1.col(0);
    }


    // Correction //////////////////////////
    // zPlane
    proj_zPlane = projectOntoPlane(head1.col(0), clav1.col(2));
    proj_zPlane.dot(clav1.col(1)) < 0 ? isFront = true : isFront = false;
    angle = acos( proj_zPlane.dot(clav1.col(0)) ) * rad2deg;
    if (angle > 45 && isFront) { // right side
        // cout << "right: " << angle << endl;
        Eigen::Matrix3d rot = Eigen::AngleAxisd( (angle-60) * deg2rad, clav1.col(2) ).toRotationMatrix();
        head1.col(0) = rot * head1.col(0);
    } 
    else if (angle > 45 && !isFront) { // left side
        // cout << "left: " << angle << endl;
        Eigen::Matrix3d rot = Eigen::AngleAxisd( -(angle-60) * deg2rad, clav1.col(2) ).toRotationMatrix();
        head1.col(0) = rot * head1.col(0);
    }
    // yPlane
    proj_yPlane = projectOntoPlane(head1.col(0), clav1.col(1));
    proj_yPlane.dot(clav1.col(2)) > 0 ? isUp = true : isUp = false;
    angle = acos( proj_yPlane.dot(clav1.col(0)) ) * rad2deg;
    if (angle > 30 && isUp) { // right side
        // cout << "down: " << angle << endl;
        Eigen::Matrix3d rot = Eigen::AngleAxisd( (angle-30) * deg2rad, clav1.col(1) ).toRotationMatrix();
        head1.col(0) = rot * head1.col(0);
    } 
    else if (angle > 30 && !isUp) { // left side
        // cout << "up: " << angle << endl;
        Eigen::Matrix3d rot = Eigen::AngleAxisd( -(angle-30) * deg2rad, clav1.col(1) ).toRotationMatrix();
        head1.col(0) = rot * head1.col(0);
    }

    // xPlane
    Eigen::RowVector3d neck2head = ( C1.row(22)-C1.row(7) ).normalized();
    if (neck2head.dot(clav1.col(2)) < 0) {
        // cout << "reverse" << endl;
        C1.row(22) = C1.row(7) + ( vQ[0].matrix() * bV.row(0).transpose() ).transpose();
    }
    else {
        proj_xPlane = projectOntoPlane(neck2head, clav1.col(0));
        proj_xPlane.dot(clav1.col(1)) < 0 ? isFront = true : isFront = false;
        angle = acos( proj_xPlane.dot(clav1.col(2)) ) * rad2deg;
        if (angle > 30+45 && isFront) {
            // cout << "front: " << angle << endl;
            Eigen::Matrix3d rot = Eigen::AngleAxisd( -(angle-45) * deg2rad, clav1.col(0) ).toRotationMatrix();
            C1.row(22) = C1.row(7) + ( rot * neck2head.transpose() * bL(0) ).transpose();
        }
        else if (!isFront) {
            // cout << "back!!!!!: " << angle << endl;
            C1.row(22) = C1.row(7) + ( vQ[0].matrix() * bV.row(0).transpose() ).transpose();
        }
    }
    // cout << C1.row(22) << endl;

    head1.col(2) = (C1.row(22) - C1.row(7)).normalized();
    head1.col(1) = head1.col(2).cross(head1.col(0)).normalized();
    if (head1.col(1).dot(clav1.col(1)) < 0) {
        // cout << "axis flip" << endl;
        head1.col(1) = -head1.col(1);
    }
        
    head1.col(2) = head1.col(0).cross(head1.col(1)).normalized();
    headQ = Eigen::Quaterniond(head1 * head0.transpose()); headQ.normalize();

    vQ[2] = chestQ;
    vQ[3] = clavQ;
    vQ[8] = clavQ;
    vQ[9] = clavQ;
    vQ[0] = headQ;
}

void PoseCorrector::LEGS_JointsCorrection()
{
    Eigen::RowVector3d hipR2hipL = C1.row(12) - C1.row(15);
    
    // @@ Knee joints ----------------------------------------
    // Matching //////////////////////////
    if ( C1(13,2) > -1000 && C1(16,2) < -1000) { // kneeL detected, kneeR not detected
        C1.row(16) = C1.row(15) + ( vQ[18].matrix() * (C1.row(13)-C1.row(12)).transpose() ).transpose();
    }
    else if ( C1(13,2) < -1000 && C1(16,2) > -1000 ) { // kneeL not detected, kneeR detected
        C1.row(13) = C1.row(12) + ( vQ[14].matrix() * (C1.row(16)-C1.row(15)).transpose() ).transpose();
    }
    else if ( C1(13,2) < -1000 && C1(16,2) < -1000 ) { // both not detected
        if ( C1_prev.rows() == 28 && C1_prev.row(13)(2) > -1000 && C1_prev.row(16)(2)) { // use previous frame data
            C1.row(13) = C1_prev.row(13) + root_trans;
            C1.row(16) = C1_prev.row(16) + root_trans;
        } else {
            C1.row(13) = C1.row(12) + ( vQ[14].matrix() * bV.row(15).transpose() ).transpose();
            C1.row(16) = C1.row(15) + ( vQ[18].matrix() * bV.row(19).transpose() ).transpose();
        }
    }

    // Correction //////////////////////////
    // Left //////////////////////////
    Eigen::RowVector3d hipL2kneeL = C1.row(13) - C1.row(12);
    if (hipL2kneeL.dot(root1.col(2)) > 0) { // Up
        C1.row(13) = C1.row(12) + ( vQ[14].matrix() * bV.row(15).transpose() ).transpose();
    }
    else if (hipL2kneeL.dot(root1.col(0)) < 0) {  // In
        C1.row(13) = C1.row(12) + projectOntoPlane(hipL2kneeL, root1.col(0));
    } 
    else { // Down
        proj_yPlane = projectOntoPlane(hipL2kneeL, root1.col(1)).normalized();
        angle = acos( proj_yPlane.dot(root1.col(0)) ) * rad2deg; // |->
        if (angle < (90-10))  {
            Eigen::Matrix3d rot = Eigen::AngleAxisd(-(angle-(90-10))*deg2rad, root1.col(1)).toRotationMatrix();
            C1.row(13) = C1.row(12) + (rot * hipL2kneeL.transpose()).transpose();
        }
    }

    hipL2kneeL = C1.row(13) - C1.row(12);
    proj_zPlane = projectOntoPlane(hipL2kneeL, root1.col(2)).normalized();
    proj_zPlane.dot(root1.col(1)) < 0 ? isFront = true : isFront = false;
    proj_xPlane = projectOntoPlane(hipL2kneeL, root1.col(0)).normalized();
    angle = acos( proj_xPlane.dot(root1.col(1)) ) * rad2deg;
    if ( !isFront && angle < (90-10) ) { // Back
        Eigen::Matrix3d rot = Eigen::AngleAxisd((angle-(90-10))*deg2rad, root1.col(0)).toRotationMatrix();
        C1.row(13) = C1.row(12) + (rot * hipL2kneeL.transpose()).transpose();
    }
    else if ( isFront && angle > (90+15) ) { // Front
        Eigen::Matrix3d rot = Eigen::AngleAxisd((angle-(90+15))*deg2rad, root1.col(0)).toRotationMatrix();
        C1.row(13) = C1.row(12) + (rot * hipL2kneeL.transpose()).transpose();
    }

    // Right //////////////////////////
    Eigen::RowVector3d hipR2kneeR = C1.row(16) - C1.row(15);
    if (hipR2kneeR.dot(root1.col(2)) > 0) { // Up
        C1.row(16) = C1.row(15) + ( vQ[18].matrix() * bV.row(19).transpose() ).transpose();
    }
    else if (hipR2kneeR.dot(root1.col(0)) > 0) { // In
        C1.row(16) = C1.row(15) + projectOntoPlane(hipR2kneeR, -root1.col(0));
    } 
    else {
        proj_yPlane = projectOntoPlane(hipR2kneeR, root1.col(1)).normalized();
        angle = acos( proj_yPlane.dot(root1.col(0)) ) * rad2deg; // |->
        if ( angle > (90+10) ) { // Out
            Eigen::Matrix3d rot = Eigen::AngleAxisd(-(angle-(90+10))*deg2rad, root1.col(1)).toRotationMatrix();
            C1.row(16) = C1.row(15) + (rot * hipR2kneeR.transpose()).transpose();
        }
    }
    
    hipR2kneeR = C1.row(16) - C1.row(15);
    proj_zPlane = projectOntoPlane(hipR2kneeR, root1.col(2)).normalized();
    proj_zPlane.dot(root1.col(1)) < 0 ? isFront = true : isFront = false;
    proj_xPlane = projectOntoPlane(hipR2kneeR, root1.col(0)).normalized();
    angle = acos( proj_xPlane.dot(root1.col(1)) ) * rad2deg;
    if (!isFront && angle < (90-10)) { // Back
        Eigen::Matrix3d rot = Eigen::AngleAxisd((angle-(90-10))*deg2rad, root1.col(0)).toRotationMatrix();
        C1.row(16) = C1.row(15) + (rot * hipR2kneeR.transpose()).transpose();
    }
    else if (isFront && angle > (90+15)) { // Front
        Eigen::Matrix3d rot = Eigen::AngleAxisd((angle-(90+15))*deg2rad, root1.col(0)).toRotationMatrix();
        C1.row(16) = C1.row(15) + (rot * hipR2kneeR.transpose()).transpose();
    }

    // if L/R bones have same direction, then correct them
    hipL2kneeL = C1.row(13) - C1.row(12);
    hipR2kneeR = C1.row(16) - C1.row(15);
    Eigen::RowVector3d L = projectOntoPlane( hipL2kneeL, root1.col(0) );
    Eigen::RowVector3d R = projectOntoPlane( hipR2kneeR, root1.col(0) );
    double angleL = acos(L.dot(-root1.col(2)) / L.norm()) * rad2deg;
    double angleR = acos(R.dot(-root1.col(2)) / R.norm()) * rad2deg;
    L = projectOntoPlane( hipL2kneeL, root1.col(2) );
    R = projectOntoPlane( hipR2kneeR, root1.col(2) );
    if (angleL > 10 && angleR > 10 && L.dot(R) > 0) {
        C1.row(13) = C1.row(12) + ( vQ[14].matrix() * bV.row(15).transpose() ).transpose();
        C1.row(16) = C1.row(15) + ( vQ[18].matrix() * bV.row(19).transpose() ).transpose();
    }
    // Rotation //////////////////////////
    set_rotation(C.row(13)-C.row(12), C1.row(13)-C1.row(12), vQ[14], vQ[15]);
    set_rotation(C.row(16)-C.row(15), C1.row(16)-C1.row(15), vQ[18], vQ[19]);
    
    // @@ Ankle joints correction ----------------------------------------
    // Matching //////////////////////////
    if ( C1(14,2) < -1000 ) {
        if ( C1_prev.rows() == 28 && C1_prev.row(14)(2) > -1000) { // use previous frame data
            C1.row(14)  = C1_prev.row(14)  + root_trans;
        } else {
            C1.row(14) = C1.row(13) + ( vQ[15].matrix() * bV.row(16).transpose() ).transpose();
        }
    }
    else if ( C1(17,2) < -1000) {
        if ( C1_prev.rows() == 28 && C1_prev.row(17)(2) > -1000) { // use previous frame data
            C1.row(17) = C1_prev.row(10) + root_trans;
        } else {
            C1.row(17) = C1.row(16) + ( vQ[19].matrix() * bV.row(20).transpose() ).transpose();
        }
    }
    
    // Correction //////////////////////////
    // Left //////////////////////////
    Eigen::RowVector3d kneeL2ankleL = C1.row(14) - C1.row(13);
    q.setFromTwoVectors( vQ[15].matrix().col(2), (C1.row(12) - C1.row(13)).normalized() );
    vQ_corr[15] = (q * vQ[15]).normalized();
    if (kneeL2ankleL.dot(vQ_corr[15].matrix().col(2)) > 0) { // prevent ankle from going upward
        C1.row(14) = C1.row(13) + ( vQ[15].matrix() * bV.row(16).transpose() ).transpose();
    }
    else if (kneeL2ankleL.dot(vQ_corr[15].matrix().col(0)) < 0) { // prevent ankle from crossing hip
        C1.row(14) = C1.row(13) + projectOntoPlane(kneeL2ankleL, vQ_corr[15].matrix().col(0));
    }
    kneeL2ankleL = C1.row(14) - C1.row(13);
    hipL2kneeL = C1.row(13) - C1.row(12);
    if ( kneeL2ankleL.normalized()(2) > hipL2kneeL.normalized()(2) ) { // prevent ankle from rising over knee
        C1.row(14) = C1.row(13) + hipL2kneeL.normalized() * bL(16);
    } else {
        kneeL2ankleL = C1.row(14) - C1.row(13);
        angle = acos(kneeL2ankleL.dot(hipL2kneeL) / (kneeL2ankleL.norm() * hipL2kneeL.norm())) * rad2deg;
        if (angle > 45) {
            Eigen::Matrix3d rot = Eigen::AngleAxisd(-(angle-45)*M_PI/180.0, vQ_corr[19].matrix().col(0)).toRotationMatrix();
            C1.row(14) = C1.row(13) + ( rot * kneeL2ankleL.transpose() ). transpose();
        }
    }
    
    // Right //////////////////////////
    Eigen::RowVector3d kneeR2ankleR = C1.row(17) - C1.row(16);
    q.setFromTwoVectors( vQ[19].matrix().col(2), (C1.row(15) - C1.row(16)).normalized() );
    vQ_corr[19] = (q * vQ[19]).normalized();
    if (kneeR2ankleR.dot(vQ_corr[19].matrix().col(2)) > 0) {
        C1.row(17) = C1.row(16) + ( vQ[19].matrix() * bV.row(20).transpose() ).transpose();
    }
    else if (kneeR2ankleR.dot(-vQ_corr[19].matrix().col(0)) < 0) {
        C1.row(17) = C1.row(16) + projectOntoPlane(kneeR2ankleR, -vQ_corr[15].matrix().col(0));
    }
    kneeR2ankleR = C1.row(17) - C1.row(16);
    hipR2kneeR = C1.row(16) - C1.row(15);
    if ( kneeR2ankleR.normalized()(2) > hipR2kneeR.normalized()(2) ) { // prevent ankle from rising over knee 
        C1.row(17) = C1.row(16) + hipR2kneeR.normalized() * bL(20);
    } else {
        kneeR2ankleR = (C1.row(17) - C1.row(16));
        angle = acos(kneeR2ankleR.dot(hipR2kneeR) / (kneeR2ankleR.norm() * hipR2kneeR.norm())) * rad2deg;
        if (angle > 45) { 
            Eigen::Matrix3d rot = Eigen::AngleAxisd(-(angle-45)*deg2rad, vQ_corr[19].matrix().col(0)).toRotationMatrix();
            C1.row(17) = C1.row(16) + ( rot * kneeR2ankleR.transpose() ). transpose();
        }
    }
    
    // Rotation //////////////////////////
    set_rotation(C.row(14)-C.row(13), C1.row(14)-C1.row(13), vQ[15], vQ[16]);
    set_rotation(C.row(17)-C.row(16), C1.row(17)-C1.row(16), vQ[19], vQ[20]);    
}

void PoseCorrector::ARMS_JointsCorrection()
{
    // @@ Shoulder joints ----------------------------------------
    // Matching //////////////////////////
    if ( C1(4,2) > -1000 && C1(9,2) < -1000) {
        C1.row(9) = C1.row(7) + ( C1.row(7) - C1.row(4) );
    }
    else if ( C1(4,2) < -1000 && C1(9,2) > -1000 ) {
        C1.row(4) = C1.row(7) + ( C1.row(7) - C1.row(9) );
    }
    else if ( C1(4,2) < -1000 && C1(9,2) < -1000) {
        if ( C1_prev.rows() == 28 && C1_prev.row(4)(2) > -1000 && C1_prev.row(9)(2) > -1000) { // use previous frame data
            C1.row(4) = C1_prev.row(4) + clav_trans;
            C1.row(9) = C1_prev.row(9) + clav_trans;
        }
        else {
            C1.row(4) = C1.row(3) + ( vQ[3] * bV.row(4).transpose() ).transpose();
            C1.row(9) = C1.row(8) + ( vQ[9] * bV.row(10).transpose() ).transpose();
        }
    }

    // Correction //////////////////////////
    // Left //////////////////////////
    // Front/Back (10/5)
    Eigen::RowVector3d clavL2shoulL = C1.row(4) - C1.row(3);
    proj_zPlane = projectOntoPlane(clavL2shoulL, clav1.col(2)).normalized();
    proj_zPlane.dot(clav1.col(1)) < 0 ? isFront = true : isFront = false;
    angle = acos( clav1.col(0).dot(proj_zPlane) ) * rad2deg;

    if (isFront && angle > 10) {
        Eigen::Matrix3d rot = Eigen::AngleAxisd((angle-10)*deg2rad, clav1.col(2)).toRotationMatrix();   
        C1.row(4) = C1.row(3) + (rot * clavL2shoulL.transpose()).transpose();
    }
    else if (!isFront && angle > 5) {
        Eigen::Matrix3d rot = Eigen::AngleAxisd(-(angle-5)*deg2rad, clav1.col(2)).toRotationMatrix();
        C1.row(4) = C1.row(3) + (rot * clavL2shoulL.transpose()).transpose();
    }
    
    // Up/Down (10/5)
    clavL2shoulL = C1.row(4) - C1.row(3);
    proj_yPlane = projectOntoPlane(clavL2shoulL, clav1.col(1)).normalized();
    proj_yPlane.dot(clav1.col(2)) > 0 ? isUp = true : isUp = false;
    angle = acos( clav1.col(0).dot(proj_yPlane) ) * rad2deg;
    
    if (isUp && angle > 10) {
        Eigen::Matrix3d rot = Eigen::AngleAxisd((angle-10)*deg2rad, clav1.col(1)).toRotationMatrix();    
        C1.row(4) = C1.row(3) + (rot * clavL2shoulL.transpose()).transpose();
    }
    else if (!isUp && angle > 30) {
        Eigen::Matrix3d rot = Eigen::AngleAxisd(-(angle-30)*deg2rad, clav1.col(1)).toRotationMatrix();
        C1.row(4) = C1.row(3) + (rot * clavL2shoulL.transpose()).transpose();
    }

    // Right //////////////////////////
    // Front/Back (10/5)
    Eigen::RowVector3d clavR2shoulR = C1.row(9) - C1.row(8);
    proj_zPlane = projectOntoPlane(clavR2shoulR, clav1.col(2)).normalized();
    proj_zPlane.dot(clav1.col(1)) < 0 ? isFront = true : isFront = false;
    angle = acos(-clav1.col(0).dot(proj_zPlane) ) * rad2deg;
    
    if (isFront && angle > 10) {
        Eigen::Matrix3d rot = Eigen::AngleAxisd(-(angle-10)*deg2rad, clav1.col(2)).toRotationMatrix();   
        C1.row(9) = C1.row(8) + (rot * clavR2shoulR.transpose()).transpose();
    }
    else if (!isFront && angle > 5) {
        Eigen::Matrix3d rot = Eigen::AngleAxisd((angle-5)*deg2rad, clav1.col(2)).toRotationMatrix();
        C1.row(9) = C1.row(8) + (rot * clavR2shoulR.transpose()).transpose();
    }
    
    // Up/Down (10/5)
    clavR2shoulR = C1.row(9) - C1.row(8);
    proj_yPlane = projectOntoPlane(clavR2shoulR, clav1.col(1)).normalized();
    proj_yPlane.dot(clav1.col(2)) > 0 ? isUp = true : isUp = false;
    angle = acos( -clav1.col(0).dot(proj_yPlane) ) * rad2deg;
    
    if (isUp && angle > 10) {
        Eigen::Matrix3d rot = Eigen::AngleAxisd(-(angle-10)*deg2rad, clav1.col(1)).toRotationMatrix();    
        C1.row(9) = C1.row(8) + (rot * clavR2shoulR.transpose()).transpose();
    }
    else if (!isUp && angle > 30) {
        Eigen::Matrix3d rot = Eigen::AngleAxisd((angle-30)*deg2rad, clav1.col(1)).toRotationMatrix();
        C1.row(9) = C1.row(8) + (rot * clavR2shoulR.transpose()).transpose();
    }
    
    // Rotation //////////////////////////
    set_rotation(C.row(4)-C.row(3), C1.row(4)-C1.row(3), vQ[3], vQ[4]);
    set_rotation(C.row(9)-C.row(8), C1.row(9)-C1.row(8), vQ[9], vQ[10]);


    // @@ Elbow joints ----------------------------------------
    // Matching //////////////////////////
    if ( C1(5,2) < -1000 ) {
        if ( C1_prev.rows() == 28 && C1_prev.row(5)(2) > -1000) { // use previous frame data
            C1.row(5)  = C1_prev.row(5)  + clav_trans;
        } else {
            C1.row(5) = C1.row(4) + ( vQ[4].matrix() * bV.row(5).transpose() ).transpose();
        }
    }
    if ( C1(10,2) < -1000) {
        if ( C1_prev.rows() == 28 && C1_prev.row(10)(2) > -1000) { // use previous frame data
            C1.row(10) = C1_prev.row(10) + clav_trans;
        } else {
            C1.row(10) = C1.row(9) + ( vQ[10].matrix() * bV.row(11).transpose() ).transpose();
        }
    }

    // Correction //////////////////////////
    bool isBackL(false), isBackR(false);
    // Left //////////////////////////
    Eigen::RowVector3d shoulL2elbowL = C1.row(5) - C1.row(4);
    Eigen::RowVector3d shoulR2shoulL = C1.row(4) - C1.row(9);
    if (shoulL2elbowL.dot(clav1.col(2)) > 0) {
        C1.row(5) = C1.row(4) + ( vQ[4].matrix() * bV.row(5).transpose() ).transpose();
    } else {
        proj_yPlane = projectOntoPlane(shoulL2elbowL, clav1.col(1)).normalized();
        angle = acos( shoulR2shoulL.dot(proj_yPlane) / (shoulR2shoulL.norm() )) * rad2deg;
        if (angle > 85) { // In
            Eigen::Matrix3d rot = Eigen::AngleAxisd(-(angle-80)*deg2rad, clav1.col(1)).toRotationMatrix();
            C1.row(5) = C1.row(4) + (rot * shoulL2elbowL.transpose()).transpose();
            
            // Front/Back
            shoulL2elbowL = C1.row(5) - C1.row(4);
            proj_xPlane = projectOntoPlane(shoulL2elbowL, clav1.col(0)).normalized();
            proj_xPlane.dot(clav1.col(1)) < 0 ? isFront = true : isFront = false;
            if(!isFront) isBackL = true;
            angle = acos( clav1.col(2).dot(proj_xPlane) ) * rad2deg;

            if (isFront && angle > (180-40)) {
                // cout << "isFront - L" << endl;
                rot = Eigen::AngleAxisd(-(angle-(180-40))*deg2rad, shoulR2shoulL.normalized()).toRotationMatrix();
                C1.row(5) = C1.row(4) + (rot * shoulL2elbowL.transpose()).transpose();
            }
            else if (!isFront && angle < 180-40) {
                // cout << "isBack - L" << endl;
                rot = Eigen::AngleAxisd((angle-(180-40))*deg2rad, shoulR2shoulL.normalized()).toRotationMatrix();
                C1.row(5) = C1.row(4) + (rot * shoulL2elbowL.transpose()).transpose();
            }
        }
    }


    // // Right //////////////////////////
    Eigen::RowVector3d shoulR2elbowR = C1.row(10) - C1.row(9);
    Eigen::RowVector3d shoulL2shoulR = C1.row(9) - C1.row(4);
    if (shoulR2elbowR.dot(clav1.col(2)) > 0) {
        C1.row(10) = C1.row(9) + ( vQ[10].matrix() * bV.row(11).transpose() ).transpose();
    } else {
        proj_yPlane = projectOntoPlane(shoulR2elbowR, clav1.col(1)).normalized();
        angle = acos(shoulL2shoulR.dot(proj_yPlane) / (shoulL2shoulR.norm() )) * rad2deg;
        if (angle > 85) { // In
            Eigen::Matrix3d rot = Eigen::AngleAxisd((angle-80)*deg2rad, clav1.col(1)).toRotationMatrix();
            C1.row(10) = C1.row(9) + (rot * shoulR2elbowR.transpose()).transpose();
            
            // Front/Back
            shoulR2elbowR = C1.row(10) - C1.row(9);
            proj_xPlane = projectOntoPlane(shoulR2elbowR, clav1.col(0)).normalized();
            proj_xPlane.dot(clav1.col(1)) < 0 ? isFront = true : isFront = false;
            if(!isFront) isBackR = true;
            angle = acos( clav1.col(2).dot(proj_xPlane) ) * rad2deg;
            
            if (isFront && angle > 180-40) {
                rot = Eigen::AngleAxisd(-(angle-(180-40))*deg2rad, shoulR2shoulL.normalized()).toRotationMatrix();
                C1.row(10) = C1.row(9) + (rot * shoulR2elbowR.transpose()).transpose();
            }
            else if (!isFront && angle < 180-40) {
                rot = Eigen::AngleAxisd((angle-(180-40))*deg2rad, shoulR2shoulL.normalized()).toRotationMatrix();
                C1.row(10) = C1.row(9) + (rot * shoulR2elbowR.transpose()).transpose();
            }
        }
    }
    
    // Rotation //////////////////////////
    set_rotation(C.row(5)-C.row(4), C1.row(5)-C1.row(4), vQ[4], vQ[5]);
    set_rotation(C.row(10)-C.row(9), C1.row(10)-C1.row(9), vQ[10], vQ[11]);

    // @@ Wrist joints ----------------------------------------
    // Matching //////////////////////////
    if ( C1(6,2) < -1000 ) {
        if ( C1_prev.rows() == 28 && C1_prev.row(6)(2) > -1000) { // use previous frame data
            C1.row(6)  = C1_prev.row(6)  + clav_trans;
        } else {
            C1.row(6) = C1.row(5) + ( vQ[5].matrix() * bV.row(6).transpose() ).transpose();
        }
    }
    else if ( C1(11,2) < -1000) {
        if ( C1_prev.rows() == 28 && C1_prev.row(11)(2) > -1000) { // use previous frame data
            C1.row(11) = C1_prev.row(11) + clav_trans;
        } else {
            C1.row(11) = C1.row(10) + ( vQ[11].matrix() * bV.row(12).transpose() ).transpose();
        }
    }

    // Correction //////////////////////////
    double zL_in_limit(45), zR_in_limit(45);
    if (isBackL) zL_in_limit = 90;
    if (isBackR) zR_in_limit = 90;
    // Left //////////////////////////
    Eigen::RowVector3d elbowL2wristL = C1.row(6) - C1.row(5);
    Eigen::RowVector3d elbowL2shoulL = C1.row(4) - C1.row(5);
    q.setFromTwoVectors(vQ[5].matrix().col(2), elbowL2shoulL.normalized());
    vQ_corr[5] = (q * vQ[5]).normalized();

    // zPlane limit
    proj_zPlane = projectOntoPlane(elbowL2wristL, clav1.col(2)).normalized();
    proj_zPlane.dot(clav1.col(1)) < 0 ? isFront = true : isFront = false;
    angle = acos( proj_zPlane.dot(clav1.col(0)) ) * rad2deg;
    if (isFront && angle < zL_in_limit) {
        // cout << "zfrontL-L: " << angle << endl;
        Eigen::Matrix3d rot = Eigen::AngleAxisd((angle-zL_in_limit)*deg2rad, clav1.col(2)).toRotationMatrix();
        C1.row(6) = C1.row(5) + (rot * elbowL2wristL.transpose()).transpose();
    } else if (isFront && angle > 135) {
        // cout << "zfrontL-R: " << angle << endl;
        Eigen::Matrix3d rot = Eigen::AngleAxisd((angle-135)*deg2rad, clav1.col(2)).toRotationMatrix();
        C1.row(6) = C1.row(5) + (rot * elbowL2wristL.transpose()).transpose();
    } else if (!isFront) {
        // cout << "zbackL: " << angle << endl;
        C1.row(6) = C1.row(5) + ( vQ[5].matrix() * bV.row(6).transpose() ).transpose();
    }
    proj_zPlane = projectOntoPlane(C1.row(5)-C1.row(6), vQ_corr[5].matrix().col(2)).normalized();
    angle = acos( proj_zPlane.dot(vQ_corr[5].matrix().col(1)) );
    vQ_corr[5] = Eigen::AngleAxisd(-angle, vQ_corr[5].matrix().col(2)) * vQ_corr[5];

    // xPlane limit
    elbowL2wristL = C1.row(6) - C1.row(5);
    proj_xPlane = projectOntoPlane(elbowL2wristL, vQ_corr[5].matrix().col(0)).normalized();
    proj_xPlane.dot(vQ_corr[5].matrix().col(1)) < 0 ? isFront = true : isFront = false;
    angle = acos(elbowL2shoulL.dot(proj_xPlane) / (elbowL2shoulL.norm())) * rad2deg;
    if (isFront && angle < 45) {
        // cout << "xfrontL: " << angle << endl;
        Eigen::Matrix3d rot = Eigen::AngleAxisd(-(angle-45)*deg2rad, vQ_corr[5].matrix().col(0)).toRotationMatrix();
        C1.row(6) = C1.row(5) + (rot * elbowL2wristL.transpose()).transpose();
    }
    // else if (!isFront) {
    //     cout << "xbackL: " << angle << endl;
    //     C1.row(6) = C1.row(5) -elbowL2shoulL.normalized() * bL(6);
    // }
    proj_zPlane = projectOntoPlane(C1.row(5)-C1.row(6), vQ_corr[5].matrix().col(2)).normalized();
    angle = acos( proj_zPlane.dot(vQ_corr[5].matrix().col(1)) );
    vQ_corr[5] = Eigen::AngleAxisd(-angle, vQ_corr[5].matrix().col(2)) * vQ_corr[5];
    
    
    // Right //////////////////////////
    Eigen::RowVector3d elbowR2wristR = C1.row(11) - C1.row(10);
    Eigen::RowVector3d elbowR2shoulR = C1.row(9)  - C1.row(10);
    q.setFromTwoVectors(vQ[11].matrix().col(2), elbowR2shoulR.normalized());
    vQ_corr[11] = (q * vQ[11]).normalized();

    // zPlane limit
    proj_zPlane = projectOntoPlane(elbowR2wristR, clav1.col(2)).normalized();
    proj_zPlane.dot(clav1.col(1)) < 0 ? isFront = true : isFront = false;
    angle = acos( proj_zPlane.dot(clav1.col(0)) ) * rad2deg;
    if (isFront && angle < zR_in_limit) {
        // cout << "zfrontR-L: " << angle << endl;
        Eigen::Matrix3d rot = Eigen::AngleAxisd((angle-zR_in_limit)*deg2rad, clav1.col(2)).toRotationMatrix();
        C1.row(11) = C1.row(10) + (rot * elbowR2wristR.transpose()).transpose();
    } else if (isFront && angle > 135) {
        // cout << "zfrontR-R: " << angle << endl;
        Eigen::Matrix3d rot = Eigen::AngleAxisd((angle-135)*deg2rad, clav1.col(2)).toRotationMatrix();
        C1.row(11) = C1.row(10) + (rot * elbowR2wristR.transpose()).transpose();
    } else if (!isFront) {
        // cout << "zbackR: " << angle << endl;
        C1.row(11) = C1.row(10) + ( vQ[11].matrix() * bV.row(12).transpose() ).transpose();
    }
    proj_zPlane = projectOntoPlane(C1.row(10)-C1.row(11), vQ_corr[11].matrix().col(2)).normalized();
    angle = acos( proj_zPlane.dot(vQ_corr[11].matrix().col(1)) );
    vQ_corr[11] = Eigen::AngleAxisd(angle, vQ_corr[11].matrix().col(2)) * vQ_corr[11];
        
    
    // xPlane limit
    elbowR2wristR = C1.row(11) - C1.row(10);
    proj_xPlane = projectOntoPlane(elbowR2wristR, vQ_corr[11].matrix().col(0)).normalized();
    proj_xPlane.dot(vQ_corr[11].matrix().col(1)) < 0 ? isFront = true : isFront = false;
    angle = acos(elbowR2shoulR.dot(proj_xPlane) / (elbowR2shoulR.norm())) * rad2deg;
    if (isFront && angle < 45) {
        // cout << "xfrontR: " << angle << endl;
        Eigen::Matrix3d rot = Eigen::AngleAxisd(-(angle-45)*deg2rad, vQ_corr[11].matrix().col(0)).toRotationMatrix();
        C1.row(11) = C1.row(10) + (rot * elbowR2wristR.transpose()).transpose();
    }
    // else if (!isFront) {
    //     cout << "xbackR: " << angle << endl;
    //     C1.row(11) = C1.row(10) -elbowR2shoulR.normalized() * bL(12);
    // }
    proj_zPlane = projectOntoPlane(C1.row(10)-C1.row(11), vQ_corr[11].matrix().col(2)).normalized();
    angle = acos( proj_zPlane.dot(vQ_corr[11].matrix().col(1)) );
    vQ_corr[11] = Eigen::AngleAxisd(angle, vQ_corr[11].matrix().col(2)) * vQ_corr[11];

    // Rotation //////////////////////////
    set_rotation(C.row(6)-C.row(5), C1.row(6)-C1.row(5), vQ[5], vQ[6]);
    set_rotation(C.row(11)-C.row(10), C1.row(11)-C1.row(10), vQ[11], vQ[12]);



}

void PoseCorrector::LIMBS_JointsCorrection()
{   
    // Foot L,R
    C1.row(20) = C1.row(14) + ( vQ[16].matrix() * bV.row(17).transpose() ).transpose();
    C1.row(21) = C1.row(17) + ( vQ[20].matrix() * bV.row(21).transpose() ).transpose();
    vQ[17] = vQ[16];
    vQ[21] = vQ[20];


    // @@ Hand joints ----------------------------------------
    // Matching //////////////////////////
    // if ( C1(18,2) < -1000 ) {
    //     if ( C1_prev.rows() == 28 && C1_prev(18,2) > -1000) { // use previous frame data
    //         C1.row(18)  = C1_prev.row(18)  + clav_trans;
    //     } else {
    //         C1.row(18) = C1.row(6) + ( vQ[6].matrix() * bV.row(7).transpose() ).transpose();
    //     }
    // }
    // else if ( C1(19,2) < -1000) {
    //     if ( C1_prev.rows() == 28 && C1_prev(11,2) > -1000) { // use previous frame data
    //         C1.row(19) = C1_prev.row(19) + clav_trans;
    //     } else {
    //         C1.row(19) = C1.row(11) + ( vQ[12].matrix() * bV.row(13).transpose() ).transpose();
    //     }
    // }

    C1.row(18) = C1.row(6) + ( vQ[6].matrix() * bV.row(7).transpose() ).transpose();
    C1.row(19) = C1.row(11) + ( vQ[12].matrix() * bV.row(13).transpose() ).transpose();    
    vQ[7] = vQ[6];
    vQ[13] = vQ[12];

    // if ( (C1.row(18)-isocenter).squaredNorm() < 10000 && C1(18,2) > isocenter(2)
    //   || (C1.row(19)-isocenter).squaredNorm() < 10000 && C1(19,2) > isocenter(2) ) {
    //     Eigen::Matrix3d handL;
    //     handL.col(0) = clav1.col(2);
    //     handL.col(1) = clav1.col(0);
    //     handL.col(2) = clav1.col(1);
    //     q.setFromTwoVectors( handL.col(2), (C1.row(6) - C1.row(18)).normalized() );
    //     Eigen::Matrix3d rot = Eigen::AngleAxisd(30*deg2rad, handL.col(2)).toRotationMatrix();
    //     handL = rot * q * handL;
    //     Eigen::Matrix3d handR;
    //     handR.col(0) = -clav1.col(2);
    //     handR.col(1) = -clav1.col(0);
    //     handR.col(2) = clav1.col(1);
    //     q.setFromTwoVectors( handR.col(2), (C1.row(11) - C1.row(19)).normalized() );
    //     handR = q * handR;
    //     vQ[7] = handL;
    //     vQ[13] = handR;
    // }
}

void PoseCorrector::HEAD_xyCorrection()
{   
    C1(27,2) = C1(7,2);
    Eigen::RowVector3d xyTrans = C1.row(27) - C1.row(7);
    C1 = C1.rowwise() + xyTrans;
}

void PoseCorrector::use_default_Larm()
{
    // cout << id << " use default Larm" << endl;
    C1.row(4) = C1.row(3) + ( vQ[3] * bV.row(4).transpose() ).transpose();
    set_rotation(C.row(4)-C.row(3), C1.row(4)-C1.row(3), vQ[3], vQ[4]);
    C1.row(5) = C1.row(4) + ( vQ[4] * bV.row(5).transpose() ).transpose();
    set_rotation(C.row(5)-C.row(4), C1.row(5)-C1.row(4), vQ[4], vQ[5]);
    C1.row(6) = C1.row(5) + ( vQ[5] * bV.row(6).transpose() ).transpose();
    set_rotation(C.row(6)-C.row(5), C1.row(6)-C1.row(5), vQ[5], vQ[6]);
    C1.row(18) = C1.row(6) + ( vQ[3] * bV.row(4).transpose() ).transpose();
    vQ[7] = vQ[6];
}

void PoseCorrector::use_default_Rarm()
{
    // cout << id << " use default Rarm" << endl;
    C1.row(9) = C1.row(8) + ( vQ[9] * bV.row(10).transpose() ).transpose();
    set_rotation(C.row(9)-C.row(8), C1.row(9)-C1.row(8), vQ[9], vQ[10]);
    C1.row(10) = C1.row(9) + ( vQ[10] * bV.row(11).transpose() ).transpose();
    set_rotation(C.row(10)-C.row(9), C1.row(10)-C1.row(9), vQ[10], vQ[11]);
    C1.row(11) = C1.row(10) + ( vQ[11] * bV.row(12).transpose() ).transpose();
    set_rotation(C.row(11)-C.row(10), C1.row(11)-C1.row(10), vQ[11], vQ[12]);
    C1.row(19) = C1.row(11) + ( vQ[9] * bV.row(10).transpose() ).transpose();
    vQ[13] = vQ[12];
}