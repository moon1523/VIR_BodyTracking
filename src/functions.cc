#include "functions.hh"

MatrixXd GenerateBonePoints(MatrixXd C, MatrixXi BE, double interval){
    MatrixXd output=C;
    for(int i=0;i<BE.rows();i++){
        Vector3d c0 = C.row(BE(i,0)).transpose();
        Vector3d c1 = C.row(BE(i,1)).transpose();
        Vector3d dir = (c1-c0).normalized();
        double l = (c1-c0).norm();
        int num = floor(l/interval);
        double interval1 = l/(double)num;
        int prevSize = output.rows();
        output.conservativeResize(prevSize+num-1,NoChange);
        for(double intvl=1;intvl<num;intvl++)
            output.row(prevSize+intvl-1) = (c0+dir*interval1*intvl).transpose();
    }
    return output;
}
map<tuple<int,int,int>,vector<int>> GenerateGrid(MatrixXd V, double size){
    map<tuple<int,int,int>,vector<int>> grid;
    for(int i=0;i<V.rows();i++){
        int x = floor(V(i,0)*size+0.5);
        int y = floor(V(i,1)*size+0.5);
        int z = floor(V(i,2)*size+0.5);
        auto key = make_tuple(x,y,z);
        if(grid.find(key)==grid.end()) grid[key]={};
        grid[key].push_back(i);
    }
    return grid;
}

map<int, map<int, double>> GenerateBarycentricCoord(MatrixXd V_f, MatrixXi T_f, MatrixXd V){
    map<tuple<int,int,int>,vector<int>> grid = GenerateGrid(V);
    map<int, map<int, double>> baryCoords;
    double epsl(-1e-5);
    for(int n=0;n<T_f.rows();n++){
        vector<Vector3d> tet;
        Vector3d max(-1e10,-1e10,-1e10), min(1e10,1e10,1e10);
        for(int e=0;e<4;e++){
            tet.push_back(V_f.row(T_f(n,e)).transpose());
            max(0) = max(0)>tet[e](0)? max(0):tet[e](0);
            max(1) = max(1)>tet[e](1)? max(1):tet[e](1);
            max(2) = max(2)>tet[e](2)? max(2):tet[e](2);
            min(0) = min(0)<tet[e](0)? min(0):tet[e](0);
            min(1) = min(1)<tet[e](1)? min(1):tet[e](1);
            min(2) = min(2)<tet[e](2)? min(2):tet[e](2);
        }
        double invVol6 = 1./(tet[1]-tet[0]).cross(tet[2]-tet[0]).dot(tet[3]-tet[0]);
        int i_max = floor(max(0)+0.5);int i_min = floor(min(0)+0.5);
        int j_max = floor(max(1)+0.5);int j_min = floor(min(1)+0.5);
        int k_max = floor(max(2)+0.5);int k_min = floor(min(2)+0.5);
        for(int i=i_min;i<i_max+1;i++){
            for(int j=j_min;j<j_max+1;j++){
                for(int k=k_min;k<k_max+1;k++){
                    auto key = make_tuple(i,j,k);
                    for(int idx:grid[key]){
                        if(baryCoords.find(idx)!=baryCoords.end()) continue;
                        Vector3d v = V.row(idx).transpose();
                        double b0 = (tet[1]-v).cross(tet[2]-v).dot(tet[3]-v)*invVol6;
                        if(b0<epsl) continue;
                        double b1 = (v-tet[0]).cross(tet[2]-tet[0]).dot(tet[3]-tet[0])*invVol6;
                        if(b1<epsl) continue;
                        double b2 = (tet[1]-tet[0]).cross(v-tet[0]).dot(tet[3]-tet[0])*invVol6;
                        if(b2<epsl) continue;
                        double b3 = (tet[1]-tet[0]).cross(tet[2]-tet[0]).dot(v-tet[0])*invVol6;
                        if(b3<epsl) continue;
                        map<int, double> bary;
                        bary[T_f(n,0)] = b0; bary[T_f(n,1)] = b1; bary[T_f(n,2)] = b2; bary[T_f(n,3)] = b3;
                        baryCoords[idx] = bary;
                    }
                }
            }
        }
        cout<<"\rGenerating barycentric coord..."<<baryCoords.size()<<"/"<<V.rows()<<"       "<<flush;
        //if((int)baryCoords.size()==V.rows()) break;
    }cout<<endl;
    if((int)baryCoords.size()!=V.rows()){
        cout<<"Check if all the vertices are in frame model!!"<<endl; exit(100);
    }
    return baryCoords;
}
void PrintBaryCoords(string file, map<int, map<int, double>> &baryCoords){
    ofstream ofs(file); ofs<<baryCoords.size()<<endl;
    for(size_t i=0;i<baryCoords.size();i++){
        ofs<<baryCoords[i].size()<<" ";
        for(auto iter:baryCoords[i])
            ofs<<iter.first<<" "<<iter.second<<" ";
        ofs<<endl;
    }ofs.close();
}
map<int, map<int, double>> ReadBaryFile(string file){
    map<int, map<int, double>> baryCoords;
    ifstream ifs(file);
    if(!ifs.is_open()) {
        cout<<file+" is not open!!!"<<endl;
        return baryCoords;
    }
    int num; ifs>>num;
    for(int i=0;i<num;i++){
        int count; ifs>>count;
        int vID; double w;
        map<int, double> bary;
        for(int n=0;n<count;n++){
            ifs>>vID>>w;
            bary[vID]=w;
        }
        baryCoords[i] = bary;
    }ifs.close();
    return baryCoords;
}
SparseMatrix<double> GenerateBarySparse(map<int, map<int, double>> &baryCoords, int v_size){
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;
    for(size_t i=0;i<baryCoords.size();i++){
        for(auto w:baryCoords[i]){
            triplets.push_back(T(i,w.first,w.second));
        }
    }
    SparseMatrix<double> mat(baryCoords.size(),v_size);
    mat.setFromTriplets(triplets.begin(),triplets.end());
    return mat;
}
bool CalculateScalingWeights(MatrixXd &C, MatrixXd &V,MatrixXi &T,MatrixXd &W){
    VectorXi b; b.resize(C.rows());
    MatrixXd bc = MatrixXd::Zero(C.rows(),C.rows());
    map<tuple<int,int,int>,vector<int>> grid = GenerateGrid(V);
    for(int n=0;n<C.rows();n++){
        int x = floor(C(n,0)+0.5);
        int y = floor(C(n,1)+0.5);
        int z = floor(C(n,2)+0.5);
        auto key = make_tuple(x,y,z);
        for(int i:grid[key]){
            if(fabs(C(n,0)-V(i,0))>0.01) continue;
            if(fabs(C(n,1)-V(i,1))>0.01) continue;
            if(fabs(C(n,2)-V(i,2))>0.01) continue;
            b(n) = i;
            bc(n,n) = 1;
            break;
        }
    }

    // compute BBW weights matrix
    igl::BBWData bbw_data;
    // only a few iterations for sake of demo
    bbw_data.active_set_params.max_iter = 10;
    bbw_data.verbosity = 2;
    return igl::bbw(V,T,b,bc,bbw_data,W);
}
bool CalculateScalingWeights(MatrixXd &C, MatrixXd &V,MatrixXi &T,MatrixXd &W, vector<int> eyeLV, vector<int> eyeRV){
    VectorXi b; b.resize(C.rows()-3+eyeLV.size()+eyeRV.size());
    MatrixXd bc = MatrixXd::Zero(C.rows()-3+eyeLV.size()+eyeRV.size(),C.rows()-1);
    map<tuple<int,int,int>,vector<int>> grid = GenerateGrid(V);
    int n=0;
    for(;n<C.rows()-3;n++) // exclude head joints
    {
        int x = floor(C(n,0)+0.5);
        int y = floor(C(n,1)+0.5);
        int z = floor(C(n,2)+0.5);
        auto key = make_tuple(x,y,z);
        for(int i:grid[key]){
            if(fabs(C(n,0)-V(i,0))>0.01) continue;
            if(fabs(C(n,1)-V(i,1))>0.01) continue;
            if(fabs(C(n,2)-V(i,2))>0.01) continue;
            b(n) = i;
            bc(n,n) = 1;
            break;
        }
    }
    for(int i:eyeLV){b(n) = i; bc(n++, C.rows()-3)=1; }
    for(int i:eyeRV){b(n) = i; bc(n++, C.rows()-2)=1; }
    // compute BBW weights matrix
    igl::BBWData bbw_data;
    // only a few iterations for sake of demo
    bbw_data.active_set_params.max_iter = 10;
    bbw_data.verbosity = 2;
    return igl::bbw(V,T,b,bc,bbw_data,W);
}
AngleAxisd GetRotMatrix(Vector3d from, Vector3d to) {
    Vector3d axis = from.cross(to); axis.normalize();
    double theta = acos(from.dot(to) / (from.norm()*to.norm()));
    return AngleAxisd(theta, axis);
}

void myDqs(
  const Eigen::MatrixXd & V,
  const std::vector<std::map<int, double>> & W,
  const RotationList & vQ,
  const std::vector<Vector3d> & vT,
  Eigen::MatrixXd & U)
{
  using namespace std;
//  assert(V.rows() <= W.rows());
//  assert(W.cols() == (int)vQ.size());
//  assert(W.cols() == (int)vT.size());
  // resize output
  U.resizeLike(V);

  // Convert quats + trans into dual parts
  vector<Eigen::Quaterniond> vD(vQ.size());
  for(size_t c = 0;c<vQ.size();c++)
  {
    const Eigen::Quaterniond  & q = vQ[c];
    vD[c].w() = -0.5*( vT[c](0)*q.x() + vT[c](1)*q.y() + vT[c](2)*q.z());
    vD[c].x() =  0.5*( vT[c](0)*q.w() + vT[c](1)*q.z() - vT[c](2)*q.y());
    vD[c].y() =  0.5*(-vT[c](0)*q.z() + vT[c](1)*q.w() + vT[c](2)*q.x());
    vD[c].z() =  0.5*( vT[c](0)*q.y() - vT[c](1)*q.x() + vT[c](2)*q.w());
  }

  // Loop over vertices
  const int nv = V.rows();
#pragma omp parallel for if (nv>10000)
  for(int i = 0;i<nv;i++)
  {
    Eigen::Quaterniond b0(0,0,0,0);
    Eigen::Quaterniond be(0,0,0,0);
    Eigen::Quaterniond vQ0;
    bool first(true);
    // Loop over handles
    for(auto iter:W[i])
    {
        if(first){
            b0.coeffs() = iter.second * vQ[iter.first].coeffs();
            be.coeffs() = iter.second * vD[iter.first].coeffs();
            vQ0 = vQ[iter.first];
            first = false;
            continue;
        }
        if( vQ0.dot( vQ[iter.first] ) < 0.f ){
            b0.coeffs() -= iter.second * vQ[iter.first].coeffs();
            be.coeffs() -= iter.second * vD[iter.first].coeffs();
        }else{
            b0.coeffs() += iter.second * vQ[iter.first].coeffs();
            be.coeffs() += iter.second * vD[iter.first].coeffs();
        }
    }
    Eigen::Quaterniond ce = be;
    ce.coeffs() /= b0.norm();
    Eigen::Quaterniond c0 = b0;
    c0.coeffs() /= b0.norm();
    // See algorithm 1 in "Geometric skinning with approximate dual quaternion
    // blending" by Kavan et al
    Vector3d v = V.row(i);
    Vector3d d0 = c0.vec();
    Vector3d de = ce.vec();
    Eigen::Quaterniond::Scalar a0 = c0.w();
    Eigen::Quaterniond::Scalar ae = ce.w();
    U.row(i) =  v + 2*d0.cross(d0.cross(v) + a0*v) + 2*(a0*de - ae*d0 + d0.cross(de));
  }

}
