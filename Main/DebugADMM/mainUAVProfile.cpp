#include <ADMM/Deformable.h>
#include <ADMM/DeformableVisualizer.h>
#include <Environment/EnvironmentVisualizer.h>
#include <Utils/VTKWriter.h>
#include <Utils/Utils.h>

using namespace PHYSICSMOTION;
using namespace DRAWER;
typedef FLOAT T;
DECL_MAT_VEC_MAP_TYPES_T

void addTwoCircle(std::vector<Vec3T>& waypoints,T d,T d2) {
  Vec3T pos=waypoints.back();
  //first circle
  waypoints.push_back(pos+Vec3T(d,0,0));
  waypoints.push_back(pos+Vec3T(d,-d2/2,0));
  waypoints.push_back(pos+Vec3T(0,-d2/2,0));
  waypoints.push_back(pos+Vec3T(0,0,0));
  //second circle
  waypoints.push_back(pos+Vec3T(d,0,0));
  waypoints.push_back(pos+Vec3T(d, d2/2,0));
  waypoints.push_back(pos+Vec3T(0, d2/2,0));
  waypoints.push_back(pos+Vec3T(0,0,0));
  //advance
  waypoints.push_back(pos+Vec3T(d,0,0));
}
void addCircle(Vec3T ctr,T r,T r2,MeshExact& m,int RES=16) {
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  for(int i=0; i<RES; i++) {
    T angle=(M_PI*2*i)/RES+M_PI/4;
    vss.push_back(Eigen::Matrix<double,3,1>(0,cos(angle),sin(angle))*r+ctr);
  }
  for(int i=0; i<RES; i++) {
    T angle=(M_PI*2*i)/RES+M_PI/4;
    vss.push_back(Eigen::Matrix<double,3,1>(0,cos(angle),sin(angle))*r2+ctr);
  }
  for(int i=0; i<RES; i++) {
    int v0=i;
    int v1=(i+1)%RES;
    int v2=((i+1)%RES)+RES;
    int v3=i+RES;
    iss.push_back({v0,v1,v2});
    iss.push_back({v0,v2,v3});
  }
  m=MeshExact(vss,iss);
}
void solve(int N) {
  T d=2.5,d2=5,r=1.5,r2=1.0;
  MeshExact m,m2,obs;
  std::vector<Vec3T> waypoints({Vec3T::Zero()});
  for(int i=0; i<=N; i++) {
    T alpha=(float)abs(i-N/2)/(float)(N/2);
    T coef=.15*(1-alpha)+1*alpha;
    T coef2=3*(1-alpha)+1*alpha;
    addCircle(Vec3T(d*(i+.5),-d2*coef2/2,0),r*coef,r2*coef,m);
    addCircle(Vec3T(d*(i+.5), d2*coef2/2,0),r*coef,r2*coef,m2);
    addTwoCircle(waypoints,d,d2*coef2);
    obs.addMesh(m);
    obs.addMesh(m2);
  }
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  obs.getMesh(vss,iss);
  obs=MeshExact(vss,iss);

  //solver
  Deformable<3> solver;
  solver.setObstacle(obs);
  solver.setDt(0);
  solver.setR(0);
  solver.setG(Vec3T(0,0,0));
  ASSERT_MSG(solver.setUAVTrajectory(waypoints,(int)waypoints.size()-1,(int)waypoints.size()*5),"Failed initialization!")
  solver.fix(0,0);
  solver.fix(solver.x().size()/3-1,0);
  std::cout << "Using " << solver.getUAVTraj().getTraj().getNumSegment() << " segments!" << std::endl;

  OptimizerParam param;
  param._initBeta=1e2f;
  param._initBetaX=0;
  param._tolG=1e-2f;
  param._maxIter=1e5;
  //param._debugGradientI=50;
  param._type=OptimizerParam::ADMM;
  solver.solve(param);
}
int main(int argc, char** argv) {
  for (int n : {10,20,30,40,50}) {
    std::string name = "UAVProfile" + std::to_string(n) + ".txt";
    if (exists(name))
      continue;
    if (freopen(name.c_str(), "w", stdout) == NULL) {
        // Handle error if redirection fails
        perror("Failed to redirect stdout");
        return 1;
    }
    std::cout << "Begin-testing " << n << " circles" << std::endl;
    solve(n);
    std::cout << "End-testing " << n << " circles" << std::endl;
  }
  return 0;
}
