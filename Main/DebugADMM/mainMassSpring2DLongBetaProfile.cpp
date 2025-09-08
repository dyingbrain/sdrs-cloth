#include <ADMM/Deformable.h>
#include <ADMM/DeformableVisualizer.h>
#include <Environment/EnvironmentVisualizer.h>
#include <Utils/VTKWriter.h>
#include <Utils/Utils.h>

using namespace PHYSICSMOTION;
using namespace DRAWER;
typedef FLOAT T;
DECL_MAT_VEC_MAP_TYPES_T

void createGrid(int N,Eigen::Matrix<double,2,1> a,Eigen::Matrix<double,2,1> b,MeshExact& m) {
  std::vector<Eigen::Matrix<double,2,1>> vss;
  std::vector<Eigen::Matrix<int,2,1>> iss;
  for(int i=0; i<=N; i++)
    vss.push_back((b-a)*i/N+a);
  for(int i=0; i<N; i++)
    iss.push_back(Eigen::Matrix<int,2,1>(i,i+1));
  m=MeshExact(vss,iss,true);
}
void createBox(Eigen::Matrix<double,2,1> ctr,Eigen::Matrix<double,2,1> sz,MeshExact& m) {
  std::vector<Eigen::Matrix<double,2,1>> vss;
  std::vector<Eigen::Matrix<int,2,1>> iss;
  vss.push_back({ctr[0]-sz[0],ctr[1]-sz[1]});
  vss.push_back({ctr[0]+sz[0],ctr[1]-sz[1]});
  vss.push_back({ctr[0]+sz[0],ctr[1]+sz[1]});
  vss.push_back({ctr[0]-sz[0],ctr[1]+sz[1]});
  iss.push_back({0,1});
  iss.push_back({1,2});
  iss.push_back({2,3});
  iss.push_back({3,0});
  m=MeshExact(vss,iss,true);
}
void solve(FLOAT initBeta) {
  //build mesh grid
  int N=100,D=100;
  MeshExact m,m2;
  createGrid(N,Eigen::Matrix<double,2,1>(0,0),Eigen::Matrix<double,2,1>(1,D),m);
  createBox(Eigen::Matrix<double,2,1>(0,N/2),Eigen::Matrix<double,2,1>(10,N/2+1),m2);
  //solver
  Deformable<2> solver;
  solver.setMassSpring(m);
  solver.setObstacle(m2);
  solver.setDt(0);
  solver.setK(1e3);
  solver.setB(1e3);
  solver.setCL(0.2f);
  solver.setCH(1.2f);
  solver.setG(Vec2T(0,-10));

  OptimizerParam param;
  param._initBeta=initBeta;
  param._tolG=1e-2f;
  param._maxIter=1e5;
  //param._debugGradientI=50;
  param._type=OptimizerParam::ADMM;
  solver.solve(param);
}
int main(int argc, char** argv) {
  FLOAT initBeta = 1e1f;
  while(initBeta<1e3f) {
    std::string name = "MassSpring2DLongBetaProfile" + std::to_string(initBeta) + ".txt";
    if (exists(name))
      continue;
    if (freopen(name.c_str(), "w", stdout) == NULL) {
        // Handle error if redirection fails
        perror("Failed to redirect stdout");
        return 1;
    }
    std::cout << "Begin-testing " << initBeta << " beta" << std::endl;
    solve(initBeta);
    std::cout << "End-testing " << initBeta << " beta" << std::endl;
    initBeta *= 2;
  }
  return 0;
}
