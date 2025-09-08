#include <ADMM/Deformable.h>
#include <ADMM/DeformableVisualizer.h>
#include <Environment/EnvironmentVisualizer.h>
#include <Utils/VTKWriter.h>
#include <Utils/Utils.h>

#include <TinyVisualizer/Drawer.h>
#include <TinyVisualizer/Povray.h>
#include <TinyVisualizer/Camera3D.h>
#include <TinyVisualizer/MakeMesh.h>
#include <TinyVisualizer/ShadowAndLight.h>
#include <TinyVisualizer/FirstPersonCameraManipulator.h>
#include <TinyVisualizer/CameraExportPlugin.h>
#include <TinyVisualizer/CaptureGIFPlugin.h>
#include <TinyVisualizer/ImGuiPlugin.h>

using namespace PHYSICSMOTION;
using namespace DRAWER;
typedef FLOAT T;
DECL_MAT_VEC_MAP_TYPES_T

void createGrid(int N,int M,Eigen::Matrix<double,3,1> dN,Eigen::Matrix<double,3,1> dM,MeshExact& m) {
#define ID(I,J) (I)*(M+1)+(J)
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  for(int i=0; i<=N; i++)
    for(int j=0; j<=M; j++) {
      vss.push_back(dN*(i-N/2.0)+dM*(j-M/2.0));
    }
  for(int i=0; i<N; i++)
    for(int j=0; j<M; j++) {
      iss.push_back(Eigen::Matrix<int,3,1>(ID(i,j),ID(i+1,j),ID(i+1,j+1)));
      iss.push_back(Eigen::Matrix<int,3,1>(ID(i,j),ID(i+1,j+1),ID(i,j+1)));
    }
  m.init<double>(vss,iss);
#undef ID
}
void createBox(Eigen::Matrix<double,3,1> ctr,Eigen::Matrix<double,3,1> sz,MeshExact& m) {
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  BBoxExact bb(sz[0],sz[1],sz[2]);
  bb.getMesh(vss,iss);
  for(auto& v:vss)
    v+=ctr;
  m=MeshExact(vss,iss,true);
}
void solve(int N,int M=5) {
  //build mesh grid
  MeshExact m,m2;
  Eigen::Matrix<double,3,1> ext;
  createGrid(N,M,Eigen::Matrix<double,3,1>(0,0,1),Eigen::Matrix<double,3,1>(1,0,0),m);
  createBox(Eigen::Matrix<double,3,1>(0,0,0),ext=Eigen::Matrix<double,3,1>(M/2.0+1,M/2.0+1,N/2.0+1),m2);
  //solver
  Deformable<3> solver;
  solver.setMassSpring(m);
  solver.setObstacle(m2);
  solver.setDt(0);//(1e-1f);
  solver.setK(1e3);
  solver.setB(1e3);
  solver.setCL(0.2f);
  solver.setCH(1.2f);
  solver.setG(Vec3T(0,0,-10));

  OptimizerParam param;
  param._initBeta=5e2f;
  param._initBetaX=5e3f;
  param._tolG=1e-2f;
  param._maxIter=2e5;
  //param._debugGradientI=50;
  param._type=OptimizerParam::ADMM;
  solver.solve(param);
}
int main(int argc, char** argv) {
  for (int n : {25,50,75,100,125}) {
    std::string name = "MassSpring3DLongProfile" + std::to_string(n) + ".txt";
    if (exists(name))
      continue;
    if (freopen(name.c_str(), "w", stdout) == NULL) {
        // Handle error if redirection fails
        perror("Failed to redirect stdout");
        return 1;
    }
    std::cout << "Begin-testing " << n << " grids" << std::endl;
    solve(n);
    std::cout << "End-testing " << n << " grids" << std::endl;
  }
  return 0;
}
