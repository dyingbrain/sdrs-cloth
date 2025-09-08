#include <ADMM/Deformable.h>
#include <ADMM/DeformableVisualizer.h>
#include <Environment/EnvironmentVisualizer.h>
#include <Utils/VTKWriter.h>
#include <Utils/Utils.h>

#include <TinyVisualizer/Drawer.h>
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

void createGrid(int N,int M,MeshExact& m) {
#define ID(I,J) (I)*(M+1)+(J)
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  for(int i=0; i<=N; i++)
    for(int j=0; j<=M; j++)
      vss.push_back(Eigen::Matrix<double,3,1>(i,j,0));
  for(int i=0; i<N; i++)
    for(int j=0; j<M; j++) {
      iss.push_back(Eigen::Matrix<int,3,1>(ID(i,j),ID(i+1,j),ID(i+1,j+1)));
      iss.push_back(Eigen::Matrix<int,3,1>(ID(i,j),ID(i+1,j+1),ID(i,j+1)));
    }
  m=MeshExact(vss,iss,true);
#undef ID
}
void solve(FLOAT initBeta) {
  //build mesh grid
  MeshExact m,m2;
  createGrid(64,12,m);
  m2.translate(MeshExact::Vec3T(-16,-10,0));

  //solver
  Deformable<2> solver;
  solver.setARAP2D(m);
  solver.fix([&](const Vec2T& pos) {
    if(pos[1]>=0)
      return pos[0]<0.5f;
    else
      return pos[0]<-15.5f || pos[0]>15.5f;
  },0);
  solver.setDt(0);
  solver.setK(1e3);
  solver.setG(Vec2T(0,-20));

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
    std::string name = "ARAP2DLongBetaProfile" + std::to_string(initBeta) + ".txt";
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
