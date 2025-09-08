#include <ADMM/Deformable.h>
#include <ADMM/DeformableVisualizer.h>
#include <Environment/EnvironmentVisualizer.h>

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

#define ID(I,J) (I)*(M+1)+(J)
void createGrid(int N,int M,MeshExact& m) {
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
}
int main(int argc,char** argv) {
  //build mesh grid
  MeshExact m;
  int N=16,M=2;
  createGrid(N,M,m);

  //solver
  bool sim=false;
  Deformable<2> solver;
  solver.setARAP2D(m);
  solver.setDt(0);//(1e-1f);
  for(int j=0; j<=M; j++)
    solver.fix(ID(0,j),1e4);
  solver.setK(1e3);
  solver.addL1(ID(N,M),Vec2T(8,8),1000);

  //draw
  Drawer drawer(argc,argv);
  std::shared_ptr<CaptureGIFPlugin> ss(new CaptureGIFPlugin(GLFW_KEY_4,"screenshot.gif",drawer.FPS(),true));
  drawer.addPlugin(std::shared_ptr<Plugin>(new CameraExportPlugin(GLFW_KEY_2,GLFW_KEY_3,"camera.dat")));
  drawer.addPlugin(std::shared_ptr<Plugin>(new CaptureGIFPlugin(GLFW_KEY_1,"record.gif",drawer.FPS())));

  auto s=visualizeDeformable(solver);
  drawer.addShape(s);
  drawer.addCamera2D(20);
  drawer.setKeyFunc([&](GLFWwindowPtr wnd,int key,int scan,int action,int mods,bool captured) {
    if(captured)
      return;
    else if(key==GLFW_KEY_R && action==GLFW_PRESS)
      sim=!sim;
  });
  drawer.setFrameFunc([&](std::shared_ptr<SceneNode>& node) {
    if(sim) {
      OptimizerParam param;
      param._initBeta=1e2f;
      param._tolG=1e-2f;
      param._maxIter=1e6;
      param._debugGradientI=50;
      param._type=OptimizerParam::ADMM;
      solver.solve(param);
      updateDeformable(s,solver);
      if(solver.dt()==0)
        sim=false;
    }
  });
  drawer.mainLoop();
  return 0;
}
