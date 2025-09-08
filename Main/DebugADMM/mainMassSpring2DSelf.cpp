#include <ADMM/Deformable.h>
#include <ADMM/DeformableVisualizer.h>
#include <Environment/EnvironmentVisualizer.h>

#include <TinyVisualizer/Drawer.h>
#include <TinyVisualizer/Camera2D.h>
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

void createGrid(int N,MeshExact& m) {
  std::vector<Eigen::Matrix<double,2,1>> vss;
  std::vector<Eigen::Matrix<int,2,1>> iss;
  for(int i=0; i<=N; i++)
    vss.push_back(Eigen::Matrix<double,2,1>(i,0));
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
int main(int argc,char** argv) {
  //build mesh grid
  int N=16;
  MeshExact m,m2,m3;
  createGrid(N,m);
  createGrid(N,m2);
  m2.translate(MeshExact::Vec3T(1,1,0));
  createBox(Eigen::Matrix<double,2,1>(0,-10),Eigen::Matrix<double,2,1>(10,1),m3);
  //solver
  bool sim=false;
  Deformable<2> solver;
  solver.setMassSpring(m);
  solver.setMassSpring(m2);
  solver.setObstacle(m3);
  solver.setDt(0);//(1e-1f);
  solver.fix(solver.findClosestVid(Deformable<2>::VecNT(0,0)),0);
  solver.fix(solver.findClosestVid(Deformable<2>::VecNT(1,1)),0);
  solver.setK(1e3);
  solver.setB(1e2);
  solver.setCL(0.2f);
  solver.setCH(1.2f);
  solver.setG(Vec2T(0,-10));

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
      param._maxIter=1e4;
      //param._debugGradientI=50;
      param._type=OptimizerParam::NEWTON;
      solver.solve(param);
      updateDeformable(s,solver);
      if(solver.dt()==0)
        sim=false;
    }
  });
  drawer.mainLoop();
  return 0;
}
