#include <ADMM/Deformable.h>
#include <ADMM/DeformableVisualizer.h>
#include <Environment/EnvironmentVisualizer.h>
#include <Utils/VTKWriter.h>
#include <Utils/Utils.h>

#include <TinyVisualizer/Drawer.h>
#include <TinyVisualizer/Povray.h>
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
int main(int argc,char** argv) {
  //build mesh grid
  int N=100,D=100;
  MeshExact m,m2;
  createGrid(N,Eigen::Matrix<double,2,1>(0,0),Eigen::Matrix<double,2,1>(1,D),m);
  createBox(Eigen::Matrix<double,2,1>(0,N/2),Eigen::Matrix<double,2,1>(10,N/2+1),m2);
  //solver
  bool sim=false;
  Deformable<2> solver;
  solver.setMassSpring(m);
  solver.setObstacle(m2);
  solver.setDt(0);
  solver.setK(1e3);
  solver.setB(1e3);
  solver.setCL(0.2f);
  solver.setCH(1.2f);
  solver.setG(Vec2T(0,-10));

  //draw
  Drawer drawer(argc,argv);
  std::shared_ptr<CaptureGIFPlugin> ss(new CaptureGIFPlugin(GLFW_KEY_4,"screenshot.gif",drawer.FPS(),true));
  drawer.addPlugin(std::shared_ptr<Plugin>(new CameraExportPlugin(GLFW_KEY_2,GLFW_KEY_3,"camera.dat")));
  drawer.addPlugin(std::shared_ptr<Plugin>(new CaptureGIFPlugin(GLFW_KEY_1,"record.gif",drawer.FPS())));
  std::shared_ptr<CompositeShape> s=visualizeDeformable(solver,10);
  drawer.addShape(s);
  drawer.clearLight();
  drawer.addCamera2D(20);
  s->getChild(0)->setColorAmbient(GL_LINES,1,.5f,0);
  s->getChild(s->numChildren()-1)->setColorAmbient(GL_LINES,.5f,.5f,.5f);
  drawer.setKeyFunc([&](GLFWwindowPtr wnd,int key,int scan,int action,int mods,bool captured) {
    if(captured)
      return;
    else if(key==GLFW_KEY_R && action==GLFW_PRESS)
      sim=!sim;
  });
  drawer.mainLoop();

  recreate(std::filesystem::path("MassSpring2DLong"));
  int outputIter=0;
  OptimizerParam param;
  param._initBeta=1e2f;
  param._tolG=1e-2f;
  param._maxIter=1e5;
  //param._debugGradientI=50;
  param._type=OptimizerParam::ADMM;
  solver.setCB([&]() {
    for(int i=0; i<(int)m.vss().size(); i++)
      m.vssNonConst()[i].template segment<2>(0)=solver.x().template segment<2>(i*2).template cast<MeshExact::T>();
    VTKWriter<double> os("MassSpring2DLong","MassSpring2DLong/frame"+std::to_string(outputIter++)+".vtk",true);
    m.writeVTK(os,MeshExact::Mat3X4T::Identity());
    updateDeformable(s,solver);
    drawer.nextFrame();
    return sim;
  });
  solver.solve(param);
  return 0;
}
