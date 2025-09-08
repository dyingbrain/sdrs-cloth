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

void createGrid(int N,MeshExact& m,T rand=0) {
#define ID(I,J) (I)*(N+1)+(J)
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  for(int i=0; i<=N; i++)
    for(int j=0; j<=N; j++) {
      vss.push_back(Eigen::Matrix<double,3,1>(i,j,0));
      if(rand>0)
        vss.back()+=Eigen::Matrix<double,3,1>::Random()*(double)rand;
    }
  for(int i=0; i<N; i++)
    for(int j=0; j<N; j++) {
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
int main(int argc,char** argv) {
  //build mesh grid
  int N=16;
  MeshExact m,m2;
  createGrid(N,m);
  createBox(Eigen::Matrix<double,3,1>(N/2,N/2,-5),Eigen::Matrix<double,3,1>(10,10,1),m2);
  //solver
  bool sim=false;
  Deformable<3> solver;
  solver.setMassSpring(m);
  solver.setObstacle(m2);
  solver.setDt(0);//(1e-1f);
  solver.fix(0,1e4);
  solver.fix(N,1e4);
  solver.setK(1e3);
  solver.setB(1e2);
  solver.setCL(0.2f);
  solver.setCH(1.2f);
  solver.setG(Vec3T(0,0,-10));

  //draw
  Drawer drawer(argc,argv);
  std::shared_ptr<CaptureGIFPlugin> ss(new CaptureGIFPlugin(GLFW_KEY_4,"screenshot.gif",drawer.FPS(),true));
  drawer.addPlugin(std::shared_ptr<Plugin>(new CameraExportPlugin(GLFW_KEY_2,GLFW_KEY_3,"camera.dat")));
  drawer.addPlugin(std::shared_ptr<Plugin>(new CaptureGIFPlugin(GLFW_KEY_1,"record.gif",drawer.FPS())));
#define USE_LIGHT
#ifdef USE_LIGHT
  drawer.addLightSystem(2048,20);
  drawer.getLight()->lightSz(10);
  drawer.getLight()->addLight(Eigen::Matrix<GLfloat,3,1>(2,2,2),
                              Eigen::Matrix<GLfloat,3,1>(1,1,1),
                              Eigen::Matrix<GLfloat,3,1>(1,1,1),
                              Eigen::Matrix<GLfloat,3,1>(0,0,0));
  drawer.getLight()->autoAdjust(true);
#endif

  auto s=visualizeDeformable(solver);
  drawer.addShape(s);
  drawer.addCamera3D(90,Eigen::Matrix<GLfloat,3,1>(0,0,1));
  drawer.getCamera3D()->setManipulator(std::shared_ptr<CameraManipulator>(new FirstPersonCameraManipulator(drawer.getCamera3D())));
  drawer.addPlugin(std::shared_ptr<Plugin>(new ImGuiPlugin([&]() {
    drawer.getCamera3D()->getManipulator()->imGuiCallback();
  })));
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
