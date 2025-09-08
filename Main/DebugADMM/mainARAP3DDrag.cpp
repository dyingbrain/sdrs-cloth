#include <ADMM/Deformable.h>
#include <ADMM/DeformableVisualizer.h>
#include <Environment/EnvironmentVisualizer.h>
#include <Environment/EnvironmentUtils.h>

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

int main(int argc,char** argv) {
  //build mesh box
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  addBox(vss,iss,Eigen::Matrix<double,3,1>(0,0,0),Eigen::Matrix<double,3,1>(10,3,3));

  //solver
  bool sim=false;
  Deformable<3> solver;
  solver.setARAP3D(MeshExact(vss,iss,true));
  solver.setDt(0);//(1e-1f);
  solver.fix([&](const Eigen::Matrix<T,3,1>& v) {
    return v[0]<0.1f;
  },0);
  solver.setK(1e3);
  solver.addL1(solver.findClosestVid(Vec3T(10,9,9)),Vec3T(5,9,9),1000);

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

  auto s=visualizeDeformable<3>(solver);
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
