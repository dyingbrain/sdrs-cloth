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
int main(int argc,char** argv) {
  int N=20;
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
  bool sim=false;
  Deformable<3> solver;
  solver.setObstacle(obs);
  solver.setDt(0);
  solver.setR(0);
  solver.setG(Vec3T(0,0,0));
  ASSERT_MSG(solver.setUAVTrajectory(waypoints,(int)waypoints.size()-1,(int)waypoints.size()*5),"Failed initialization!")
  solver.fix(0,0);
  solver.fix(solver.x().size()/3-1,0);
  std::cout << "Using " << solver.getUAVTraj().getTraj().getNumSegment() << " segments!" << std::endl;

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

  std::shared_ptr<CompositeShape> s=visualizeDeformable(solver);
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
  drawer.mainLoop();

  recreate(std::filesystem::path("UAV"));
  {
    VTKWriter<double> os("Base","UAV/obs.vtk",true);
    obs.writeVTK(os,MeshExact::Mat3X4T::Identity());
  }
  int outputIter=0;
  OptimizerParam param;
  param._initBeta=1e2f;
  param._initBetaX=0;
  param._tolG=1e-2f;
  param._maxIter=1e5;
  //param._debugGradientI=50;
  param._type=OptimizerParam::ADMM;
  solver.setCB([&]() {
    Deformable<3>::Vec CPs=solver.getTrajSubd()*solver.x();
    solver.getUAVTraj().writeTrajVTK("UAV/frame"+std::to_string(outputIter++)+".vtk",&CPs,0);
    updateDeformable(s,solver);
    drawer.nextFrame();
    return sim;
  });
  solver.solve(param);
  return 0;
}
