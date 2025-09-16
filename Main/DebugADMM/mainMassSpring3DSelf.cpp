#include <ADMM/Deformable.h>
#include <ADMM/DeformableVisualizer.h>
#include <Environment/EnvironmentVisualizer.h>
#include <Utils/VTKWriter.h>

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

void createGrid(int N,MeshExact& m,double rand=0) {
#define ID(I,J) (I)*(N+1)+(J)
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  for(int i=0; i<=N; i++)
    for(int j=0; j<=N; j++) {
      vss.push_back(Eigen::Matrix<double,3,1>(i,j,0));
      if(rand>0)
        vss.back()+=Eigen::Matrix<double,3,1>::Random()*rand;
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
void testStretch(int argc,char** argv,double alpha) {
  //build mesh grid
  int N=16;
  MeshExact m,m2,m3;
  createGrid(N,m);
  createGrid(N,m2);
  m2.translate(Vec3T(1,0,1));
  createBox(Eigen::Matrix<double,3,1>(N/2,N/2,-5),Eigen::Matrix<double,3,1>(10,10,1),m3);
  //solver
  bool sim=false;
  Deformable<3> solver;
  solver.setK(1e3);
  solver.setB(1e2);
  solver.setCL(1-alpha);
  solver.setCH(1+alpha);
  solver.setMassSpring(m);
  //solver.setMassSpring(m2);
  //solver.setObstacle(m3);
  solver.setDt(0);//(1e-1f);
  solver.fix(solver.findClosestVid(Deformable<3>::VecNT(0,0,0)),0);
  solver.fix(solver.findClosestVid(Deformable<3>::VecNT(0,N,0)),0);
  //solver.fix(solver.findClosestVid(Deformable<3>::VecNT(1,0,1)),1e4);
  //solver.fix(solver.findClosestVid(Deformable<3>::VecNT(1,N,1)),1e4);
  solver.setG(Vec3T(0,0,-100));

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
  drawer.mainLoop();

  OptimizerParam param;
  param._initBeta=1e2f;
  param._tolG=1e-2f;
  param._maxIter=1e5;
  param._printI=1;
  param._debugGradientI=1;
  param._type=OptimizerParam::DIRECT_NEWTON;
  solver.setCB([&]() {
    updateDeformable(s,solver);
    drawer.nextFrame();
    return sim;
  });
  solver.solve(param);
  for(int i=0; i<(int)m.vss().size(); i++)
    m.vssNonConst()[i]=solver.x().template segment<3>(i*3).template cast<MeshExact::T>();
  VTKWriter<double> os("MassSpring3DSelf","MassSpringStretch"+std::to_string(alpha)+".vtk",true);
  m.writeVTK(os,MeshExact::Mat3X4T::Identity());
}
int main(int argc,char** argv) {
  testStretch(argc,argv,0.4f);
  testStretch(argc,argv,0.2f);
  testStretch(argc,argv,0.1f);
  return 0;
}
