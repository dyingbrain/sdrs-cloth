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
int main(int argc,char** argv) {
  //build mesh grid
  int N=100,M=10;
  MeshExact m,m2;
  Eigen::Matrix<double,3,1> ext;
  createGrid(N,M,Eigen::Matrix<double,3,1>(0,0,1),Eigen::Matrix<double,3,1>(1,0,0),m);
  createBox(Eigen::Matrix<double,3,1>(0,0,0),ext=Eigen::Matrix<double,3,1>(M/2.0+1,M/2.0+1,N/2.0+1),m2);
  //solver
  bool sim=false;
  Deformable<3> solver;
  solver.setMassSpring(m);
  solver.setObstacle(m2);
  solver.setDt(0);//(1e-1f);
  solver.setK(1e3);
  solver.setB(1e3);
  solver.setCL(0.2f);
  solver.setCH(1.2f);
  solver.setMargin(0.05f);
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

  recreate(std::filesystem::path("MassSpring3DLong"));
  int outputIter=0;
  OptimizerParam param;
  param._initBeta=5e2f;
  param._initBetaX=5e3f;
  param._tolG=1e-2f;
  param._maxIter=1e10;
  //param._debugGradientI=50;
  param._type=OptimizerParam::ADMM;
  solver.setCB([&]() {
    for(int i=0; i<(int)m.vss().size(); i++)
      m.vssNonConst()[i]=solver.x().template segment<3>(i*3).template cast<MeshExact::T>();
    VTKWriter<double> os("MassSpring3DLong","MassSpring3DLong/frame"+std::to_string(outputIter++)+".vtk",true);
    m.writeVTK(os,MeshExact::Mat3X4T::Identity());
    updateDeformable(s,solver);
    drawer.nextFrame();
    return sim;
  });
  solver.solve(param);
  return 0;
}
