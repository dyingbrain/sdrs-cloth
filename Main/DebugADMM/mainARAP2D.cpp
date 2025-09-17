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
int main(int argc,char** argv) {
  //build mesh grid
  MeshExact m,m2;
  createGrid(32,6,m);
  createGrid(32,6,m2);
  m2.translate(MeshExact::Vec3T(-16,-10,0));

  //solver
  bool sim=false;
  Deformable<2> solver;
  solver.setARAP2D(m);
  //solver.setARAP2D(m2);
  solver.fix([&](const Vec2T& pos) {
    if(pos[1]>=0)
      return pos[0]<0.5f;
    else
      return pos[0]<-15.5f || pos[0]>15.5f;
  },0);
  solver.setDt(0);
  solver.setK(1e3);
  solver.setG(Vec2T(0,-20));

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
  drawer.mainLoop();

  recreate(std::filesystem::path("ARAP2D"));
  int outputIter=0;
  OptimizerParam param;
  param._initBeta=1e2f;
  param._tolG=1e-2f;
  param._maxIter=5e4;
  param._printI=1;
  //param._debugGradientI=50;
  param._type=OptimizerParam::DIRECT_NEWTON;
  solver.setCB([&]() {
    for(int i=0; i<(int)m.vss().size(); i++)
      m.vssNonConst()[i].template segment<2>(0)=solver.x().template segment<2>(i*2).template cast<MeshExact::T>();
    VTKWriter<double> os("ARAP2D","ARAP2D/frame"+std::to_string(outputIter++)+".vtk",true);
    m.writeVTK(os,MeshExact::Mat3X4T::Identity());
    updateDeformable(s,solver);
    drawer.nextFrame();
    return sim;
  });
  solver.solve(param);
  return 0;
}
