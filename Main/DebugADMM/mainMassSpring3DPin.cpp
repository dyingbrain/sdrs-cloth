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

void createGrid(int sz,int theta,Eigen::Matrix<double,3,1> ctr,MeshExact& m) {
#define ID(I,J) (I)*(sz+1)+(J)
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  auto R=Eigen::AngleAxis<double>(theta*EIGEN_PI/180.0,Eigen::Matrix<double,3,1>::UnitZ()).toRotationMatrix();
  for(int i=0; i<=sz; i++)
    for(int j=0; j<=sz; j++) {
      Eigen::Matrix<double,3,1> v(i-sz/2.0, j-sz/2.0, 0);
      vss.push_back(R*v+ctr);
    }
  for(int i=0; i<sz; i++)
    for(int j=0; j<sz; j++) {
      iss.push_back(Eigen::Matrix<int,3,1>(ID(i,j),ID(i+1,j),ID(i+1,j+1)));
      iss.push_back(Eigen::Matrix<int,3,1>(ID(i,j),ID(i+1,j+1),ID(i,j+1)));
    }
  m.init<double>(vss,iss);
#undef ID
}
void createPin(Eigen::Matrix<double,3,1> ctr,Eigen::Matrix<double,3,1> sz,MeshExact& m) {
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  vss.push_back(Eigen::Matrix<double,3,1>(-sz[0],-sz[1],-sz[2])+ctr);
  vss.push_back(Eigen::Matrix<double,3,1>( sz[0],-sz[1],-sz[2])+ctr);
  vss.push_back(Eigen::Matrix<double,3,1>( sz[0], sz[1],-sz[2])+ctr);
  vss.push_back(Eigen::Matrix<double,3,1>(-sz[0], sz[1],-sz[2])+ctr);
  vss.push_back(Eigen::Matrix<double,3,1>(0,0,sz[2])+ctr);
  iss.push_back(Eigen::Matrix<int,3,1>(0,1,4));
  iss.push_back(Eigen::Matrix<int,3,1>(1,2,4));
  iss.push_back(Eigen::Matrix<int,3,1>(2,3,4));
  iss.push_back(Eigen::Matrix<int,3,1>(3,0,4));
  //iss.push_back(Eigen::Matrix<int,3,1>(0,1,2));
  //iss.push_back(Eigen::Matrix<int,3,1>(0,2,3));
  m.init<double>(vss,iss);
}
void solve(int ClothSize,int argc,char** argv,bool vis,OptimizerParam::TYPE type,const std::string& typeString) {
  //build mesh grid
  int BaseSize=1;
  int PinSeparation=10;
  int lmt=((ClothSize/4+PinSeparation-1)/PinSeparation)*PinSeparation;
  MeshExact m,m2,m3;
  //solver
  bool sim=false;
  Deformable<3> solver;
  createGrid(ClothSize,0,Eigen::Matrix<double,3,1>(0,0,5),m);
  solver.setMassSpring(m);
  createPin(Eigen::Matrix<double,3,1>(0,0,-1),Eigen::Matrix<double,3,1>(BaseSize/2.0,BaseSize/2.0,2),m2);
  for(int x=-lmt;x<=lmt;x+=PinSeparation) {
    for(int y=-lmt;y<=lmt;y+=PinSeparation) {
      createPin(Eigen::Matrix<double,3,1>(x,y,-1),Eigen::Matrix<double,3,1>(BaseSize/2.0,BaseSize/2.0,2),m3);
      m2.addMesh(m3);
    }
  }
  m2.init<double>(m2.vss(),m2.iss());
  solver.setObstacle(m2);
  solver.setDt(0);//(1e-1f);
  solver.setK(1e3);
  solver.setB(1e3);
  solver.setCL(0.2f);
  solver.setCH(1.2f);
  solver.setG(Vec3T(0,0,-10));
  {
    VTKWriter<double> os("Obstacle", "Obstacle.vtk", true);
    m2.writeVTK(os, MeshExact::Mat3X4T::Identity());
  }

  std::shared_ptr<Drawer> drawer;
  std::shared_ptr<CompositeShape> s;
  if(vis) {
    //draw
    drawer.reset(new Drawer(argc,argv));
    std::shared_ptr<CaptureGIFPlugin> ss(new CaptureGIFPlugin(GLFW_KEY_4,"screenshot.gif",drawer->FPS(),true));
    drawer->addPlugin(std::shared_ptr<Plugin>(new CameraExportPlugin(GLFW_KEY_2,GLFW_KEY_3,"camera.dat")));
    drawer->addPlugin(std::shared_ptr<Plugin>(new CaptureGIFPlugin(GLFW_KEY_1,"record.gif",drawer->FPS())));
#define USE_LIGHT
#ifdef USE_LIGHT
    drawer->addLightSystem(2048,20);
    drawer->getLight()->lightSz(10);
    drawer->getLight()->addLight(Eigen::Matrix<GLfloat,3,1>(2,2,2),
                                 Eigen::Matrix<GLfloat,3,1>(1,1,1),
                                 Eigen::Matrix<GLfloat,3,1>(1,1,1),
                                 Eigen::Matrix<GLfloat,3,1>(0,0,0));
    drawer->getLight()->autoAdjust(true);
#endif

    s=visualizeDeformable(solver);
    drawer->addShape(s);
    drawer->addCamera3D(90,Eigen::Matrix<GLfloat,3,1>(0,0,1));
    drawer->getCamera3D()->setManipulator(std::shared_ptr<CameraManipulator>(new FirstPersonCameraManipulator(drawer->getCamera3D())));
    drawer->addPlugin(std::shared_ptr<Plugin>(new ImGuiPlugin([&]() {
      drawer->getCamera3D()->getManipulator()->imGuiCallback();
    })));
    drawer->setKeyFunc([&](GLFWwindowPtr wnd,int key,int scan,int action,int mods,bool captured) {
      if(captured)
        return;
      else if(key==GLFW_KEY_R && action==GLFW_PRESS)
        sim=!sim;
    });
    drawer->mainLoop();
  } else {
    std::string name = "MassSpring3DPin"+typeString+std::to_string(ClothSize) + ".txt";
    if (exists(name))
      return;
    if (freopen(name.c_str(), "w", stdout) == NULL) {
      // Handle error if redirection fails
      perror("Failed to redirect stdout");
      return;
    }
  }

  int outputIter=0;
  OptimizerParam param;
  param._initBeta=5e2f;
  param._initBetaX=5e3f;
  param._tolG=1e-2f;
  param._maxIter=1e5;
  //param._debugGradientI=50;
  param._type=type;
  if(vis) {
    recreate(std::filesystem::path("MassSpring3DPin" + typeString));
    solver.setCB([&]() {
      for(int i=0; i<(int)m.vss().size(); i++)
        m.vssNonConst()[i]=solver.x().template segment<3>(i*3).template cast<MeshExact::T>();
      VTKWriter<double> os("MassSpring3DPin","MassSpring3DPin"+typeString+"/frame"+std::to_string(outputIter++)+".vtk",true);
      m.writeVTK(os,MeshExact::Mat3X4T::Identity());
      updateDeformable(s, solver);
      drawer->nextFrame();
      return sim;
    });
  }
  solver.solve(param);
}
int main(int argc,char** argv) {
  for (int n : {30, 60, 90, 120}) {
    solve(n,argc,argv,false,OptimizerParam::TYPE::ADMM,"ADMM");
    solve(n,argc,argv,false,OptimizerParam::TYPE::GD,"GD");
    solve(n,argc,argv,false,OptimizerParam::TYPE::NEWTON,"NEWTON");
  }
  return 0;
}