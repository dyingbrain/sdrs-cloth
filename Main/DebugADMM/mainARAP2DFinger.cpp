#include <ADMM/Deformable.h>
#include <ADMM/DeformableVisualizer.h>
#include <Environment/EnvironmentUtils.h>
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

void createBall(int N,double r,double sz,MeshExact& m) {
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  std::vector<Eigen::Matrix<double,2,1>> seeds;
  for(int i=0; i<N; i++) {
    vss.push_back(Eigen::Matrix<double,3,1>(cos(M_PI*2*i/N),sin(M_PI*2*i/N),0)*r);
    iss.push_back(Eigen::Matrix<int,3,1>(i,(i+1)%N,0));
  }
  generateTriMesh(vss,iss,seeds,sz);
  m.init(vss,iss,true);
}
void createFinger(int N,double thick,double length,double pitLength,double pitDepth,double sz,bool flip,MeshExact& m) {
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  std::vector<Eigen::Matrix<double,2,1>> seeds;
  //vss
  vss.push_back(Eigen::Matrix<double,3,1>(0,0,0));
  vss.push_back(Eigen::Matrix<double,3,1>(N*length,0,0));
  vss.push_back(Eigen::Matrix<double,3,1>(N*length,thick,0));
  for(int i=N-1; i>0; i--) {
    vss.push_back(Eigen::Matrix<double,3,1>(i*length+pitLength/2,thick,0));
    vss.push_back(Eigen::Matrix<double,3,1>(i*length,thick-pitDepth,0));
    vss.push_back(Eigen::Matrix<double,3,1>(i*length-pitLength/2,thick,0));
  }
  vss.push_back(Eigen::Matrix<double,3,1>(0,thick,0));
  if(flip)
    for(auto& v:vss)
      v[1]*=-1;
  //iss
  for(int i=0; i<(int)vss.size(); i++)
    iss.push_back(Eigen::Matrix<int,3,1>(i,(i+1)%(int)vss.size(),0));
  //seed
  generateTriMesh(vss,iss,seeds,sz);
  m.init(vss,iss,true);
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
  int N=3;
  double radius=2;
  MeshExact m1,m2,ball,obs;
  double thick=1,length=2,pitLength=1,pitDepth=0.6,sz=.3f;
  //finger
  {
    createFinger(N,thick,length,pitLength,pitDepth,sz,false,m1);
    VTKWriter<double> os("Finger","fingerA.vtk",true);
    m1.writeVTK(os,MeshExact::Mat3X4T::Identity());
  }
  {
    createFinger(N,thick,length,pitLength,pitDepth,sz,true,m2);
    m2.translate(MeshExact::Vec3T(0,(thick+radius)*2,0));
    VTKWriter<double> os("Finger","fingerB.vtk",true);
    m2.writeVTK(os,MeshExact::Mat3X4T::Identity());
  }
  //ball
  {
    createBall(32,radius*0.99,sz,ball);
    ball.translate(MeshExact::Vec3T(length*(N-1),thick+radius,0));
  }
  //obs
  {
    createBox(Eigen::Matrix<double,2,1>(sz,thick+radius),Eigen::Matrix<double,2,1>(sz,radius*0.98),obs);
  }

  //solver
  bool sim=false;
  Deformable<2> solver;
  solver.setK(1e3);
  solver.setARAP2D(m1);
  solver.setARAP2D(m2);
  solver.setK(1e1);
  solver.setARAP2D(ball);
  solver.setObstacle(obs);
  solver.fix([&](const Vec2T& pos) {
    return pos[0]<=0.01f;
  },0);
  solver.setDt(0);
  solver.setMargin(.02f);
  solver.setCollCoef(1e3);
  solver.setG(Vec2T(0,0));
  for(int i=1; i<N; i++)
    solver.addL1({solver.findClosestVid(Vec2T(i*length-pitLength/2,thick)),
                  solver.findClosestVid(Vec2T(i*length+pitLength/2,thick))},3e3);
  for(int i=1; i<N; i++)
    solver.addL1({solver.findClosestVid(Vec2T(i*length-pitLength/2,thick+radius*2)),
                  solver.findClosestVid(Vec2T(i*length+pitLength/2,thick+radius*2))},3e3);

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

  recreate(std::filesystem::path("Finger2D"));
  int outputIter=0;
  OptimizerParam param;
  param._initBeta=1e3f;
  param._initBetaX=1e3f;
  param._tolG=1e-2f;
  param._maxIter=1e5;
  //param._debugGradientI=50;
  param._type=OptimizerParam::ADMM;
  solver.setCB([&]() {
    int off=0;
    //fingerA
    {
      for(int i=0; i<(int)m1.vss().size(); i++,off++)
        m1.vssNonConst()[i].template segment<2>(0)=solver.x().template segment<2>(off*2).template cast<MeshExact::T>();
      VTKWriter<double> os("Finger","Finger2D/fingerA"+std::to_string(outputIter)+".vtk",true);
      m1.writeVTK(os,MeshExact::Mat3X4T::Identity());
    }
    //fingerB
    {
      for(int i=0; i<(int)m2.vss().size(); i++,off++)
        m2.vssNonConst()[i].template segment<2>(0)=solver.x().template segment<2>(off*2).template cast<MeshExact::T>();
      VTKWriter<double> os("Finger","Finger2D/fingerB"+std::to_string(outputIter)+".vtk",true);
      m2.writeVTK(os,MeshExact::Mat3X4T::Identity());
    }
    //ball
    {
      for(int i=0; i<(int)ball.vss().size(); i++,off++)
        ball.vssNonConst()[i].template segment<2>(0)=solver.x().template segment<2>(off*2).template cast<MeshExact::T>();
      VTKWriter<double> os("Ball","Finger2D/ball"+std::to_string(outputIter)+".vtk",true);
      ball.writeVTK(os,MeshExact::Mat3X4T::Identity());
    }
    //cable
    {
      off=0;
      VTKWriter<double> os("Base","Finger2D/cable"+std::to_string(outputIter)+".vtk",true);
      std::vector<Eigen::Matrix<double,3,1>> vss;
      std::vector<Eigen::Matrix<int,2,1>> iss;
      Eigen::Matrix<double,3,1> pos;
      pos.setZero();
      for(auto l1:solver.getL1ss()) {
        pos.template segment<2>(0)=solver.x().template segment<2>(l1.first[0]*2).template cast<double>();
        vss.push_back(pos);
        pos.template segment<2>(0)=solver.x().template segment<2>(l1.first[1]*2).template cast<double>();
        vss.push_back(pos);
        iss.push_back(Eigen::Matrix<int,2,1>(off,off+1));
        off+=2;
      }
      os.appendPoints(vss.begin(),vss.end());
      os.appendCells(iss.begin(),iss.end(),VTKWriter<double>::LINE);
    }
    outputIter++;
    //draw
    updateDeformable(s,solver);
    drawer.nextFrame();
    return sim;
  });
  solver.solve(param);
  return 0;
}
