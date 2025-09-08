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

void createFinger(int N,double thick,double width,double length,double pitLength,double pitDepth,double sz,bool flip,MeshExact& m) {
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
  generateTriMesh(vss,iss,seeds,1);
  extrude(vss,iss,width/2);
  m=MeshExact(vss,iss,true);
}
void createBase(double radius,double sz,MeshExact& m) {
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss({Eigen::Matrix<int,3,1>(0,1,2)});
  vss.push_back(Eigen::Matrix<double,3,1>(0,2*radius,0));
  vss.push_back(Eigen::Matrix<double,3,1>(-std::sqrt(3.0)*radius,-radius,0));
  vss.push_back(Eigen::Matrix<double,3,1>( std::sqrt(3.0)*radius,-radius,0));
  extrude(vss,iss,sz);
  m.init(vss,iss,true);
  MeshExact::Mat3X4T t=MeshExact::Mat3X4T::Identity();
  t.template block<3,3>(0,0)=Eigen::AngleAxisd(M_PI/2,Eigen::Matrix<double,3,1>::UnitY()).toRotationMatrix().template cast<MeshExact::T>();
  m.transform(t);
  m.translate(MeshExact::Vec3T(sz,0,0));
}
int main(int argc,char** argv) {
  int N=3;
  double radius=2;
  MeshExact::Mat3X4T Rt[3];
  MeshExact m,m1,m2,m3,ball,obs;
  double thick=1.5,width=1.5,length=2,pitLength=1.5,pitDepth=0.6,sz=.3f;
  //finger
  createFinger(N,thick,width,length,pitLength,pitDepth,sz,false,m);
  //fingerA
  {
    m1=m;
    MeshExact::Vec3T t(0,-thick-radius,0);
    MeshExact::Mat3T R=Eigen::AngleAxisd(M_PI*0,Eigen::Matrix<double,3,1>::UnitX()).toRotationMatrix().template cast<MeshExact::T>();
    Rt[0].template block<3,3>(0,0)=R;
    Rt[0].col(3)=R*t;
    m1.transform(Rt[0]);
    //write
    VTKWriter<double> os("Finger","fingerA.vtk",true);
    m1.writeVTK(os,MeshExact::Mat3X4T::Identity());
  }
  //fingerB
  {
    m2=m;
    MeshExact::Vec3T t(0,-thick-radius,0);
    MeshExact::Mat3T R=Eigen::AngleAxisd(M_PI*120/180,Eigen::Matrix<double,3,1>::UnitX()).toRotationMatrix().template cast<MeshExact::T>();
    Rt[1].template block<3,3>(0,0)=R;
    Rt[1].col(3)=R*t;
    m2.transform(Rt[1]);
    VTKWriter<double> os("Finger","fingerB.vtk",true);
    m2.writeVTK(os,MeshExact::Mat3X4T::Identity());
  }
  //fingerC
  {
    m3=m;
    MeshExact::Vec3T t(0,-thick-radius,0);
    MeshExact::Mat3T R=Eigen::AngleAxisd(M_PI*240/180,Eigen::Matrix<double,3,1>::UnitX()).toRotationMatrix().template cast<MeshExact::T>();
    Rt[2].template block<3,3>(0,0)=R;
    Rt[2].col(3)=R*t;
    m3.transform(Rt[2]);
    VTKWriter<double> os("Finger","fingerC.vtk",true);
    m3.writeVTK(os,MeshExact::Mat3X4T::Identity());
  }
  //ball
  {
    std::vector<Eigen::Matrix<double,3,1>> vss;
    std::vector<Eigen::Matrix<int,3,1>> iss;
    addSphere(vss,iss,Eigen::Matrix<double,3,1>(length*(N-1),0,0),radius*0.99,4);
    ball.init(vss,iss,true);
    VTKWriter<double> os("Ball","ball.vtk",true);
    ball.writeVTK(os,MeshExact::Mat3X4T::Identity());
  }
  //obs
  {
    createBase(radius*0.9,sz,obs);
  }

  //solver
  int offV1=0,offV2=0,offV3=0,offVB=0;
  int offT1=0,offT2=0,offT3=0,offTB=0;
  bool sim=false;
  Deformable<3> solver;
  solver.setK(1e3);
  //fingerA
  solver.setARAP3D(m1);
  offV1=solver.x().size()/3;
  offT1=solver.getTss().size();
  //fingerB
  solver.setARAP3D(m2);
  offV2=solver.x().size()/3;
  offT2=solver.getTss().size();
  //fingerC
  solver.setARAP3D(m3);
  offV3=solver.x().size()/3;
  offT3=solver.getTss().size();
  solver.setK(1e1);
  //ball
  solver.setARAP3D(ball,91,0.5);
  offVB=solver.x().size()/3;
  offTB=solver.getTss().size();
  solver.setObstacle(obs);
  solver.fix([&](const Vec3T& pos) {
    return pos[0]<=0.01f;
  },0);
  solver.setDt(0);
  solver.setMargin(.02f);
  solver.setCollCoef(1e3);
  solver.setG(Vec3T(0,0,0));
  for(int d=0; d<3; d++) {
    for(int i=1; i<N; i++) {
      solver.addL1({solver.findClosestVid(Rt[d].template block<3,3>(0,0).template cast<T>()*Vec3T(i*length-pitLength/2,thick,-width/2)+Rt[d].col(3).template cast<T>()),
                    solver.findClosestVid(Rt[d].template block<3,3>(0,0).template cast<T>()*Vec3T(i*length+pitLength/2,thick,-width/2)+Rt[d].col(3).template cast<T>())},3e4);
      solver.addL1({solver.findClosestVid(Rt[d].template block<3,3>(0,0).template cast<T>()*Vec3T(i*length-pitLength/2,thick, width/2)+Rt[d].col(3).template cast<T>()),
                    solver.findClosestVid(Rt[d].template block<3,3>(0,0).template cast<T>()*Vec3T(i*length+pitLength/2,thick, width/2)+Rt[d].col(3).template cast<T>())},3e4);
    }
  }

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
  drawer.addCamera3D(90,Eigen::Matrix<GLfloat,3,1>(1,0,0));
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

  recreate(std::filesystem::path("Finger3D"));
  {
    VTKWriter<double> os("Base","Finger3D/obs.vtk",true);
    obs.writeVTK(os,MeshExact::Mat3X4T::Identity());
  }
  int outputIter=0;
  OptimizerParam param;
  param._initBeta=5e4f;
  param._initBetaX=5e4f;
  param._tolG=1e-2f;
  param._maxIter=3e5;
  //param._debugGradientI=50;
  param._type=OptimizerParam::ADMM;
  solver.setCB([&]() {
    //fingerA
    {
      std::vector<Eigen::Matrix<double,3,1>> vss;
      std::vector<Eigen::Matrix<int,4,1>> tss;
      for(int i=0; i<offV1; i++)
        vss.push_back(solver.x().template segment<3>(i*3).template cast<double>());
      for(int i=0; i<offT1; i++)
        tss.push_back(solver.getTss()[i]._indices-Eigen::Matrix<int,4,1>::Constant(0));
      VTKWriter<double> os("Finger","Finger3D/fingerA"+std::to_string(outputIter)+".vtk",true);
      os.appendPoints(vss.begin(),vss.end());
      os.appendCells(tss.begin(),tss.end(),VTKWriter<double>::TETRA);
    }
    //fingerB
    {
      std::vector<Eigen::Matrix<double,3,1>> vss;
      std::vector<Eigen::Matrix<int,4,1>> tss;
      for(int i=offV1; i<offV2; i++)
        vss.push_back(solver.x().template segment<3>(i*3).template cast<double>());
      for(int i=offT1; i<offT2; i++)
        tss.push_back(solver.getTss()[i]._indices-Eigen::Matrix<int,4,1>::Constant(offV1));
      VTKWriter<double> os("Finger","Finger3D/fingerB"+std::to_string(outputIter)+".vtk",true);
      os.appendPoints(vss.begin(),vss.end());
      os.appendCells(tss.begin(),tss.end(),VTKWriter<double>::TETRA);
    }
    //fingerC
    {
      std::vector<Eigen::Matrix<double,3,1>> vss;
      std::vector<Eigen::Matrix<int,4,1>> tss;
      for(int i=offV2; i<offV3; i++)
        vss.push_back(solver.x().template segment<3>(i*3).template cast<double>());
      for(int i=offT2; i<offT3; i++)
        tss.push_back(solver.getTss()[i]._indices-Eigen::Matrix<int,4,1>::Constant(offV2));
      VTKWriter<double> os("Finger","Finger3D/fingerC"+std::to_string(outputIter)+".vtk",true);
      os.appendPoints(vss.begin(),vss.end());
      os.appendCells(tss.begin(),tss.end(),VTKWriter<double>::TETRA);
    }
    //fingerB
    {
      std::vector<Eigen::Matrix<double,3,1>> vss;
      std::vector<Eigen::Matrix<int,4,1>> tss;
      for(int i=offV3; i<offVB; i++)
        vss.push_back(solver.x().template segment<3>(i*3).template cast<double>());
      for(int i=offT3; i<offTB; i++)
        tss.push_back(solver.getTss()[i]._indices-Eigen::Matrix<int,4,1>::Constant(offV3));
      VTKWriter<double> os("Finger","Finger3D/ball"+std::to_string(outputIter)+".vtk",true);
      os.appendPoints(vss.begin(),vss.end());
      os.appendCells(tss.begin(),tss.end(),VTKWriter<double>::TETRA);
    }
    //cable
    {
      int off=0;
      VTKWriter<double> os("Base","Finger3D/cable"+std::to_string(outputIter)+".vtk",true);
      std::vector<Eigen::Matrix<double,3,1>> vss;
      std::vector<Eigen::Matrix<int,2,1>> iss;
      Eigen::Matrix<double,3,1> pos;
      pos.setZero();
      for(auto l1:solver.getL1ss()) {
        vss.push_back(solver.x().template segment<3>(l1.first[0]*3).template cast<double>());
        vss.push_back(solver.x().template segment<3>(l1.first[1]*3).template cast<double>());
        iss.push_back(Eigen::Matrix<int,2,1>(off,off+1));
        off+=2;
      }
      os.appendPoints(vss.begin(),vss.end());
      os.appendCells(iss.begin(),iss.end(),VTKWriter<double>::LINE);
    }
    outputIter++;
    updateDeformable(s,solver);
    drawer.nextFrame();
    return sim;
  });
  solver.solve(param);
  return 0;
}
