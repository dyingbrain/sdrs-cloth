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

typedef std::unordered_map<int,std::pair<T,Deformable<3>::VecNT>> FIX;
FIX move(const FIX& fix,double t,typename Deformable<3>::VecNT dir,double rad=4) {
  FIX ret=fix;
  for (auto& p:ret) {
    p.second.second+=std::sin(t)*rad*dir;
  }
  return ret;
}
void createSkirt(int N1,int N2,double r0,double r1,MeshExact& m) {
#define ID(I,J) ((I)%N1)+(J)*N1
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  //vss
  double dr=2*M_PI*r1/N2;
  for(double r=r0;r<=r1+dr;r+=dr) {
    for(int d=0;d<N1;d++) {
      double angle=2*M_PI*(d+0.5)/N1;
      vss.push_back(Eigen::Matrix<double,3,1>(std::cos(angle),std::sin(angle),0)*r);
    }
  }
  //iss
  int NRing=(int)vss.size()/N1;
  for(int d=0;d<N1;d++) {
    for(int d2=0;d2<NRing-1;d2++) {
      iss.push_back(Eigen::Matrix<int,3,1>(ID(d,d2),ID(d,d2+1),ID(d+1,d2+1)));
      iss.push_back(Eigen::Matrix<int,3,1>(ID(d,d2),ID(d+1,d2+1),ID(d+1,d2)));
    }
  }
  m.init<double>(vss,iss);
#undef ID
}
void writeObj(const std::filesystem::path& path,const MeshExact& m) {
  std::ofstream os(path);
  for(auto v:m.vss())
    os << "v " << v[0] << " " << v[1] << " " << v[2] << std::endl;
  for(auto i:m.iss())
    os << "f " << i[0]+1 << " " << i[1]+1 << " " << i[2]+1 << std::endl;
}
int main(int argc,char** argv) {
  //build mesh grid
  MeshExact m;
  int N1=64,N2=32;
  createSkirt(N1,N2,1,10,m);
  //solver
  bool sim=false;
  Deformable<3> solver;
  solver.setK(1e4);
  solver.setB(1e2);
  solver.setCL(0.2);
  solver.setCH(1.5);
  solver.setMargin(0.005);
  solver.setG(Vec3T(0,0,-10));
  solver.setMassSpring(m);
  solver.setDt(0.05);
  solver.fix([](const Deformable<3>::VecNT& v){
    return v.norm()<=1.1;
  },1e4);

  //fix
  std::unordered_map<int,std::pair<T,Deformable<3>::VecNT>> fix=solver.getFix();

  //draw
  {
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
    drawer.mainLoop();
  }


  //Avoid override
  if(exists("MassSpring3DSkirt")) {
    std::cout << "Path already exists." << std::endl;
    return 0;
  }
  //Run the simulation
  int dirId=0;
  double lastSwitchTime=0,speed=1;
  std::vector<typename Deformable<3>::VecNT> dirs;
  dirs.push_back(Deformable<3>::VecNT({1,0,0}).normalized());
  dirs.push_back(Deformable<3>::VecNT({0,1,0}).normalized());
  dirs.push_back(Deformable<3>::VecNT({0,0,1}).normalized());
  dirs.push_back(Deformable<3>::VecNT({1,0,1}).normalized());
  dirs.push_back(Deformable<3>::VecNT({0,1,1}).normalized());
  dirs.push_back(Deformable<3>::VecNT({1,1,0}).normalized());
  recreate("MassSpring3DSkirt");
  double time=0;
  int outputIter=0;
  while(outputIter<100000) {
    OptimizerParam param;
    param._initBeta=1e2f;
    param._tolG=1e-2f;
    param._maxIter=1e4;
    param._printI=1;
    param._type=OptimizerParam::DIRECT_NEWTON;
    solver.getFix()=move(fix,time*speed,dirs[dirId]);
    solver.solve(param);
    {
      for(int i=0; i<(int)m.vss().size(); i++)
        m.vssNonConst()[i]=solver.x().template segment<3>(i*3).template cast<MeshExact::T>();
      VTKWriter<double> os("MassSpring3DSkirt","MassSpring3DSkirt/frame"+std::to_string(outputIter)+".vtk",true);
      writeObj("MassSpring3DSkirt/frame"+std::to_string(outputIter)+".obj",m);
      m.writeVTK(os,MeshExact::Mat3X4T::Identity());
    }
    if(solver.dt()==0)
      sim=false;
    time+=solver.dt();
    outputIter++;
    if((time-lastSwitchTime)*speed>M_PI*2) {
      lastSwitchTime=time;
      dirId=(dirId+1)%dirs.size();
    }
  }
  return 0;
}
