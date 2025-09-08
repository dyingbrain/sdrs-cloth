#include <ADMM/Deformable.h>
#include <ADMM/DeformableVisualizer.h>
#include <Environment/EnvironmentUtils.h>
#include <Environment/EnvironmentVisualizer.h>
#include <Utils/VTKWriter.h>
#include <Utils/Utils.h>
#include "HSV2RGB.h"

#include <TinyVisualizer/Drawer.h>
#include <TinyVisualizer/Povray.h>
#include <TinyVisualizer/Camera2D.h>
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

void createAgentsCircle(Vec& x,Vec& xT,T r,int N,T rand=0.1f) {
  x.resize(N*2);
  xT.resize(N*2);
  //(M_PI*2*R/N)*0.5f=r;
  T R=2*r*N/(M_PI*2);
  for(int i=0; i<N; i++) {
    T angle=M_PI*2*i/N,angle2=angle+M_PI;
    x.template segment<2>(i*2)=Vec2T(cos(angle),sin(angle))*R;
    xT.template segment<2>(i*2)=Vec2T(cos(angle2),sin(angle2))*R;
  }
  x+=Vec::Random(x.size())*rand*r;
  xT+=Vec::Random(x.size())*rand*r;
}
void createCircle(MeshExact& m,T rObs,int N) {
  std::vector<Eigen::Matrix<double,2,1>> vss;
  std::vector<Eigen::Matrix<int,2,1>> iss;
  for(int i=0; i<N; i++) {
    T angle=M_PI*2*i/N;
    vss.push_back(Eigen::Matrix<double,2,1>(cos(angle),sin(angle))*rObs);
  }
  for(int i=0; i<N; i++)
    iss.push_back(Eigen::Matrix<int,2,1>(i,(i+1)%N));
  m=MeshExact(vss,iss);
}
int main(int argc,char** argv) {
  Vec x,xT;
  int N=300;
  T r=1,rObs=5,margin=r*0.1f;
  MeshExact obs;
  createAgentsCircle(x,xT,r,N);
  createCircle(obs,rObs,N);
  obs.translate(MeshExact::Vec3T(0,rObs/2,0));
  r*=0.5f;
  //solver
  bool sim=false;
  Deformable<2> solver;
  solver.setAgents(x,xT,r);
  solver.setObstacle(obs);
  solver.setDt(0);
  solver.setCollCoef(1e2f);
  solver.setMargin(margin);
  solver.setG(Vec2T(0,0));
  //Mesh
  std::vector<Eigen::Matrix<int,3,1>> iss;
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<std::vector<Eigen::Matrix<double,3,1>>> pss(N);
  addCircle(vss,iss,Eigen::Matrix<double,3,1>::Zero(),solver.actualR(),32);
  MeshExact m(vss,iss);

  //draw
  Drawer drawer(argc,argv);
  std::shared_ptr<CaptureGIFPlugin> ss(new CaptureGIFPlugin(GLFW_KEY_4,"screenshot.gif",drawer.FPS(),true));
  drawer.addPlugin(std::shared_ptr<Plugin>(new CameraExportPlugin(GLFW_KEY_2,GLFW_KEY_3,"camera.dat")));
  drawer.addPlugin(std::shared_ptr<Plugin>(new CaptureGIFPlugin(GLFW_KEY_1,"record.gif",drawer.FPS())));
  std::shared_ptr<CompositeShape> s=visualizeDeformable(solver,1);
  drawer.addShape(s);
  drawer.clearLight();
  drawer.addCamera2D(20);
  s->getChild(0)->setColorAmbient(GL_LINES,1,.5f,0);
  s->getChild(s->numChildren()-1)->setColorAmbient(GL_LINES,.5f,.5f,.5f);
  drawer.setKeyFunc([&](GLFWwindowPtr wnd,int key,int scan,int action,int mods,bool captured) {
    if(captured)
      return;
    else if(key==GLFW_KEY_R && action==GLFW_PRESS)
      sim=!sim;
  });
  drawer.mainLoop();

  recreate(std::filesystem::path("Navigation"));
  {
    VTKWriter<double> os("Base","Navigation/obs.vtk",true);
    obs.writeVTK(os,MeshExact::Mat3X4T::Identity());
  }
  int outputIter=0;
  OptimizerParam param;
  param._initAlpha=0.001;
  param._initBeta=1e3f;
  param._initBetaX=1e3f;
  param._tolG=0;
  param._maxIter=1e5;
  //param._debugGradientI=50;
  param._ensureMonotonic=false;
  param._type=OptimizerParam::ADMM;
  solver.setCB([&]() {
    //write agents
    {
      std::vector<Eigen::Matrix<double,3,1>> colors;
      MeshExact::Mat3X4T t=MeshExact::Mat3X4T::Identity();
      VTKWriter<double> os("Agent","Navigation/agent"+std::to_string(outputIter)+".vtk",true);
      for(int i=0; i<N; i++) {
        hsv in({360.*i/N,1,1});
        rgb out=hsv2rgb(in);
        for(int j=0; j<(int)m.vss().size(); j++)
          colors.push_back({out.r,out.g,out.b});
        t.col(3).template segment<2>(0)=solver.x().template segment<2>(i*2).template cast<MeshExact::T>();
        m.writeVTK(os,t);
      }
      os.appendCustomPointColorData("AgentColor",colors.begin(),colors.end());
    }
    //write path
    {
      std::vector<Eigen::Matrix<double,4,1>> colors;
      Eigen::Matrix<double,3,1> pos=Eigen::Matrix<double,3,1>::Zero();
      VTKWriter<double> os("Agent","Navigation/path"+std::to_string(outputIter)+".vtk",true);
      for(int i=0; i<N; i++) {
        pos.template segment<2>(0)=solver.x().template segment<2>(i*2).template cast<double>();
        pss[i].push_back(pos);
        //write
        if((int)pss[i].size()>1) {
          hsv in({360.*i/N,1,1});
          rgb out=hsv2rgb(in);
          os.setRelativeIndex();
          os.appendPoints(pss[i].begin(),pss[i].end());
          VTKWriter<double>::IteratorIndex<Eigen::Matrix<int,2,1>> beg(0,1,0),end((int)pss[i].size()-1,1,0);
          os.appendCells(beg,end,VTKWriter<double>::LINE,true);
          for(int j=0; j<(int)pss[i].size(); j++)
            colors.push_back({out.r,out.g,out.b,j/(double)pss[i].size()});
        }
      }
      if(!colors.empty())
        os.appendCustomPointColorData("PathColor",colors.begin(),colors.end());
    }
    outputIter++;
    updateDeformable(s,solver);
    drawer.nextFrame();
    return sim;
  });
  solver.solve(param);
  return 0;
}
