#include "SkirtDatasetCommon.h"

int main(int argc,char** argv) {
  //build mesh grid
  MeshExact m;
  int N1=64,N2=32;
  double speed=0.5;
  createSkirt(N1,N2,1,10,m);
  //solver
  bool sim=false;
  Deformable<3> solver;
  solver.setK(1e4);
  solver.setB(1e3);
  solver.setCL(0.2);
  solver.setCH(1.5);
  solver.setMargin(0.005);
  solver.setDamping(0.01);
  solver.setG(Vec3T(0,0,-10));
  solver.setMassSpring(m);
  solver.setCollCoef(0);
  solver.setDt(0.01);
  solver.fix([](const Deformable<3>::VecNT& v){
    return v.norm()<=1.1;
  },1e4);

  //fix
  std::unordered_map<int,std::pair<T,Deformable<3>::VecNT>> fix=solver.getFix();
  
  //param
  OptimizerParam param;
  param._initBeta=1e2f;
  param._tolG=1e-2f;
  param._maxIter=1e4;
  param._printI=1;
  param._type=OptimizerParam::DIRECT_NEWTON;
  
  //visualization
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
  std::string datasetName="MassSpring3DSkirtNoCollision";
  if(!exists(datasetName))
    create(datasetName);
  int dirId=0,outputIter=0,lastSwitchIter=0,frame=0;
  loadCheckpoint(solver,dirId,outputIter,lastSwitchIter,frame,datasetName);
  Deformable<3>::VecNT dir({1,0,0});
  while(true) {
    solver.getFix()=move(fix,outputIter*solver.dt()*speed,dir);
    solver.solve(param);
    //output
    if((outputIter%1)==0) {
      saveCheckpoint(solver,dirId,outputIter,lastSwitchIter,frame,datasetName);
      for(int i=0; i<(int)m.vss().size(); i++)
        m.vssNonConst()[i]=solver.x().template segment<3>(i*3).template cast<MeshExact::T>();
      VTKWriter<double> os(datasetName,datasetName+"/frame"+std::to_string(frame)+".vtk",true);
      writeObj(datasetName+"/frame"+std::to_string(frame)+".obj",m);
      m.writeVTK(os,MeshExact::Mat3X4T::Identity());
      frame++;
    }
    //update
    if(solver.dt()==0)
      sim=false;
    outputIter++;
  }
  return 0;
}
