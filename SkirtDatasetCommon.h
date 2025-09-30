#pragma once

#include <ADMM/Deformable.h>
#include <ADMM/DeformableVisualizer.h>
#include <Environment/EnvironmentVisualizer.h>
#include <Utils/VTKWriter.h>
#include <Utils/Utils.h>
#include <regex>

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
int extractCheckpoint(const std::string& filename) {
  if(!endsWith(filename,".ckpt"))
    return -1;
  std::regex pattern("frame(\\d+)"); // Matches one or more digits
  std::smatch matches;
  if(std::regex_search(filename,matches,pattern))
    return std::stoi(matches[1].str());
  else return -1;
}
void loadCheckpoint(Deformable<3>& solver,int& dirId,int& outputIter,int& lastSwitchIter,int& frame,std::string datasetName) {
  int latestCkpt=-1;
  for(auto file:files(datasetName)) {
    int ckpt=extractCheckpoint(file.string());
    if(ckpt>=latestCkpt)
      latestCkpt=ckpt;
  }
  if(latestCkpt<0)
    return;
  std::cout << "Continuing from " << "frame"+std::to_string(latestCkpt)+".ckpt" << std::endl;
  std::ifstream is(datasetName+"/frame"+std::to_string(latestCkpt)+".ckpt",std::ios::binary);
  readBinaryData(dirId,is);
  readBinaryData(outputIter,is);
  readBinaryData(lastSwitchIter,is);
  readBinaryData(frame,is);
  readBinaryData(solver.x(),is);
  readBinaryData(solver.xL(),is);
  readBinaryData(solver.xLL(),is);
}
void saveCheckpoint(const Deformable<3>& solver,int dirId,int outputIter,int lastSwitchIter,int frame,std::string datasetName) {
  std::ofstream os(datasetName+"/frame"+std::to_string(frame)+".ckpt",std::ios::binary);
  writeBinaryData(dirId,os);
  writeBinaryData(outputIter,os);
  writeBinaryData(lastSwitchIter,os);
  writeBinaryData(frame,os);
  writeBinaryData(solver.x(),os);
  writeBinaryData(solver.xL(),os);
  writeBinaryData(solver.xLL(),os);
}
