#include <Environment/EnvironmentUtils.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  std::vector<Eigen::Matrix<int,4,1>> ess;
  addBox(vss,iss,Eigen::Matrix<double,3,1>(0,0,0),Eigen::Matrix<double,3,1>(10,3,3));
  generateTetMesh(vss,iss,"tet.mesh");
  readTetMesh(vss,ess,"tet.mesh");
  writeTetVTK(vss,ess,"tet.vtk");
  return 0;
}
