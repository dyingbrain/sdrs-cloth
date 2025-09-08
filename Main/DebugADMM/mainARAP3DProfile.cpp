#include <ADMM/Deformable.h>
#include <ADMM/DeformableVisualizer.h>
#include <Environment/EnvironmentVisualizer.h>
#include <Environment/EnvironmentUtils.h>
#include <Utils/VTKWriter.h>
#include <Utils/Utils.h>

using namespace PHYSICSMOTION;
using namespace DRAWER;
typedef FLOAT T;
DECL_MAT_VEC_MAP_TYPES_T

void solve(FLOAT sz) {
  //build mesh box
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  addBox(vss,iss,Eigen::Matrix<double,3,1>(0,0,0),Eigen::Matrix<double,3,1>(10,3,3));
  MeshExact m(vss,iss,true);

  //solver
  Deformable<3> solver;
  solver.setARAP3D(m,91,sz);
  solver.setDt(0);
  solver.fix([&](const Eigen::Matrix<T,3,1>& v) {
    return v[0]<0.1f;
  },0);
  solver.setK(1e3);
  solver.setCollCoef(0);
  solver.setG(Vec3T(0,0,-50));
  m=solver.buildSelfCollisionMesh();
  std::cout << "Generated " << solver.getTss().size() << " tets" << std::endl;

  OptimizerParam param;
  param._initBeta=1e2f;
  param._tolG=1e-2f;
  param._maxIter=1e5;
  param._debugGradientI=50;
  param._type=OptimizerParam::ADMM;
  solver.solve(param);
}
int main(int argc, char** argv) {
  for (FLOAT sz : {1.0f, 0.75f, 0.5f, 0.25f, 0.125f}) {
    std::string name = "ARAP3DProfile"+std::to_string(sz) + ".txt";
    if (exists(name))
      continue;
    if (freopen(name.c_str(), "w", stdout) == NULL) {
      // Handle error if redirection fails
      perror("Failed to redirect stdout");
      return 1;
    }
    std::cout << "Begin-testing " << sz << " size" << std::endl;
    solve(sz);
    std::cout << "End-testing " << sz << " size" << std::endl;
  }
  return 0;
}