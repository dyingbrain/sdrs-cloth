#include <ADMM/Deformable.h>
#include <ADMM/DeformableVisualizer.h>
#include <Environment/EnvironmentUtils.h>
#include <Environment/EnvironmentVisualizer.h>
#include <Utils/VTKWriter.h>
#include <Utils/Utils.h>
#include "HSV2RGB.h"

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
void solve(int N) {
  Vec x,xT;
  T r=1,rObs=5,margin=r*0.1f;
  MeshExact obs;
  createAgentsCircle(x,xT,r,N);
  createCircle(obs,rObs,N);
  obs.translate(MeshExact::Vec3T(0,rObs/2,0));
  r*=0.5f;
  //solver
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

  OptimizerParam param;
  param._initAlpha=0.001;
  param._initBeta=1e3f;
  param._initBetaX=1e3f;
  param._tolG=0;
  param._maxIter=1e5;
  //param._debugGradientI=50;
  param._ensureMonotonic=false;
  param._type=OptimizerParam::ADMM;
  solver.solve(param);
}
int main(int argc, char** argv) {
  for (int n : {50, 100, 200, 300, 400}) {
    std::string name = "NavigationProfile" + std::to_string(n) + ".txt";
    if (exists(name))
      continue;
    if (freopen(name.c_str(), "w", stdout) == NULL) {
      // Handle error if redirection fails
      perror("Failed to redirect stdout");
      return 1;
    }
    std::cout << "Begin-testing " << n << " agents" << std::endl;
    solve(n);
    std::cout << "End-testing " << n << " agents" << std::endl;
  }
  return 0;
}
