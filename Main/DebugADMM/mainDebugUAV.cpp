#include <ADMM/ADMMUAVTrajectory.h>

using namespace PHYSICSMOTION;
typedef FLOAT T;
DECL_MAT_VEC_MAP_TYPES_T

template <int N,int RES=32>
void debugUAV() {
  T rand=5;
  ADMMUAVTrajectory<N> traj;
  std::vector<typename ADMMUAVTrajectory<N>::VecNT> waypoints;
  for(int i=0; i<5; i++)
    waypoints.push_back(ADMMUAVTrajectory<N>::VecNT::Random()*rand);
  traj.setNodes(waypoints,5,10,true);
  traj.writeTrajVTK("path"+std::to_string(N)+"D.vtk",NULL,RES);
}
int main(int argc,char** argv) {
  debugUAV<2>();
  debugUAV<3>();
  return 0;
}
