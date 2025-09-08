#ifndef ADMM_UAV_TRAJECTORY_H
#define ADMM_UAV_TRAJECTORY_H

#include <Environment/MeshExact.h>
#include <SIPCollision/ThetaTrajectory.h>

namespace PHYSICSMOTION {
template <int N>
class ADMMUAVTrajectory {
 public:
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  typedef Eigen::Matrix<T,N,1> VecNT;
  typedef Eigen::Matrix<T,N,N> MatNT;
  typedef Eigen::Matrix<T,N,-1> MatNXT;
  typedef Eigen::Matrix<T,-1,N> MatXNT;
  typedef Eigen::Triplet<T,int> STrip;
  typedef ParallelVector<STrip> STrips;
  typedef Eigen::SparseMatrix<T,0,int> SMatT;
 public:
  //construction
  void setNodes(const std::vector<VecNT>& waypoints,int order,int numSegment,bool debug=false);
  bool setNodesExhaustive(const std::vector<VecNT>& waypoints,const MeshExact& obs,T margin,int order,int numSegmentBeg,int numSegmentEnd);
  bool intersect(const MeshExact& obs,T margin) const;
  void setNodes(int order,int numSegment);
  void buildSmoothness(int d);
  //getter
  Vec getCP() const;
  int numCP() const;
  SMatT getNodeToCP() const;
  SMatT getCPCons() const;
  SMatT getSmoothH() const;
  SMatT subdivide(int RES) const;
  const ThetaTrajectory<T>& getTraj() const;
  //writer
  void writeTrajVTK(const std::string& path,const Vec* pssIn,int RES) const;
  void writeNodesVTK(const std::string& path) const;
 private:
  bool intersect(int segId,const MeshExact& obs,T margin) const;
  //data
  MatXNT _nodes,_CPs;
  ThetaTrajectory<T> _traj;
  SMatT _nodeToCP,_CPCons,_smoothH;
};
}

#endif
