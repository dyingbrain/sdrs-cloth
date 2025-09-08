#ifndef COLLISION_DETECTOR_H
#define COLLISION_DETECTOR_H

#include "Optimizer.h"
#include <Environment/MeshExact.h>
#include <Utils/ParallelVector.h>
#include <deque>

namespace PHYSICSMOTION {
template <int N,int M>
class CollisionDetector : public CollisionChecker<N> {
 public:
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  typedef Eigen::Matrix<T,N,1> VecNT;
  typedef unsigned long long ID;
  typedef std::deque<ID> STACK;
  static const ID INVALID_ID=-1;
  CollisionDetector(int nCP,const MeshExact& obs);
  CollisionDetector(const MeshExact& mesh,const MeshExact& obs);
  void generateBroadBF(const Vec& x,const Vec& x2,T eps,T epsObs,bool noNeighbor,bool self,bool obs);
  void generateBroad(const Vec& x,const Vec& x2,T eps,T epsObs,bool noNeighbor,bool parallel,bool self,bool obs);
  Eigen::Matrix<int,2,1> generate(const Vec& x,const Vec& x2,Optimizer<N>& opt,bool noNeighbhor,bool parallel,int restoreSlot) override;
  Eigen::Matrix<int,2,1> remove(const Vec& x,Optimizer<N>& opt,T margin) override;
  std::string info(const Optimizer<N>& opt) const override;
  //getter
  void extractPos(const MeshExact& m,Vec& x) const;
  const std::vector<ID>& selfCollisions() const;
  const std::vector<ID>& obsCollisions() const;
  Eigen::Matrix<int,M*2,1> selfCollisionId(ID c) const;
  Eigen::Matrix<int,M,1> obsCollisionId(ID c) const;
  VecNT obsCollisionPos(ID c,int j) const;
  //ID
  static void separate(ID pair,int& l,int& r);
  static ID combine(ID l,ID r);
 private:
  void generateSelf(ID nodePair,const Vec& x,bool noNeighbor,bool parallel);
  void generateObs(ID nodePair,const Vec& x,T epsObs,bool parallel);
  void updateBVH(const Vec& x,const Vec& x2,T eps,bool parallel);
  void updateBVH(const Vec& x,const Vec& x2,T eps,int i);
  //mesh
  std::vector<int> _layerOffsets;
  std::vector<Node<int,BBoxExact>> _meshBVH;
  std::vector<Eigen::Matrix<int,M,1>> _meshIss;
  //obstacle
  MeshExact _obs;
  //pairs
  std::vector<ID> _nodePairs;
  ParallelVector<ID> _pssSelf,_pssObs;
  std::vector<STACK> _sss;
  int _nProc;
};
}

#endif
