#ifndef SPATIAL_HASH_H
#define SPATIAL_HASH_H

#include "CollisionDetector.h"

namespace PHYSICSMOTION {
template <int N>
class CollisionDetector<N,1> : public CollisionChecker<N> {
 public:
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  typedef Eigen::Matrix<T,N,1> VecNT;
  typedef unsigned long long ID;
  typedef std::deque<ID> STACK;
  static const ID INVALID_ID=-1;
  struct SpatialHashNode {
    BBoxExact _bb;
    VecNT _ctr;
    T _radius;
  };
  CollisionDetector(const MeshExact& obs);
  void generateBroadBF(const Vec& x,const Vec& x2,T eps,T epsObs,bool self,bool obs);
  void generateBroad(const Vec& x,const Vec& x2,T eps,T epsObs,bool self,bool obs);
  Eigen::Matrix<int,2,1> generate(const Vec& x,const Vec& x2,Optimizer<N>& opt,bool noNeighbhor,bool parallel,int restoreSlot) override;
  Eigen::Matrix<int,2,1> remove(const Vec& x,Optimizer<N>& opt,T margin) override;
  std::string info(const Optimizer<N>& opt) const override;
  //getter
  void extractPos(const MeshExact& m,Vec& x) const;
  const std::vector<ID>& selfCollisions() const;
  const std::vector<ID>& obsCollisions() const;
  Eigen::Matrix<int,2,1> selfCollisionId(ID c) const;
  Eigen::Matrix<int,1,1> obsCollisionId(ID c) const;
  VecNT obsCollisionPos(ID c,int j) const;
  //ID
  static void separate(ID pair,int& l,int& r);
  static ID combine(ID l,ID r);
 private:
  int hashOff(const VecNT& pt) const;
  Eigen::Matrix<int,N,1> hash(const VecNT& pt) const;
  BBoxExact getBB(const Vec& x,const Vec& x2,int id,T eps) const;
  SpatialHashNode getNode(const Vec& x,const Vec& x2,int id,T eps) const;
  void reduce(std::function<void(SpatialHashNode&,SpatialHashNode&)> op);
  void generateObs(int vid,const BBoxExact& bb);
  void buildHash(const Vec& x,const Vec& x2,T eps);
  //obstacle
  MeshExact _obs;
  //pairs
  ParallelVector<ID> _pssSelf,_pssObs;
  //spatial hash
  std::vector<int> _indices,_offsetsInv,_offsets,_starts,_ends;
  std::vector<SpatialHashNode> _nodes,_nodesBkg;
  Eigen::Matrix<int,N,1> _nrCell,_stride;
  BBoxExact _bb;
  T _R,_invR;
};
}

#endif
