#ifndef COLLISION_SELF_H
#define COLLISION_SELF_H

#include "OptimizerTerm.h"
#include "CollisionMatrix.h"

namespace PHYSICSMOTION {
template <int N,int M>
class CollisionDetector;
template <int N=3,int M=3>
class CollisionSelf : public CLogx, public OptimizerTerm {
 public:
  typedef CLogx Penalty;
  typedef unsigned long long ID;
  typedef Eigen::Matrix<T,N,1> VecNT;
  typedef Eigen::Matrix<T,N+1,1> VecNdT;
  typedef Eigen::Matrix<T,N,N> MatNT;
  typedef Eigen::Matrix<T,N,N+1> MatNNdT;
  typedef Eigen::Matrix<T,N+1,N+1> MatNdT;
  typedef Eigen::Matrix<T,N,-1> MatNXT;
  typedef Eigen::Matrix<T,N*M*2,N*M*2> MatNMMT;
  typedef Eigen::Matrix<T,N*M*2,1> VecNMMT;
  typedef Eigen::Matrix<T,N*M*2+1,1> VecNMMdT;
  CollisionSelf(T r=0,T x0=1e-2f,T coef=1e-2f);
  const std::unordered_map<ID,int>& terms() const;
  void insertCollisions(const CollisionDetector<N,M>& detector);
  int removeCollisions(const std::unordered_set<int>& deleteHash);
  int removeCollisionsByDistance(T margin);
  int removeCollisionsByEnergy(T thres);
  VecM y0() override;
  VecCM y0() const override;
  VecCM G0() const override;
  int n() const override;
  std::shared_ptr<OptimizerTerm> copy() const override;
  T evalG(bool calcG,bool initL,SMatT* H,int y0Off) override;
  T evalGDirect(bool calcG,SMatT* H,int y0Off,bool projPSD) override;
  bool updateY(T betaY,T beta,T tolG) override;
  bool updateZ(T tolG) override;
  void reset(int mask) override;
  void save(int id,int mask) override;
  void load(int id,int mask) override;
  void debugEnergy(int m,T tolG);
  T eps() const;
 private:
  //helper
  void initializePlane(int i);
  bool energyYd(const VecNMMdT& yd,T& E,VecNMMdT* G,CollisionMatrix<N,M*2>* H,int i) const;
  bool energyYDirect(const VecNMMT& y,T& E,VecNMMT* G,MatNMMT* H,int i,bool projPSD) const;
  bool energyZ(const VecNT& z,T& E,VecNT* G,MatNT* H,int i) const;
  bool energyZd(const VecNdT& z,T& E,VecNdT* G,MatNdT* H,int i) const;
  //data
  Vec _d0,_Gd0;
  MatNXT _z;
  std::vector<T> _alphaY,_alphaZ;
  //param
  std::unordered_map<ID,int> _terms;
  T _r,_coef;
};
}

#endif
