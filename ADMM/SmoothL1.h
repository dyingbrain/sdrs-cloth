#ifndef SMOOTH_L1_H
#define SMOOTH_L1_H

#include "OptimizerTerm.h"

namespace PHYSICSMOTION {
//beta/2*||d-Ax||^2+Lambda^T(Ax-d)+
//betaY/2*||d-Ax||^2+
//k/2*\sqrt{||d-n*l||^2+\eps}
template <int N=3>
class SmoothL1 : public OptimizerTerm {
 public:
  typedef Eigen::Matrix<T,N,1> VecNT;
  typedef Eigen::Matrix<T,N,N> MatNT;
  typedef Eigen::Matrix<T,N,-1> MatNXT;
  void addL1(std::array<int,2> id,T k,T eps=1e-3f);
  int n() const override;
  std::shared_ptr<OptimizerTerm> copy() const override;
  T evalG(bool calcG,bool initL,SMatT* H,int y0Off) override;
  T evalGDirect(bool calcG,SMatT* H,int y0Off,bool projPSD) override;
  bool updateY(T betaY,T beta,T tolG) override;
  void reset(int mask) override;
  void debugEnergy(int M);
 private:
  //helper
  bool energyY(const VecNT& y,T& E,VecNT* G,MatNT* H,int i) const;
  //data
  std::vector<T> _alpha;
  //param
  std::vector<T> _k,_eps;
};
}

#endif
