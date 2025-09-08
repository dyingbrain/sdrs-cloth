#ifndef ADMM_H
#define ADMM_H

#include "Optimizer.h"

namespace PHYSICSMOTION {
//beta/2*||Ax-y||^2+Lambda^T(Ax-y)+betaX/2*||x-xLast||^2+
//1/2*x^TH*x+g
template <int N>
class ADMM : public Optimizer<N> {
 public:
  using typename Optimizer<N>::T;
  DECL_MAT_VEC_MAP_TYPES_T
  using typename Optimizer<N>::VecNT;
  using typename Optimizer<N>::MatNT;
  using typename Optimizer<N>::STrip;
  using typename Optimizer<N>::STrips;
  using typename Optimizer<N>::SMatT;
  using typename Optimizer<N>::LinearConstraint;
  ADMM();
  void assemble(int nX,const LinearConstraint& cons= {}) override;
  void optimize(const OptimizerParam& param) override;
 private:
  void assembleLHS(int nX);
  T eval(Vec* G);
  void updateX();
  bool updateY(T tolG);
  void updateL();
  //data
  Vec _xLast,_RHS;
  T _beta,_betaX,_betaY;
  //matrix A
  SMatT _A,_LHS;
  Eigen::SparseLU<SMatT> _LHSInv;
};
}

#endif
