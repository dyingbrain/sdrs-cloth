#ifndef GRADIENT_DESCEND_H
#define GRADIENT_DESCEND_H

#include "Optimizer.h"

namespace PHYSICSMOTION {
template <int N>
class GradientDescend : public Optimizer<N> {
 public:
  using typename Optimizer<N>::T;
  DECL_MAT_VEC_MAP_TYPES_T
  using typename Optimizer<N>::VecNT;
  using typename Optimizer<N>::MatNT;
  using typename Optimizer<N>::STrip;
  using typename Optimizer<N>::STrips;
  using typename Optimizer<N>::SMatT;
  using typename Optimizer<N>::LinearConstraint;
  void optimize(const OptimizerParam& param) override;
 protected:
  Vec assembleX();
  bool updateZ(const Vec& x,T tolG);
  T evalGD(const Vec& x,Vec* G,SMatT* H=NULL);
  void debugGradient(const Vec& x);
  void lineSearch(Vec& x,T& E,const Vec& D,T& alpha,const OptimizerParam& param);
};
}

#endif
