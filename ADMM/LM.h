#ifndef LM_H
#define LM_H

#include "Newton.h"

namespace PHYSICSMOTION {
template <int N>
class LM : public GradientDescend<N> {
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
  T evalGD(const Vec& x,Vec* G,SMatT* H = NULL);
  void solveSchur(Vec& x2,const Vec& x,const Vec& G,const Vec& H,T alpha);
  void debugGradient(const Vec& x);
};
}

#endif
