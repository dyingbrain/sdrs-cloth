#ifndef DIRECT_NEWTON_H
#define DIRECT_NEWTON_H

#include "Newton.h"

namespace PHYSICSMOTION {
template <int N>
class DirectNewton : public Newton<N> {
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
  T evalGD(const Vec& x,Vec* G,SMatT* H=NULL,bool projPSD=true);
  void debugGradient(const Vec& x);
};
}

#endif
