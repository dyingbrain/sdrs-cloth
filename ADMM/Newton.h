#ifndef NEWTON_H
#define NEWTON_H

#include "GradientDescend.h"

namespace PHYSICSMOTION {
template <int N>
class Newton : public GradientDescend<N> {
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
  Eigen::SparseLU<SMatT> _invH;
};
}

#endif
