#ifndef MASS_SPRING_H
#define MASS_SPRING_H

#include "OptimizerTerm.h"

namespace PHYSICSMOTION {
//beta/2*||d-Ax||^2+Lambda^T(Ax-d)+
//betaY/2*||d-Ax||^2+
//k/2*||d-n*l||^2-log(n^Td-lL)-log(lH-||d||)
template <int N=3>
class MassSpring : public CLogx, public OptimizerTerm {
 public:
  typedef CLogx Penalty;
  typedef Eigen::Matrix<T,N,1> VecNT;
  typedef Eigen::Matrix<T,N,N> MatNT;
  typedef Eigen::Matrix<T,N,-1> MatNXT;
  void addSpring(std::array<int,2> id,T k,T l,T lL,T lH);
  int n() const override;
  std::shared_ptr<OptimizerTerm> copy() const override;
  T evalG(bool calcG,bool initL,SMatT* H,int y0Off) override;
  bool updateY(T betaY,T beta,T tolG) override;
  bool updateZ(T tolG) override;
  void reset(int mask) override;
  void save(int id,int mask) override;
  void load(int id,int mask) override;
  void debugEnergy(int M,T tolG);
 private:
  //helper
  bool energyY(const VecNT& y,T& E,VecNT* G,MatNT* H,int i) const;
  //data
  MatNXT _z;
  std::vector<T> _alpha;
  //param
  std::vector<T> _k,_l,_lL,_lH;
};
}

#endif
