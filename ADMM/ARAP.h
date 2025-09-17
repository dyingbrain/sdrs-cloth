#ifndef ARAP_H
#define ARAP_H

#include "OptimizerTerm.h"

namespace PHYSICSMOTION {
//\sum_i\sum_{j\neq i}-log(n^T(d_j-d_0))+k/2*||F*invF0-R||_F^2
template <int N=3>
class ARAP : public CLogx, public OptimizerTerm {
 public:
  typedef CLogx Penalty;
  typedef Eigen::Matrix<T,N,1> VecVT;
  typedef Eigen::Matrix<T,N,N> MatVT;
  typedef Eigen::Matrix<T,N,-1> MatVXT;
  typedef Eigen::Matrix<T,N*(N+1),1> VecYT;
  typedef Eigen::Matrix<T,N*(N+1),-1> MatYXT;
  typedef Eigen::Matrix<T,N*(N+2),1> VecYDT;
  typedef Eigen::Matrix<T,N*(N+2),N*(N+2)> MatYDT;
  typedef Eigen::Matrix<T,N*N,-1> MatRXT;
  void addElement(std::array<int,N+1> id,T k,const MatVT& F0);
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
  void debugEnergy(int M,T tolG);
 private:
  //helper
  bool energyYD(const VecYDT& yd,T& E,VecYDT* G,MatYDT* H,int i) const;
  bool energyYDDirect(const VecYDT& yd,T& E,VecYDT* G,MatYDT* H,int i,bool projPSD) const;
  bool energyZ(const VecVT& z,T& E,VecVT* G,MatVT* H,int tetId,int fid) const;
  VecVT center(int tetId) const;
  VecVT normal(int tetId,int fid) const;
  T minEdge(const MatVT& F0) const;
  //data
  MatRXT _R;
  MatVXT _d0,_z[N+1],_Gd0;
  std::vector<T> _alphaY,_alphaZ[N+1];
  //param
  std::vector<MatVT> _invF0;
  std::vector<T> _k;
  bool _updateR=true;
  bool _debug=false;
};
}

#endif
