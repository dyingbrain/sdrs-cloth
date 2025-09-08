#ifndef COLLISION_MATRIX_H
#define COLLISION_MATRIX_H

#include <Eigen/Dense>
#include <Utils/Epsilon.h>
#include <Utils/Pragma.h>

namespace PHYSICSMOTION {
//each block is a*I+n*n^T*b
template <int N,int M>
class CollisionMatrix {
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  typedef Eigen::Matrix<T,N,1> VecNT;
  typedef Eigen::Matrix<T,N,N> MatNT;
  typedef Eigen::Map<VecNT> VecNTM;
  typedef Eigen::Map<const VecNT> VecNTCM;
  typedef Eigen::Matrix<T,N*M+1,1> RHS;
  typedef Eigen::Matrix<T,2,M> Mat2XMT;
 public:
  CollisionMatrix();
  void setZero();
  void setN(const VecNT& n);
  void addIdentity(int m,T coef);
  void addNNT(int m,T coef);
  void solve(T shift,RHS& G) const;
  void solveDense(T shift,RHS& G) const;
  T maxCoeff() const;
  MatT toDense() const;
  static void debug();
 private:
  void solveBlock(Vec2T ab,T shift,VecNT& b) const;
  //data
  Mat2XMT _abPss;
  T _nLen2,_ad0;
  VecNT _n;
};
}

#endif
