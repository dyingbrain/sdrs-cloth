#include "CollisionMatrix.h"
#include <Utils/DebugGradient.h>

namespace PHYSICSMOTION {
template <int N,int M>
CollisionMatrix<N,M>::CollisionMatrix() {
  setZero();
}
template <int N,int M>
void CollisionMatrix<N,M>::setZero() {
  _abPss.setZero();
  _n.setZero();
  _ad0=0;
}
template <int N,int M>
void CollisionMatrix<N,M>::setN(const VecNT& n) {
  _n=n;
  _nLen2=n.squaredNorm();
}
template <int N,int M>
void CollisionMatrix<N,M>::addIdentity(int m,T coef) {
  if(m<_abPss.cols())
    _abPss(0,m)+=coef;
  else _ad0+=coef;
}
template <int N,int M>
void CollisionMatrix<N,M>::addNNT(int m,T coef) {
  _abPss(1,m)+=coef;
  _ad0+=coef;
}
template <int N,int M>
void CollisionMatrix<N,M>::solve(T shift,RHS& G) const {
  VecNT tmp;
  T ad0=_ad0;
  //accumulate d0 system
  T& Gd0=G.data()[_abPss.cols()*N];
  for(int i=0; i<_abPss.cols(); i++) {
    //accumulate LHS
    solveBlock(_abPss.col(i),shift,tmp=_n);
    ad0-=tmp.dot(_n)*(_abPss(1,i)*_abPss(1,i));
    //accumulate RHS
    VecNTCM Gp(G.data()+i*N);
    solveBlock(_abPss.col(i),shift,tmp=Gp);
    Gd0-=tmp.dot(_n)*_abPss(1,i);
  }
  //solve d0 system
  Gd0/=(ad0+shift);
  //solve other system
  for(int i=0; i<_abPss.cols(); i++) {
    VecNTM Gp(G.data()+i*N);
    Gp-=_n*(Gd0*_abPss(1,i));
    solveBlock(_abPss.col(i),shift,tmp=Gp);
    Gp=tmp;
  }
}
template <int N,int M>
void CollisionMatrix<N,M>::solveDense(T shift,RHS& G) const {
  MatT m=toDense();
  m.diagonal().array()+=shift;
  G=m.template cast<double>().ldlt().solve(G.template cast<double>()).template cast<T>();
}
template <int N,int M>
typename CollisionMatrix<N,M>::T CollisionMatrix<N,M>::maxCoeff() const {
  return std::max<T>(_ad0,(_abPss.row(0)+_abPss.row(1)).maxCoeff());
}
template <int N,int M>
typename CollisionMatrix<N,M>::MatT CollisionMatrix<N,M>::toDense() const {
  MatT ret;
  ret.setZero(N*_abPss.cols()+1,N*_abPss.cols()+1);
  for(int i=0; i<_abPss.cols(); i++) {
    ret.template block<N,N>(i*N,i*N)=MatNT::Identity()*_abPss(0,i)+_n*_n.transpose()*_abPss(1,i);
    ret.template block<1,N>(_abPss.cols()*N,i*N)=_n.transpose()*_abPss(1,i);
    ret.template block<N,1>(i*N,_abPss.cols()*N)=_n*_abPss(1,i);
  }
  ret(_abPss.cols()*N,_abPss.cols()*N)=_ad0;
  return ret;
}
template <int N,int M>
void CollisionMatrix<N,M>::debug() {
  CollisionMatrix<N,M> cm;
  cm.setN(VecNT::Random());
  for(int i=0; i<=M; i++)
    cm.addIdentity(i,rand()/(T)RAND_MAX);
  for(int i=0; i<M; i++)
    cm.addNNT(i,rand()/(T)RAND_MAX);

  RHS G,GRef;
  G.setRandom();
  GRef=G;

  DEFINE_NUMERIC_DELTA_T(T)
  T shift=rand()/(T)RAND_MAX;
  cm.solve(shift,G);
  cm.solveDense(shift,GRef);
  DEBUG_GRADIENT("GInv",G.norm(),(G-GRef).norm())
}
//helper
template <int N,int M>
void CollisionMatrix<N,M>::solveBlock(Vec2T ab,T shift,VecNT& b) const {
  ab[0]+=shift;
  T InvA=1/ab[0],BByA=ab[1]*InvA;
  T InvB=-BByA*InvA/(1+_nLen2*BByA);
  b=InvA*b+_n*(_n.dot(b)*InvB);
}
//instance 2D
template class CollisionMatrix<2,1>;
template class CollisionMatrix<2,2>;
template class CollisionMatrix<2,4>;
template class CollisionMatrix<2,6>;
template class CollisionMatrix<2,12>;
//instance 3D
template class CollisionMatrix<3,1>;
template class CollisionMatrix<3,2>;
template class CollisionMatrix<3,3>;
template class CollisionMatrix<3,6>;
template class CollisionMatrix<3,12>;
}
