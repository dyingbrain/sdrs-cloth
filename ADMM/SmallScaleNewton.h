#ifndef SMALL_SCALE_NEWTON_H
#define SMALL_SCALE_NEWTON_H

#include "CollisionMatrix.h"
#include <Eigen/Dense>
#include <iostream>

namespace PHYSICSMOTION {
template <int N,typename MAT_TYPE>
struct SmallScaleNewton {
  typedef FLOAT T;
  typedef Eigen::Matrix<T,N,1> VecNT;
  typedef Eigen::Matrix<T,N,N> MatNT;
  typedef std::function<bool(const VecNT&,T&,VecNT*,MAT_TYPE*)> Energy;
  template <typename BARRIER>
  static bool energySoft(const BARRIER& b,const VecNT& x,T& E,VecNT* G,MAT_TYPE* H) {
    //z is in the unit sphere
    T lenX=x.norm(),D=0,DD=0;
    VecNT pos=x/lenX;
    E+=b.template eval<FLOAT>(1+b._x0-lenX,G?&D:NULL,H?&DD:NULL,0,1);
    if(!isfinite(E))
      return false;
    pos=x/lenX;
    if(G)
      *G-=D*pos;
    if(H)
      *H+=-D*(MatNT::Identity()-pos*pos.transpose())/lenX+DD*pos*pos.transpose();
    return true;
  }
  static bool optimize(T& alpha,VecNT& x,T& E,VecNT& G,MAT_TYPE& H,Energy energy,T tolG,int maxIter=1) {
    int iter=0;
    MAT_TYPE H2;
    VecNT D,G2,x2;
    T E2,alphaDec=0.5,alphaInc=3.0,alphaMax=1e10;
    alpha=std::max<T>(alpha,Epsilon<T>::defaultEps());
    if(alpha>=alphaMax)
      alpha=1;//try again
    //assemble
    if(!energy(x,E,&G,&H))
      return false;
    //main loop
    while(alpha<alphaMax && iter<maxIter) {
      //search direction
      solve(H,alpha,D=G);
      //termination
      if(G.cwiseAbs().maxCoeff()<=tolG)
        return true;
      //test
      x2=x-D;
      if(energy(x2,E2,&G2,&H2) && E2<E) {
        alpha=std::max<T>(alpha*alphaDec,Epsilon<T>::defaultEps());
        x=x2;
        E=E2;
        G=G2;
        H=H2;
        iter++;
      } else {
        alpha*=alphaInc;
      }
    }
    return alpha<alphaMax;
  }
  static bool optimizeOnSphere(T& alpha,VecNT& x,T& E,VecNT& G,MAT_TYPE& H,Energy energy,T tolG,int maxIter=1) {
    int iter=0;
    MAT_TYPE H2;
    VecNT invHN,D,G2,GPrj,x2;
    T E2,alphaDec=0.5,alphaInc=3.0,alphaMax=1e10;
    alpha=std::max<T>(alpha,Epsilon<T>::defaultEps());
    if(alpha>=alphaMax)
      alpha=1;//try again
    //assemble
    x.normalize();//project
    if(!energy(x,E,&G,&H))
      return false;
    //main loop
    while(alpha<alphaMax && iter<maxIter) {
      //search direction
      solve(H,alpha,invHN=x);
      solve(H,alpha,D=G);
      D-=invHN*invHN.dot(G)/invHN.dot(x);
      //termination
      GPrj=G-G.dot(x)*x;
      if(GPrj.cwiseAbs().maxCoeff()<=tolG)
        return true;
      //test
      x2=x-D;
      x2.normalize();
      if(energy(x2,E2,&G2,&H2) && E2<E) {
        alpha=std::max<T>(alpha*alphaDec,Epsilon<T>::defaultEps());
        x=x2;
        E=E2;
        G=G2;
        H=H2;
        iter++;
      } else {
        alpha*=alphaInc;
      }
    }
    return alpha<alphaMax;
  }
  template <int A,int B>
  static void solve(const CollisionMatrix<A,B>& H,T shift,VecNT& x) {
    H.solve(shift,x);
  }
  static void solve(const MatNT& H,T shift,VecNT& x) {
    x=(H+MatNT::Identity()*shift).template cast<double>().ldlt().solve(x.template cast<double>()).template cast<T>();
  }
  template <int A,int B>
  static T maxCoeff(const CollisionMatrix<A,B>& H) {
    return H.maxCoeff();
  }
  static T maxCoeff(const MatNT& H) {
    return H.cwiseAbs().maxCoeff();
  }
};
}

#endif
