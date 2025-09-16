#include "SmoothL1.h"
#include "SmallScaleNewton.h"
#include <Utils/SparseUtils.h>
#include <Utils/DebugGradient.h>

namespace PHYSICSMOTION {
template <int N>
void SmoothL1<N>::addL1(std::array<int,2> id,T k,T eps) {
  addBlockId<T>(_E,n()*N,id[0]*N,N, 1);
  addBlockId<T>(_E,n()*N,id[1]*N,N,-1);
  _k.push_back(k);
  _eps.push_back(eps);
}
template <int N>
int SmoothL1<N>::n() const {
  return (int)_k.size();
}
template <int N>
std::shared_ptr<OptimizerTerm> SmoothL1<N>::copy() const {
  return std::shared_ptr<OptimizerTerm>(new SmoothL1<N>);
}
template <int N>
typename SmoothL1<N>::T SmoothL1<N>::evalG(bool calcG,bool initL,SMatT* H,int y0Off) {
  _evalgOnly=true;
  T g=0;
  if(calcG)
    _G.resize(N,n());
  if(H)
    _HBlks.resize(n());
  OMP_PARALLEL_FOR_
  for(int i=0; i<n(); i++) {
    T E;
    MatNT HBlk;
    VecNT y=_Ax.col(i),G;
    energyY(y,E,calcG?&G:NULL,H?&HBlk:NULL,i);
    parallelAdd(g,E);
    if(calcG)
      _G.col(i)=G;
    if(H) {
      _HBlks[i]._blk=HBlk;
      _HBlks[i]._nY=N;
    }
  }
  if(calcG && initL)
    initializeL();
  if(H)
    assembleHessian(*H,y0Off);
  return g;
}
template <int N>
typename SmoothL1<N>::T SmoothL1<N>::evalGDirect(bool calcG,SMatT* H,int y0Off,bool projPSD) {
  return evalG(calcG,false,H,y0Off);
}
template <int N>
bool SmoothL1<N>::updateY(T betaY,T beta,T tolG) {
  _evalgOnly=false;
  bool succ=true;
  _betaY=betaY;
  _beta=beta;
  _yLast=_y;
  if((int)_alpha.size()<n())
    _alpha.resize(n(),1);
  OMP_PARALLEL_FOR_
  for(int i=0; i<n(); i++) {
    T E;
    MatNT H;
    VecNT y=_y.col(i),G;
    SmallScaleNewton<N,MatNT> opt;
    auto energyFunc=[&](const VecNT& x,T& E,VecNT* G,MatNT* H)->bool {
      return energyY(x,E,G,H,i);
    };
    if(!opt.optimize(_alpha[i],y,E,G,H,energyFunc,tolG))
      succ=false;
    _y.col(i)=y;
  }
  return succ;
}
template <int N>
void SmoothL1<N>::reset(int mask) {
  OptimizerTerm::reset(mask);
  if(mask&MASK_ALPHA)
    _alpha.clear();
}
template <int N>
void SmoothL1<N>::debugEnergy(int M) {
  _y.setRandom(N,M);
  _yLast.setRandom(N,M);
  _L.setRandom(N,M);
  _Ax.setRandom(N,M);

  _k.resize(M);
  _eps.resize(M);

  for(int i=0; i<M; i++) {
    _k[i]=std::rand()/(T)RAND_MAX;
    _eps[i]=1e-2f*std::rand()/(T)RAND_MAX;
  }

  _beta=std::rand()/(T)RAND_MAX;
  _betaY=std::rand()/(T)RAND_MAX;

  DEFINE_NUMERIC_DELTA_T(T)
  for(bool evalgOnly: {
        false,true
      }) {
    std::cout << "DebugEnergyY evalgOnly=" << evalgOnly << std::endl;
    for(int i=0; i<M; i++) {
      T E,E2;
      MatNT H;
      VecNT y,y2,dy,G,G2;
      y=_y.col(i);
      dy.setRandom();
      _evalgOnly=evalgOnly;
      if(energyY(y,E,&G,&H,i)) {
        y2=y+dy*DELTA;
        energyY(y2,E2,&G2,NULL,i);
        DEBUG_GRADIENT("G",G.dot(dy),G.dot(dy)-(E2-E)/DELTA)
        DEBUG_GRADIENT("H",(H*dy).norm(),((G2-G)/DELTA-H*dy).norm())
      }
    }
  }
}
//helper
template <int N>
bool SmoothL1<N>::energyY(const VecNT& y,T& E,VecNT* G,MatNT* H,int i) const {
  //energy
  E=0;
  T D=0;
  VecNT DY;
  if(!_evalgOnly) {
    E+=_beta*(y-_Ax.col(i)).squaredNorm()/2;
    E+=_L.col(i).dot(_Ax.col(i)-y);
    E+=_betaY*(y-_yLast.col(i)).squaredNorm()/2;
  }
  D=sqrt(y.squaredNorm()+_eps[i]);
  E+=_k[i]*D;
  //gradient
  if(G) {
    G->setZero();
    if(!_evalgOnly) {
      *G+=_beta*(y-_Ax.col(i));
      *G-=_L.col(i);
      *G+=_betaY*(y-_yLast.col(i));
    }
    DY=y/D;
    *G+=_k[i]*DY;
  }
  //hessian
  if(H) {
    if(!_evalgOnly) {
      *H=MatNT::Identity();
      *H*=_beta+_betaY;
    } else H->setZero();
    *H+=(MatNT::Identity()-DY*DY.transpose())*_k[i]/D;
  }
  return true;
}
//instance
template class SmoothL1<2>;
template class SmoothL1<3>;
}
