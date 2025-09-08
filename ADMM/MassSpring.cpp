#include "MassSpring.h"
#include "SmallScaleNewton.h"
#include <Utils/SparseUtils.h>
#include <Utils/DebugGradient.h>

namespace PHYSICSMOTION {
template <int N>
void MassSpring<N>::addSpring(std::array<int,2> id,T k,T l,T lL,T lH) {
  addBlockId<T>(_E,n()*N,id[0]*N,N, 1);
  addBlockId<T>(_E,n()*N,id[1]*N,N,-1);
  _k.push_back(k);
  _l.push_back(l);
  _lL.push_back(lL);
  _lH.push_back(lH);
}
template <int N>
int MassSpring<N>::n() const {
  return (int)_k.size();
}
template <int N>
std::shared_ptr<OptimizerTerm> MassSpring<N>::copy() const {
  return std::shared_ptr<OptimizerTerm>(new MassSpring<N>);
}
template <int N>
typename MassSpring<N>::T MassSpring<N>::evalG(bool calcG,bool initL,SMatT* H,int y0Off) {
  if(initL) {
    //initialize activation range
    for(int i=0; i<n(); i++)
      if(i==0)
        _x0=(double)_lL[i];
      else _x0=(double)std::min<T>(_x0,_lL[i]);
  }
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
bool MassSpring<N>::updateY(T betaY,T beta,T tolG) {
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
bool MassSpring<N>::updateZ(T tolG) {
  if(_z.cols()!=n())
    _z.resize(N,n());
  OMP_PARALLEL_FOR_
  for(int i=0; i<n(); i++)
    _z.col(i)=_y.col(i).normalized();
  return true;
}
template <int N>
void MassSpring<N>::reset(int mask) {
  OptimizerTerm::reset(mask);
  if(mask&MASK_Z)
    _z.resize(N,0);
  if(mask&MASK_ALPHA)
    _alpha.clear();
}
template <int N>
void MassSpring<N>::save(int id,int mask) {
  OptimizerTerm::save(id,mask);
  auto ptr=std::dynamic_pointer_cast<MassSpring<N>>(_saved[id]);
  if(mask&MASK_Z)
    ptr->_z=_z;
}
template <int N>
void MassSpring<N>::load(int id,int mask) {
  OptimizerTerm::load(id,mask);
  auto ptr=std::dynamic_pointer_cast<MassSpring<N>>(_saved[id]);
  if(mask&MASK_Z)
    _z=ptr->_z;
}
template <int N>
void MassSpring<N>::debugEnergy(int M,T tolG) {
  _y.setRandom(N,M);
  _yLast.setRandom(N,M);
  _L.setRandom(N,M);
  _Ax.setRandom(N,M);

  _k.resize(M);
  _l.resize(M);
  _lL.resize(M);
  _lH.resize(M);

  _x0=1;
  updateZ(tolG);
  for(int i=0; i<M; i++) {
    _k[i]=std::rand()/(T)RAND_MAX;
    _l[i]=std::rand()/(T)RAND_MAX;
    _lL[i]=_z.col(i).dot(_y.col(i))*std::rand()/(T)RAND_MAX;
    _lH[i]=_z.col(i).dot(_y.col(i))*(1+std::rand()/(T)RAND_MAX);
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
bool MassSpring<N>::energyY(const VecNT& y,T& E,VecNT* G,MatNT* H,int i) const {
  //energy
  T D=0,DD=0,D2=0,DD2=0,yLen=y.norm();
  E=0;
  if(!_evalgOnly) {
    E+=_beta*(y-_Ax.col(i)).squaredNorm()/2;
    E+=_L.col(i).dot(_Ax.col(i)-y);
    E+=_betaY*(y-_yLast.col(i)).squaredNorm()/2;
  }
  E+=_k[i]*(y-_z.col(i)*_l[i]).squaredNorm()/2;
  if(_lL[i]>0)
    E+=Penalty::eval<FLOAT>(_z.col(i).dot(y)-_lL[i],G?&D:NULL,H?&DD:NULL,0,1);
  if(_lH[i]>0)
    E+=Penalty::eval<FLOAT>(_lH[i]-yLen,G?&D2:NULL,H?&DD2:NULL,0,1);
  if(!isfinite(E))
    return false;
  //gradient
  if(G) {
    G->setZero();
    if(!_evalgOnly) {
      *G+=_beta*(y-_Ax.col(i));
      *G-=_L.col(i);
      *G+=_betaY*(y-_yLast.col(i));
    }
    *G+=_k[i]*(y-_z.col(i)*_l[i]);
    if(_lL[i]>0)
      *G+=D*_z.col(i);
    if(_lH[i]>0)
      *G-=D2*y/yLen;
  }
  //hessian
  if(H) {
    *H=MatNT::Identity();
    if(_evalgOnly)
      *H*=_k[i];
    else *H*=_beta+_betaY+_k[i];
    if(_lL[i]>0)
      *H+=DD*_z.col(i)*_z.col(i).transpose();
    if(_lH[i]>0) {
      *H-=(D2/yLen)*(MatNT::Identity()-y*y.transpose()/(yLen*yLen));
      *H+=(DD2/(yLen*yLen))*y*y.transpose();
    }
  }
  return true;
}
//instance
template class MassSpring<2>;
template class MassSpring<3>;
}
