#include "ARAP.h"
#include "SmallScaleNewton.h"
#include <Utils/SparseUtils.h>

namespace PHYSICSMOTION {
template <int N>
void ARAP<N>::addElement(std::array<int,N+1> id,T k,const MatVT& F0) {
  for(int d=0; d<N+1; d++)
    addBlockId<T>(_E,(n()*(N+1)+d)*N,id[d]*N,N,1);
  _invF0.push_back(F0.inverse());
  _k.push_back(k);
}
template <int N>
typename ARAP<N>::VecM ARAP<N>::y0() {
  return VecM(_d0.data(),_d0.size());
}
template <int N>
typename ARAP<N>::VecCM ARAP<N>::y0() const {
  return VecCM(_d0.data(),_d0.size());
}
template <int N>
typename ARAP<N>::VecCM ARAP<N>::G0() const {
  return VecCM(_Gd0.data(),_Gd0.size());
}
template <int N>
int ARAP<N>::n() const {
  return (int)_k.size();
}
template <int N>
std::shared_ptr<OptimizerTerm> ARAP<N>::copy() const {
  return std::shared_ptr<OptimizerTerm>(new ARAP<N>);
}
template <int N>
typename ARAP<N>::T ARAP<N>::evalG(bool calcG,bool initL,SMatT* H,int y0Off) {
  if(initL) {
    //initialize activation range
    for(int i=0; i<n(); i++) {
      MatVT F0=_invF0[i].inverse();
      if(i==0)
        _x0=(double)minEdge(F0);
      else _x0=(double)std::min<T>(_x0,minEdge(F0));
    }
    //_x0*=0.1f;
  }
  _evalgOnly=true;
  T g=0;
  if(calcG) {
    _G.resize(N*(N+1),n());
    _Gd0.resize(N,n());
  }
  if(H)
    _HBlks.resize(n());
  OMP_PARALLEL_FOR_
  for(int i=0; i<n(); i++) {
    T E;
    MatYDT HBlk;
    VecYDT yd,G;
    yd.template segment<N*(N+1)>(0)=_Ax.col(i);
    yd.template segment<N>(N*(N+1))=_d0.col(i);
    energyYD(yd,E,calcG?&G:NULL,H?&HBlk:NULL,i);
    parallelAdd(g,E);
    if(calcG) {
      _G.col(i)=G.template segment<N*(N+1)>(0);
      _Gd0.col(i)=G.template segment<N>(N*(N+1));
    }
    if(H) {
      _HBlks[i]._blk=HBlk;
      _HBlks[i]._nY=N*(N+1);
    }
  }
  if(calcG && initL)
    initializeL();
  if(H)
    assembleHessian(*H,y0Off);
  return g;
}
template <int N>
bool ARAP<N>::updateY(T betaY,T beta,T tolG) {
  _evalgOnly=false;
  bool succ=true;
  _betaY=betaY;
  _beta=beta;
  _yLast=_y;
  if((int)_alphaY.size()<n())
    _alphaY.resize(n(),1);
  OMP_PARALLEL_FOR_
  for(int i=0; i<n(); i++) {
    T E;
    MatYDT H;
    VecYDT yd,G;
    SmallScaleNewton<N*(N+2),MatYDT> opt;
    yd.template segment<N*(N+1)>(0)=_y.col(i);
    yd.template segment<N>(N*(N+1))=_d0.col(i);
    auto energyFunc=[&](const VecYDT& x,T& E,VecYDT* G,MatYDT* H)->bool {
      return energyYD(x,E,G,H,i);
    };
    if(!opt.optimize(_alphaY[i],yd,E,G,H,energyFunc,tolG))
      succ=false;
    _y.col(i)=yd.template segment<N*(N+1)>(0);
    _d0.col(i)=yd.template segment<N>(N*(N+1));
  }
  return succ;
}
template <int N>
bool ARAP<N>::updateZ(T tolG) {
  bool succ=true;
  //initialize _d0
  if(_d0.cols()!=n()) {
    _d0.resize(N,n());
    OMP_PARALLEL_FOR_
    for(int i=0; i<n(); i++)
      _d0.col(i)=center(i);
  }
  //initialize _z,_alphaZ
  for(int d=0; d<N+1; d++) {
    if(_z[d].cols()!=n()) {
      _z[d].resize(N,n());
      OMP_PARALLEL_FOR_
      for(int i=0; i<n(); i++)
        _z[d].col(i)=normal(i,d);
    }
    if((int)_alphaZ[d].size()!=n())
      _alphaZ[d].resize(n(),1);
  }
  //initialize R
  if(_R.cols()!=n())
    _R.resize(N*N,n());
  OMP_PARALLEL_FOR_
  for(int i=0; i<n(); i++) {
    //solve 4 optimizations to compute z
    T E;
    MatVT H,F;
    VecVT z,G;
    SmallScaleNewton<N,MatVT> opt;
    for(int d=0; d<N+1; d++) {
#ifdef OPTIMIZE_ON_SPHERE
      auto energyFunc=[&](const VecVT& x,T& E,VecVT* G,MatVT* H)->bool {
        return energyZ(x,E,G,H,i,d);
      };
      if(!opt.optimizeOnSphere(_alphaZ[d][i],z=_z[d].col(i),E,G,H,energyFunc,tolG))
        succ=false;
#else
      auto energyFunc=[&](const VecVT& x,T& E,VecVT* G,MatVT* H)->bool {
        if(!energyZ(x,E,G,H,i,d))
          return false;
        if(!SmallScaleNewton<N,MatVT>::template energySoft<Penalty>(*this,x,E,G,H))
          return false;
        return true;
      };
      if(!opt.optimize(_alphaZ[d][i],z=_z[d].col(i),E,G,H,energyFunc,tolG))
        succ=false;
#endif
      _z[d].col(i)=z;
    }
    //polar decomposition to compute R
    for(int d=0; d<N; d++)
      F.col(d)=_y.col(i).template segment<N>((d+1)*N)-_y.col(i).template segment<N>(0);
    F*=_invF0[i];
    Eigen::JacobiSVD<Eigen::Matrix<double,N,N>> svd(F.template cast<double>(),Eigen::ComputeFullU|Eigen::ComputeFullV);
    Eigen::Map<MatVT>(_R.col(i).data())=(svd.matrixU()*svd.matrixV().transpose()).template cast<T>();
  }
  return succ;
}
template <int N>
void ARAP<N>::reset(int mask) {
  OptimizerTerm::reset(mask);
  if(mask&MASK_Y0)
    _d0.resize(N,0);
  if(mask&MASK_Z) {
    _R.resize(N*N,0);
    for(int i=0; i<N+1; i++)
      _z[i].resize(N,0);
  }
  if(mask&MASK_ALPHA) {
    _alphaY.clear();
    for(int i=0; i<N+1; i++)
      _alphaZ[i].clear();
  }
}
template <int N>
void ARAP<N>::save(int id,int mask) {
  OptimizerTerm::save(id,mask);
  auto ptr=std::dynamic_pointer_cast<ARAP<N>>(_saved[id]);
  if(mask&MASK_Y0)
    ptr->_d0=_d0;
  if(mask&MASK_Z) {
    ptr->_R=_R;
    for(int i=0; i<N+1; i++)
      ptr->_z[i]=_z[i];
  }
}
template <int N>
void ARAP<N>::load(int id,int mask) {
  OptimizerTerm::load(id,mask);
  auto ptr=std::dynamic_pointer_cast<ARAP<N>>(_saved[id]);
  if(mask&MASK_Y0)
    _d0=ptr->_d0;
  if(mask&MASK_Z) {
    _R=ptr->_R;
    for(int i=0; i<N+1; i++)
      _z[i]=ptr->_z[i];
  }
}
template <int N>
void ARAP<N>::debugEnergy(int M,T tolG) {
  _y.setRandom(N*(N+1),M);
  _yLast.setRandom(N*(N+1),M);
  _L.setRandom(N*(N+1),M);
  _Ax.setRandom(N*(N+1),M);

  _k.resize(M);
  _invF0.resize(M);
  for(int i=0; i<M; i++) {
    MatVT F;
    _k[i]=std::rand()/(T)RAND_MAX;
    for(int d=0; d<N; d++)
      F.col(d)=_y.col(i).template segment<N>((d+1)*N)-_y.col(i).template segment<N>(0);
    _invF0[i].setRandom();
  }

  _x0=1;
  updateZ(tolG);

  //we can now debug energyZ
  DEFINE_NUMERIC_DELTA_T(T)
  for(int i=0; i<M; i++) {
    std::cout << "DebugEnergyZ" << std::endl;
    for(int fid=0; fid<N+1; fid++) {
      T E,E2;
      MatVT H;
      VecVT z,z2,dz,G,G2;
      z=_z[fid].col(i);
      dz.setRandom();
      if(energyZ(z,E,&G,&H,i,fid)) {
        z2=z+dz*DELTA;
        energyZ(z2,E2,&G2,NULL,i,fid);
        DEBUG_GRADIENT("G",G.dot(dz),G.dot(dz)-(E2-E)/DELTA)
        DEBUG_GRADIENT("H",(H*dz).norm(),((G2-G)/DELTA-H*dz).norm())
      }
    }
  }

  _beta=std::rand()/(T)RAND_MAX;
  _betaY=std::rand()/(T)RAND_MAX;

  //we can now debug energyY
  _debug=true;
  for(bool evalgOnly: {
        false,true
      }) {
    std::cout << "DebugEnergyYD evalgOnly=" << evalgOnly << std::endl;
    for(int i=0; i<M; i++) {
      T E,E2;
      MatYDT H;
      VecYDT yd,yd2,dyd,G,G2;
      yd.template segment<N*(N+1)>(0)=_y.col(i);
      yd.template segment<N>(N*(N+1))=_d0.col(i);
      dyd.setRandom();
      _evalgOnly=evalgOnly;
      if(energyYD(yd,E,&G,&H,i)) {
        yd2=yd+dyd*DELTA;
        energyYD(yd2,E2,&G2,NULL,i);
        DEBUG_GRADIENT("G",G.dot(dyd),G.dot(dyd)-(E2-E)/DELTA)
        DEBUG_GRADIENT("H",(H*dyd).norm(),((G2-G)/DELTA-H*dyd).norm())
      }
    }
  }
  _debug=false;
}
//helper
template <typename T>
Eigen::Matrix<T,2,1> normalInner(const Eigen::Matrix<T,2,1> dir[2]) {
  Eigen::Matrix<T,2,1> n=(dir[1]-dir[0]).normalized();
  return Eigen::Matrix<T,2,1>(-n[1],n[0]);
}
template <typename T>
Eigen::Matrix<T,3,1> normalInner(const Eigen::Matrix<T,3,1> dir[3]) {
  return (dir[1]-dir[0]).cross(dir[2]-dir[0]).normalized();
}
template <int N>
bool ARAP<N>::energyYD(const VecYDT& yd,T& E,VecYDT* G,MatYDT* H,int i) const {
  //energy
  T D=0,DD=0;
  E=0;
  if(!_evalgOnly) {
    E+=_beta*(yd.template segment<N*(N+1)>(0)-_Ax.col(i)).squaredNorm()/2;
    E+=_L.col(i).dot(_Ax.col(i)-yd.template segment<N*(N+1)>(0));
    E+=_betaY*(yd.template segment<N*(N+1)>(0)-_yLast.col(i)).squaredNorm()/2;
  }
  //gradient
  if(G) {
    G->setZero();
    if(!_evalgOnly) {
      G->template segment<N*(N+1)>(0)+=_beta*(yd.template segment<N*(N+1)>(0)-_Ax.col(i));
      G->template segment<N*(N+1)>(0)-=_L.col(i);
      G->template segment<N*(N+1)>(0)+=_betaY*(yd.template segment<N*(N+1)>(0)-_yLast.col(i));
    }
  }
  //hessian
  if(H) {
    H->setZero();
    if(!_evalgOnly) {
      H->template block<N*(N+1),N*(N+1)>(0,0).setIdentity();
      *H*=_beta+_betaY;
    }
  }
  //ARAP
  Eigen::Map<const MatVT> R(_R.col(i).data());
  const auto& invF0=_invF0[i];
  for(int r=0; r<N; r++)
    for(int c=0; c<N; c++) {
      //term is k/2*|F.row(r)*invF0.col(c)-R(r,c)|^2
      int vid[4]= {0+r,N+r,N*2+r,N*3+r};
      T coef[4] = {-invF0.col(c).sum(),invF0(0,c),invF0(1,c),2>=N?0:invF0(2,c)},d=-R(r,c);
      for(int rr=0; rr<N+1; rr++)
        d+=yd[vid[rr]]*coef[rr];
      E+=d*d*_k[i]/2;
      if(G)
        for(int rr=0; rr<N+1; rr++)
          G->coeffRef(vid[rr])+=d*coef[rr]*_k[i];
      if(H)
        for(int rr=0; rr<N+1; rr++)
          for(int cc=0; cc<N+1; cc++)
            H->coeffRef(vid[rr],vid[cc])+=coef[rr]*coef[cc]*_k[i];
    }
  if(_evalgOnly && _debug) {
    MatVT F;
    for(int i2=0; i2<N; i2++)
      F.col(i2)=yd.template segment<N>((i2+1)*N)-yd.template segment<N>(0);
    MatVT DF=F*_invF0[i];
    T ERef=(DF-R).squaredNorm()*_k[i]/2;
    std::cout << "E=" << E << " ERef=" << ERef << std::endl;
  }
  //barrier
  VecVT dir;
  for(int fid=0; fid<N+1; fid++) {
    const auto& n=_z[fid].col(i);
    for(int d=1; d<=N; d++) {
      int vid=(fid+d)%(N+1),voff=vid*N;
      dir=yd.template segment<N>(voff)-yd.template segment<N>(N*(N+1));
      E+=Penalty::eval<FLOAT>(n.dot(dir),G?&D:NULL,H?&DD:NULL,0,1);
      if(!isfinite(E))
        return false;
      if(G) {
        G->template segment<N>(voff)+=D*n;
        G->template segment<N>(N*(N+1))-=D*n;
      }
      if(H) {
        H->template block<N,N>(voff,voff)+=DD*n*n.transpose();
        H->template block<N,N>(voff,N*(N+1))-=DD*n*n.transpose();
        H->template block<N,N>(N*(N+1),voff)-=DD*n*n.transpose();
        H->template block<N,N>(N*(N+1),N*(N+1))+=DD*n*n.transpose();
      }
    }
  }
  return true;
}
template <int N>
bool ARAP<N>::energyZ(const VecVT& z,T& E,VecVT* G,MatVT* H,int tetId,int fid) const {
  //energy
  if(E)
    E=0;
  if(G)
    G->setZero();
  if(H)
    H->setZero();
  for(int d=1; d<=N; d++) {
    T D=0,DD=0;
    int vid=(fid+d)%(N+1),voff=vid*N;
    VecVT dir=_y.col(tetId).template segment<N>(voff)-_d0.col(tetId);
    E+=Penalty::eval<FLOAT>(z.dot(dir),G?&D:NULL,H?&DD:NULL,0,1);
    if(G)
      *G+=D*dir;
    if(H)
      *H+=DD*dir*dir.transpose();
  }
  return isfinite(E);
}
template <int N>
typename ARAP<N>::VecVT ARAP<N>::center(int tetId) const {
  VecVT ctr=VecVT::Zero();
  for(int d=0; d<=N; d++) {
    int voff=d*N;
    ctr+=_y.col(tetId).template segment<N>(voff);
  }
  return ctr/(N+1);
}
template <int N>
typename ARAP<N>::VecVT ARAP<N>::normal(int tetId,int fid) const {
  //collect all directions
  VecVT dir[N];
  for(int d=1; d<=N; d++) {
    int vid=(fid+d)%(N+1),voff=vid*N;
    dir[d-1]=_y.col(tetId).template segment<N>(voff)-_d0.col(tetId);
  }
  //compute normal
  VecVT n=normalInner(dir);
  if(n.dot(dir[0])<0)
    n*=-1;
  return n;
}
template <int N>
typename ARAP<N>::T ARAP<N>::minEdge(const MatVT& F0) const {
  std::vector<T> elss;
  for(int d=0; d<N; d++)
    elss.push_back(F0.col(d).norm());
  for(int d=0; d<N; d++)
    elss.push_back((F0.col(d)-F0.col((d+1)%N)).norm());
  return *std::min_element(elss.begin(),elss.end());
}
//instance
template class ARAP<2>;
template class ARAP<3>;
}
