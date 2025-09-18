#include "CollisionSelf.h"
#include "CollisionDetector.h"
#include "SmallScaleNewton.h"
#include "SmallShape.h"
#include <Environment/GJK.h>
#include <Utils/DebugGradient.h>

namespace PHYSICSMOTION {
//Collision
template <int N,int M>
CollisionSelf<N,M>::CollisionSelf(T r,T x0,T coef):_r(r),_coef(coef) {
  Penalty::_x0=(double)x0;
}
template <int N,int M>
const std::unordered_map<typename CollisionSelf<N,M>::ID,int>& CollisionSelf<N,M>::terms() const {
  return _terms;
}
template <int N,int M>
void CollisionSelf<N,M>::insertCollisions(const CollisionDetector<N,M>& detector) {
  int off=0;
  Eigen::Matrix<int,M*2,1> col;
  for(ID id:detector.selfCollisions())
    if(id!=CollisionDetector<N,M>::INVALID_ID) {
      //not already exists
      if(_terms.find(id)==_terms.end()) {
        //insert
        _terms[id]=off=(int)_terms.size();
        col=detector.selfCollisionId(id);
        off*=N*M*2;
        for(int i=0; i<M*2; i++,off+=N)
          addBlockId<T>(_E,off,col[i]*N,N,1);
      }
    }
}
template <int N,int M>
int CollisionSelf<N,M>::removeCollisions(const std::unordered_set<int>& deleteHash) {
  int nOld=n();
  //fill hash
  std::unordered_map<int,ID> invTerm;
  for(auto t:_terms)
    invTerm[t.second]=t.first;
  //go ahead and remove
  if(deleteHash.empty())
    return (int)deleteHash.size();
  OptimizerTerm::removeColumns(deleteHash,N*M*2);
  //compact
  int newSize=0;
  _terms.clear();
  for(int i=0; i<nOld; i++)
    if(deleteHash.find(i)==deleteHash.end()) {
      _d0[newSize]=_d0[i];
      _z.col(newSize)=_z.col(i);
      _terms[invTerm[i]]=newSize;
      newSize++;
    }
  //resize
  _d0.conservativeResize(newSize);
  _z.conservativeResize(_z.rows(),newSize);
  ASSERT(newSize==(int)_terms.size())
  return (int)deleteHash.size();
}
template <int N,int M>
int CollisionSelf<N,M>::removeCollisionsByDistance(T margin) {
  std::vector<T> GJKResult(n());
  OMP_PARALLEL_FOR_
  for(int i=0; i<(int)GJKResult.size(); i++) {
    GJK::Vec3T pAL,pBL;
    bool intersect;
    Eigen::Matrix<T,N,M> a,b;
    SmallShape<N,M> A(a=_Ax.col(i).template segment<N*M>(0).reshaped(N,M));
    SmallShape<N,M> B(b=_Ax.col(i).template segment<N*M>(N*M).reshaped(N,M));
    GJKResult[i]=
      (T)GJK::runGJK(A,B,
                     GJK::Mat3X4T::Identity(),
                     GJK::Mat3X4T::Identity(),
                     pAL,pBL,&intersect);
  }
  //fill hash
  std::unordered_set<int> deleteHash;
  for(int i=0; i<(int)GJKResult.size(); i++)
    if(GJKResult[i]>margin)
      deleteHash.insert(i);
  return removeCollisions(deleteHash);
}
template <int N,int M>
int CollisionSelf<N,M>::removeCollisionsByEnergy(T thres) {
  updateZ(Epsilon<T>::finiteDifferenceEps());
  //Start estimation
  std::vector<T> ess(n(),0);
  OMP_PARALLEL_FOR_
  for(int i=0; i<n(); i++) {
    VecNMMT y=_Ax.col(i);
    energyYDirect(y,ess[i],NULL,NULL,i,false);
  }
  //fill hash
  std::unordered_set<int> deleteHash;
  for(int i=0; i<(int)ess.size(); i++)
    if(ess[i]<thres)
      deleteHash.insert(i);
  return removeCollisions(deleteHash);
}
template <int N,int M>
typename CollisionSelf<N,M>::VecM CollisionSelf<N,M>::y0() {
  if(_directMode)
    return VecM(NULL,0);
  else return VecM(_d0.data(),_d0.size());
}
template <int N,int M>
typename CollisionSelf<N,M>::VecCM CollisionSelf<N,M>::y0() const {
  if(_directMode)
    return VecCM(NULL,0);
  else return VecCM(_d0.data(),_d0.size());
}
template <int N,int M>
typename CollisionSelf<N,M>::VecCM CollisionSelf<N,M>::G0() const {
  return VecCM(_Gd0.data(),_Gd0.size());
}
template <int N,int M>
int CollisionSelf<N,M>::n() const {
  return (int)_terms.size();
}
template <int N,int M>
std::shared_ptr<OptimizerTerm> CollisionSelf<N,M>::copy() const {
  return std::shared_ptr<OptimizerTerm>(new CollisionSelf<N,M>(_r,Penalty::_x0));
}
template <int N,int M>
typename CollisionSelf<N,M>::T CollisionSelf<N,M>::evalG(bool calcG,bool initL,SMatT* H,int y0Off) {
  _evalgOnly=true;
  T g=0;
  if(calcG) {
    _G.resize(N*M*2,n());
    _Gd0.resize(n());
  }
  if(H)
    _HBlks.resize(n());
  //this allow updating only the new L
  int off=0;
  if(calcG && initL)
    off=_L.cols();
  OMP_PARALLEL_FOR_
  for(int i=off; i<n(); i++) {
    T E;
    VecNMMdT yd,G;
    CollisionMatrix<N,M*2> HBlk;
    yd.template segment<N*M*2>(0)=_Ax.col(i);
    yd[N*M*2]=_d0[i];
    energyYd(yd,E,&G,H?&HBlk:NULL,i);
    parallelAdd(g,E);
    if(calcG) {
      _G.col(i)=G.template segment<N*M*2>(0);
      _Gd0[i]=G[N*M*2];
    }
    if(H) {
      _HBlks[i]._blk=HBlk.toDense();
      _HBlks[i]._nY=N*M*2;
    }
  }
  if(calcG && initL)
    initializeL();
  if(H)
    assembleHessian(*H,y0Off);
  return g;
}
template <int N,int M>
typename CollisionSelf<N,M>::T CollisionSelf<N,M>::evalGDirect(bool calcG,SMatT* H,int y0Off,bool projPSD) {
  T g=0;
  if(calcG)
    _G.resize(N*M*2,n());
  if(H)
    _HBlks.resize(n());
  assert(_directMode);
  updateZ(Epsilon<T>::finiteDifferenceEps());
  //Start estimation
  OMP_PARALLEL_FOR_
  for(int i=0; i<n(); i++) {
    T E;
    VecNMMT y=_Ax.col(i),G;
    MatNMMT HBlk;
    energyYDirect(y,E,&G,H?&HBlk:NULL,i,projPSD);
    parallelAdd(g,E);
    if(calcG) {
      _G.col(i)=G;
    }
    if(H) {
      _HBlks[i]._blk=HBlk;
      _HBlks[i]._nY=N*M*2;
    }
  }
  if(H)
    assembleHessian(*H,y0Off);
  return g;
}
template <int N,int M>
bool CollisionSelf<N,M>::updateY(T betaY,T beta,T tolG) {
  _evalgOnly=false;
  bool succ=true;
  _betaY=betaY;
  _beta=beta;
  _yLast=_y;
  if((int)_alphaY.size()!=n())
    _alphaY.resize(n(),1);
  OMP_PARALLEL_FOR_
  for(int i=0; i<n(); i++) {
    T E;
    VecNMMdT yd,G;
    CollisionMatrix<N,M*2> H;
    SmallScaleNewton<N*M*2+1,CollisionMatrix<N,M*2>> opt;
    yd.template segment<N*M*2>(0)=_y.col(i);
    yd[N*M*2]=_d0[i];
    auto energyFunc=[&](const VecNMMdT& x,T& E,VecNMMdT* G,CollisionMatrix<N,M*2>* H)->bool {
      return energyYd(x,E,G,H,i);
    };
    if(!opt.optimize(_alphaY[i],yd,E,G,H,energyFunc,tolG))
      succ=false;
    _y.col(i)=yd.template segment<N*M*2>(0);
    _d0[i]=yd[N*M*2];
  }
  return succ;
}
template <int N,int M>
bool CollisionSelf<N,M>::updateZ(T tolG) {
  bool succ=true;
  //initialize _z
  if(_z.cols()!=n()) {
    int currSz=_z.cols();
    _z.conservativeResize(N,n());
    _d0.conservativeResize(n());
    OMP_PARALLEL_FOR_
    for(int j=currSz; j<n(); j++)
      initializePlane(j);
    //this means we are only updating term, then we just need initializeZ, no need to update it
    if(tolG==-1)
      return succ;
  }
  if((int)_alphaZ.size()!=n())
    _alphaZ.resize(n(),1);
  OMP_PARALLEL_FOR_
  for(int i=0; i<n(); i++) {
    if(_directMode) {
      //solve optimization to compute z and d0
      SmallScaleNewton<N+1,MatNdT> opt;
      auto energyFunc=[&](const VecNdT& x,T& E,VecNdT* G,MatNdT* H)->bool {
        if(!energyZd(x,E,G,H,i))
          return false;
        if(!SmallScaleNewton<N,MatNdT>::template energySoft<Penalty,VecNdT>(*this,x,E,G,H))
          return false;
        return true;
      };
      T E;
      MatNdT H;
      VecNdT z,G;
      z.template segment<N>(0)=_z.col(i);
      z[N]=_d0[i];
      int newtonIter=_directMode?100:1; //A heuristic value
      if(!opt.optimize(_alphaZ[i],z,E,G,H,energyFunc,tolG,newtonIter))
        succ=false;
      _z.col(i)=z.template segment<N>(0);
      _d0[i]=z[N];
    } else {
      //solve optimization to compute z
      SmallScaleNewton<N,MatNT> opt;
      auto energyFunc=[&](const VecNT& x,T& E,VecNT* G,MatNT* H)->bool {
        if(!energyZ(x,E,G,H,i))
          return false;
        if(!SmallScaleNewton<N,MatNT>::template energySoft<Penalty>(*this,x,E,G,H))
          return false;
        return true;
      };
      T E;
      MatNT H;
      VecNT z,G;
      int newtonIter=_directMode?100:1; //A heuristic value
      if(!opt.optimize(_alphaZ[i],z=_z.col(i),E,G,H,energyFunc,tolG,newtonIter))
        succ=false;
      _z.col(i)=z;
    }
  }
  return succ;
}
template <int N,int M>
void CollisionSelf<N,M>::reset(int mask) {
  OptimizerTerm::reset(mask);
  if(mask&MASK_Y0)
    _d0.resize(0);
  if(mask&MASK_Z)
    _z.resize(N,0);
  if(mask&MASK_ALPHA) {
    _alphaY.clear();
    _alphaZ.clear();
  }
}
template <int N,int M>
void CollisionSelf<N,M>::save(int id,int mask) {
  OptimizerTerm::save(id,mask);
  auto ptr=std::dynamic_pointer_cast<CollisionSelf<N,M>>(_saved[id]);
  if(mask&MASK_Y0)
    ptr->_d0=_d0;
  if(mask&MASK_Z)
    ptr->_z=_z;
}
template <int N,int M>
void CollisionSelf<N,M>::load(int id,int mask) {
  OptimizerTerm::load(id,mask);
  auto ptr=std::dynamic_pointer_cast<CollisionSelf<N,M>>(_saved[id]);
  if(mask&MASK_Y0)
    _d0=ptr->_d0;
  if(mask&MASK_Z)
    _z=ptr->_z;
}
template <int N,int M>
void CollisionSelf<N,M>::debugEnergy(int m,T tolG) {
  for(int i=0; i<m; i++) {
    Eigen::Matrix<T,N,M*2> y;
    for(int c=0; c<M; c++)
      y.col(c).setRandom();
    for(int c=M; c<M*2; c++)
      y.col(c)=VecNT::Random()+VecNT::Constant(2);
    _y=concatCol<MatT>(_y,y.reshaped());
    _terms[i]=i;
  }
  _yLast.setRandom(N*M*2,m);
  _L.setRandom(N*M*2,m);
  _Ax.setRandom(N*M*2,m);

  _evalgOnly=true;
  updateZ(tolG);
  //we are now ready to debug energyZ
  DEFINE_NUMERIC_DELTA_T(T)
  std::cout << "DebugEnergyZ" << std::endl;
  for(int i=0; i<m; i++) {
    T E,E2;
    MatNT H;
    VecNT z,z2,dz,G,G2;
    z=_z.col(i);
    dz.setRandom();
    if(energyZ(z,E,&G,&H,i)) {
      //debug energy
      VecNMMdT yd;
      yd.template segment<N*M*2>(0)=_y.col(i);
      yd[N*M*2]=_d0[i];
      energyYd(yd,E2,NULL,NULL,i);
      DEBUG_GRADIENT("E",E,E-E2)
      //derivative energy
      z2=z+dz*DELTA;
      energyZ(z2,E2,&G2,NULL,i);
      DEBUG_GRADIENT("G",G.dot(dz),G.dot(dz)-(E2-E)/DELTA)
      DEBUG_GRADIENT("H",(H*dz).norm(),((G2-G)/DELTA-H*dz).norm())
    }
  }

  _beta=std::rand()/(T)RAND_MAX;
  _betaY=std::rand()/(T)RAND_MAX;

  //we are now ready to debug energyY
  for(bool evalgOnly: {
        false,true
      }) {
    std::cout << "DebugEnergyYD evalgOnly=" << evalgOnly << std::endl;
    for(int i=0; i<m; i++) {
      T E,E2;
      CollisionMatrix<N,M*2> H;
      VecNMMdT yd,yd2,dyd,G,G2;
      yd.template segment<N*M*2>(0)=_y.col(i);
      yd[N*M*2]=_d0[i];
      dyd.setRandom();
      _evalgOnly=evalgOnly;
      energyYd(yd,E,&G,&H,i);
      //derivative energy
      yd2=yd+dyd*DELTA;
      energyYd(yd2,E2,&G2,NULL,i);
      DEBUG_GRADIENT("G",G.dot(dyd),G.dot(dyd)-(E2-E)/DELTA)
      DEBUG_GRADIENT("H",(H.toDense()*dyd).norm(),((G2-G)/DELTA-H.toDense()*dyd).norm())
    }
  }
  updateY(_betaY,_beta,tolG);
  
  bool projPSD=false;   //Otherwise, we are using Gauss-Newton hessian, which is not exact
  bool savedDirectMode=_directMode;
  _directMode=true; //Optimize both n and d0
  //we finally debug energyYDDirect
  std::cout << "DebugEnergyYDirect" << std::endl;
  //since we are using the inverse function theorem, the finite-difference epsilon cannot be too high
  DELTA = Epsilon<double>::finiteDifferenceEps();
  for(int i=0; i<m; i++) {
    T E,E2;
    MatNMMT H;
    VecNMMT yd,yd2,dyd,G,G2;
    yd=_y.col(i);
    dyd.setRandom();
    updateZ(Epsilon<T>::finiteDifferenceEps());
    energyYDirect(yd,E,&G,&H,i,projPSD);
    //derivative energy
    yd2=yd+dyd*DELTA;
    _y.col(i)=yd2;
    updateZ(Epsilon<T>::finiteDifferenceEps());
    energyYDirect(yd2,E2,&G2,NULL,i,projPSD);
    DEBUG_GRADIENT("G",G.dot(dyd),G.dot(dyd)-(E2-E)/DELTA)
    DEBUG_GRADIENT("H",(H*dyd).norm(),((G2-G)/DELTA-H*dyd).norm())
  }
  _directMode=savedDirectMode;
}
template <int N,int M>
typename CollisionSelf<N,M>::T CollisionSelf<N,M>::eps() const {
  return (_r+_x0)*2;
}
//helper
template <int N,int M>
void CollisionSelf<N,M>::initializePlane(int i) {
  GJK::Vec3T pAL,pBL;
  bool intersect;
  Eigen::Matrix<T,N,M> a,b;
  SmallShape<N,M> A(a=_y.col(i).template segment<N*M>(0).reshaped(N,M));
  SmallShape<N,M> B(b=_y.col(i).template segment<N*M>(N*M).reshaped(N,M));
  GJK::runGJK(A,B,
              GJK::Mat3X4T::Identity(),
              GJK::Mat3X4T::Identity(),
              pAL,pBL,&intersect);
  _z.col(i)=(pAL-pBL).template segment<N>(0).template cast<T>().normalized();
  _d0[i]=-_z.col(i).dot((pAL+pBL).template segment<N>(0).template cast<T>())/2;
}
template <int N,int M>
bool CollisionSelf<N,M>::energyYd(const VecNMMdT& yd,T& E,VecNMMdT* G,CollisionMatrix<N,M*2>* H,int i) const {
  E=0;
  if(!_evalgOnly) {
    E+=_beta*(yd.template segment<N*M*2>(0)-_Ax.col(i)).squaredNorm()/2;
    E+=_L.col(i).dot(_Ax.col(i)-yd.template segment<N*M*2>(0));
    E+=_betaY*(yd.template segment<N*M*2>(0)-_yLast.col(i)).squaredNorm()/2;
  }
  if(G) {
    G->setZero();
    if(!_evalgOnly) {
      G->template segment<N*M*2>(0)+=_beta*(yd.template segment<N*M*2>(0)-_Ax.col(i));
      G->template segment<N*M*2>(0)-=_L.col(i);
      G->template segment<N*M*2>(0)+=_betaY*(yd.template segment<N*M*2>(0)-_yLast.col(i));
    }
  }
  if(H) {
    H->setZero();
    H->setN(_z.col(i));
    if(!_evalgOnly)
      for(int i2=0; i2<M*2; i2++)
        H->addIdentity(i2,_beta+_betaY);
  }
  VecNT pos;
  int off=0,d0Off=N*M*2;
  T D=0,DD=0,d0=yd[d0Off];
  //positive shape
  for(int r=0; r<M; r++,off+=N) {
    pos=yd.template segment<N>(off);
    E+=Penalty::eval<FLOAT>(_z.col(i).dot(pos)+d0-_r,G?&D:NULL,H?&DD:NULL,0,_coef);
    if(!isfinite(E))
      return false;
    if(G) {
      G->template segment<N>(off)+=D*_z.col(i);
      G->coeffRef(d0Off)+=D;
    }
    if(H)
      H->addNNT(r,DD);
  }
  //negative shape
  for(int r=0; r<M; r++,off+=N) {
    pos=yd.template segment<N>(off);
    E+=Penalty::eval<FLOAT>(-_z.col(i).dot(pos)-d0-_r,G?&D:NULL,H?&DD:NULL,0,_coef);
    if(!isfinite(E))
      return false;
    if(G) {
      G->template segment<N>(off)-=D*_z.col(i);
      G->coeffRef(d0Off)-=D;
    }
    if(H)
      H->addNNT(r+M,DD);
  }
  return true;
}
template <int N,int M>
bool CollisionSelf<N,M>::energyYDirect(const VecNMMT& y,T& E,VecNMMT* G,MatNMMT* H,int i,bool projPSD) const {
  E=0;
  if(G)
    G->setZero();
  if(H)
    H->setZero();
  T D=0,DD=0;
  MatNdT Hnd;
  VecNdT nd,pos1=VecNdT::Ones();
  nd.template segment<N>(0)=_z.col(i);
  nd[N]=_d0[i];
  Hnd.setZero();
  T En=0;
  if(!SmallScaleNewton<N,MatNdT>::template energySoft<Penalty,VecNdT>(*this,nd,En,NULL,&Hnd))
    return false;
  E+=En;
  //positive shape
  MatNNdT DDInv[M*2];
  for(int r=0,off=0;r<M;r++,off+=N) {
    pos1.template segment<N>(0)=y.template segment<N>(off);
    E+=Penalty::eval<FLOAT>(nd.dot(pos1)-_r,G?&D:NULL,H?&DD:NULL,0,_coef);
    DDInv[r]=DD*_z.col(i)*pos1.transpose();
    DDInv[r].template block<N,N>(0,0)+=D*MatNT::Identity();
    if(!isfinite(E))
      return false;
    if(G)
      G->template segment<N>(off)+=D*_z.col(i);
    if(H) {
      H->template block<N,N>(off,off)+=DD*_z.col(i)*_z.col(i).transpose();
      if(!projPSD)
        Hnd+=pos1*pos1.transpose()*DD;
    }
  }
  //negative shape
  for(int r=0,off=N*M;r<M;r++,off+=N) {
    pos1.template segment<N>(0)=y.template segment<N>(off);
    E+=Penalty::eval<FLOAT>(-nd.dot(pos1)-_r,G?&D:NULL,H?&DD:NULL,0,_coef);
    DDInv[r+M]=DD*_z.col(i)*pos1.transpose();
    DDInv[r+M].template block<N,N>(0,0)-=D*MatNT::Identity();
    if(!isfinite(E))
      return false;
    if(G)
      G->template segment<N>(off)-=D*_z.col(i);
    if(H) {
      H->template block<N,N>(off,off)+=DD*_z.col(i)*_z.col(i).transpose();
      if(!projPSD)
        Hnd+=pos1*pos1.transpose()*DD;
    }
  }
  if(H && !projPSD) {
    //Compute hessian with respect to normal and d0
    for(int i=0;i<=N;i++)
      Hnd(i,i)=std::max(Hnd(i,i),Epsilon<T>::finiteDifferenceEps());
    Hnd=Hnd.inverse().eval();
    for(int r=0,off=0;r<M*2;r++,off+=N) 
      for(int r2=0,off2=0;r2<M*2;r2++,off2+=N)
        H->template block<N,N>(off,off2)+=-DDInv[r]*Hnd*DDInv[r2].transpose();
  }
  if(H && projPSD) {
    //Compute hessian with respect to d0
    Hnd(N,N)=std::max(Hnd(N,N),Epsilon<T>::finiteDifferenceEps());
    Hnd(N,N)=1/Hnd(N,N);
    for(int r=0,off=0;r<M*2;r++,off+=N) 
      for(int r2=0,off2=0;r2<M*2;r2++,off2+=N)
        H->template block<N,N>(off,off2)+=-DDInv[r].col(N)*Hnd(N,N)*DDInv[r2].col(N).transpose();
  }
  return true;
}
template <int N,int M>
bool CollisionSelf<N,M>::energyZ(const VecNT& z,T& E,VecNT* G,MatNT* H,int i) const {
  E=0;
  if(G)
    G->setZero();
  if(H)
    H->setZero();
  VecNT pos;
  int off=0;
  T D=0,DD=0;
  //positive shape
  for(int r=0; r<M; r++,off+=N) {
    pos=_y.template block<N,1>(off,i);
    E+=Penalty::eval<FLOAT>(z.dot(pos)+_d0[i]-_r,G?&D:NULL,H?&DD:NULL,0,_coef);
    if(!isfinite(E))
      return false;
    if(G)
      *G+=D*pos;
    if(H)
      *H+=DD*pos*pos.transpose();
  }
  //negative shape
  for(int r=0; r<M; r++,off+=N) {
    pos=_y.template block<N,1>(off,i);
    E+=Penalty::eval<FLOAT>(-z.dot(pos)-_d0[i]-_r,G?&D:NULL,H?&DD:NULL,0,_coef);
    if(!isfinite(E))
      return false;
    if(G)
      *G-=D*pos;
    if(H)
      *H+=DD*pos*pos.transpose();
  }
  return true;
}
template <int N,int M>
bool CollisionSelf<N,M>::energyZd(const VecNdT& zd,T& E,VecNdT* G,MatNdT* H,int i) const {
  E=0;
  if(G)
    G->setZero();
  if(H)
    H->setZero();
  VecNdT pos1=VecNdT::Ones();
  int off=0;
  T D=0,DD=0;
  //positive shape
  for(int r=0;r<M;r++,off+=N) {
    pos1.template segment<N>(0)=_y.template block<N,1>(off,i);
    E+=Penalty::eval<FLOAT>(zd.dot(pos1)-_r,G?&D:NULL,H?&DD:NULL,0,_coef);
    if(!isfinite(E))
      return false;
    if(G)
      *G+=D*pos1;
    if(H)
      *H+=DD*pos1*pos1.transpose();
  }
  //negative shape
  for(int r=0;r<M;r++,off+=N) {
    pos1.template segment<N>(0)=_y.template block<N,1>(off,i);
    E+=Penalty::eval<FLOAT>(-zd.dot(pos1)-_r,G?&D:NULL,H?&DD:NULL,0,_coef);
    if(!isfinite(E))
      return false;
    if(G)
      *G-=D*pos1;
    if(H)
      *H+=DD*pos1*pos1.transpose();
  }
  return true;
}
//instance 2D
template class CollisionSelf<2,1>;
template class CollisionSelf<2,2>;
template class CollisionSelf<2,6>;
//instance 3D
template class CollisionSelf<3,1>;
template class CollisionSelf<3,3>;
template class CollisionSelf<3,6>;
}
