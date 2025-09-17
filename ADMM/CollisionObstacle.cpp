#include "CollisionObstacle.h"
#include "CollisionDetector.h"
#include "SmallScaleNewton.h"
#include "SmallShape.h"
#include <Environment/GJK.h>
#include <Utils/DebugGradient.h>

namespace PHYSICSMOTION {
//Collision
template <int N,int M,int MO>
CollisionObstacle<N,M,MO>::CollisionObstacle(T r,T x0,T coef):_r(r),_coef(coef) {
  Penalty::_x0=(double)x0;
}
template <int N,int M,int MO>
const std::unordered_map<typename CollisionObstacle<N,M,MO>::ID,int>& CollisionObstacle<N,M,MO>::terms() const {
  return _terms;
}
template <int N,int M,int MO>
void CollisionObstacle<N,M,MO>::insertCollisions(const CollisionDetector<N,M>& detector) {
  int off=0;
  MatNMOT yObs;
  Eigen::Matrix<int,M,1> col;
  for(ID id:detector.obsCollisions())
    if(id!=CollisionDetector<N,M>::INVALID_ID) {
      //not already exists
      if(_terms.find(id)==_terms.end()) {
        //insert
        _terms[id]=off=(int)_terms.size();
        col=detector.obsCollisionId(id);
        off*=N*M;
        for(int i=0; i<M; i++,off+=N)
          addBlockId<T>(_E,off,col[i]*N,N,1);
        for(int i=0; i<MO; i++)
          yObs.col(i)=detector.obsCollisionPos(id,i);
        _yObs.push_back(yObs);
      }
    }
}
template <int N,int M,int MO>
int CollisionObstacle<N,M,MO>::removeCollisions(T margin) {
  std::vector<T> GJKResult(n());
  OMP_PARALLEL_FOR_
  for(int i=0; i<(int)GJKResult.size(); i++) {
    GJK::Vec3T pAL,pBL;
    bool intersect;
    Eigen::Matrix<T,N,M> a;
    SmallShape<N,M> A(a=_Ax.col(i).template segment<N*M>(0).reshaped(N,M));
    SmallShape<N,MO> B(_yObs[i]);
    GJKResult[i]=
      (T)GJK::runGJK(A,B,
                     GJK::Mat3X4T::Identity(),
                     GJK::Mat3X4T::Identity(),
                     pAL,pBL,&intersect);
  }
  //fill hash
  std::unordered_map<int,ID> invTerm;
  for(auto t:_terms)
    invTerm[t.second]=t.first;
  std::unordered_set<int> deleteHash;
  for(int i=0; i<(int)GJKResult.size(); i++)
    if(GJKResult[i]>margin)
      deleteHash.insert(i);
  //go ahead and remove
  if(deleteHash.empty())
    return (int)deleteHash.size();
  OptimizerTerm::removeColumns(deleteHash,N*M);
  //compact
  int newSize=0;
  _terms.clear();
  for(int i=0; i<(int)GJKResult.size(); i++)
    if(deleteHash.find(i)==deleteHash.end()) {
      _d0[newSize]=_d0[i];
      _z.col(newSize)=_z.col(i);
      _yObs[newSize]=_yObs[i];
      _terms[invTerm[i]]=newSize;
      newSize++;
    }
  //resize
  _d0.conservativeResize(newSize);
  _z.conservativeResize(_z.rows(),newSize);
  _yObs.resize(newSize);
  ASSERT(newSize==(int)_terms.size())
  return (int)deleteHash.size();
}
template <int N,int M,int MO>
typename CollisionObstacle<N,M,MO>::VecM CollisionObstacle<N,M,MO>::y0() {
  return VecM(_d0.data(),_d0.size());
}
template <int N,int M,int MO>
typename CollisionObstacle<N,M,MO>::VecCM CollisionObstacle<N,M,MO>::y0() const {
  return VecCM(_d0.data(),_d0.size());
}
template <int N,int M,int MO>
typename CollisionObstacle<N,M,MO>::VecCM CollisionObstacle<N,M,MO>::G0() const {
  return VecCM(_Gd0.data(),_Gd0.size());
}
template <int N,int M,int MO>
int CollisionObstacle<N,M,MO>::n() const {
  return (int)_terms.size();
}
template <int N,int M,int MO>
std::shared_ptr<OptimizerTerm> CollisionObstacle<N,M,MO>::copy() const {
  return std::shared_ptr<OptimizerTerm>(new CollisionObstacle<N,M,MO>(_r,Penalty::_x0));
}
template <int N,int M,int MO>
typename CollisionObstacle<N,M,MO>::T CollisionObstacle<N,M,MO>::evalG(bool calcG,bool initL,SMatT* H,int y0Off) {
  _evalgOnly=true;
  T g=0;
  if(calcG) {
    _G.resize(N*M,n());
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
    VecNMdT yd,G;
    CollisionMatrix<N,M> HBlk;
    yd.template segment<N*M>(0)=_Ax.col(i);
    yd[N*M]=_d0[i];
    energyYd(yd,E,&G,H?&HBlk:NULL,i);
    parallelAdd(g,E);
    if(calcG) {
      _G.col(i)=G.template segment<N*M>(0);
      _Gd0[i]=G[N*M];
    }
    if(H) {
      _HBlks[i]._blk=HBlk.toDense();
      _HBlks[i]._nY=N*M;
    }
  }
  if(calcG && initL)
    initializeL();
  if(H)
    assembleHessian(*H,y0Off);
  return g;
}
template <int N,int M,int MO>
typename CollisionObstacle<N,M,MO>::T CollisionObstacle<N,M,MO>::evalGDirect(bool calcG,SMatT* H,int y0Off,bool projPSD) {
  T g=0;
  if(calcG)
    _G.resize(N*M,n());
  if(H)
    _HBlks.resize(n());
  OMP_PARALLEL_FOR_
  for(int i=0; i<n(); i++) {
    T E;
    VecNMT y=_Ax.col(i),G;
    MatNMT HBlk;
    energyYDirect(y,E,&G,H?&HBlk:NULL,i,projPSD);
    parallelAdd(g,E);
    if(calcG) {
      _G.col(i)=G;
    }
    if(H) {
      _HBlks[i]._blk=HBlk;
      _HBlks[i]._nY=N*M;
    }
  }
  if(H)
    assembleHessian(*H,y0Off);
  return g;
}
template <int N,int M,int MO>
bool CollisionObstacle<N,M,MO>::updateY(T betaY,T beta,T tolG) {
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
    VecNMdT yd,G;
    CollisionMatrix<N,M> H;
    SmallScaleNewton<N*M+1,CollisionMatrix<N,M>> opt;
    yd.template segment<N*M>(0)=_y.col(i);
    yd[N*M]=_d0[i];
    auto energyFunc=[&](const VecNMdT& x,T& E,VecNMdT* G,CollisionMatrix<N,M>* H)->bool {
      return energyYd(x,E,G,H,i);
    };
    if(!opt.optimize(_alphaY[i],yd,E,G,H,energyFunc,tolG))
      succ=false;
    _y.col(i)=yd.template segment<N*M>(0);
    _d0[i]=yd[N*M];
  }
  return succ;
}
template <int N,int M,int MO>
bool CollisionObstacle<N,M,MO>::updateZ(T tolG) {
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
    //solve optimization to compute z
    T E;
    MatNT H;
    VecNT z,G;
    SmallScaleNewton<N,MatNT> opt;
#ifdef OPTIMIZE_ON_SPHERE
    auto energyFunc=[&](const VecNT& x,T& E,VecNT* G,MatNT* H)->bool {
      return energyZ(x,E,G,H,i);
    };
    if(!opt.optimizeOnSphere(_alphaZ[i],z=_z.col(i),E,G,H,energyFunc,tolG))
      succ=false;
#else
    auto energyFunc=[&](const VecNT& x,T& E,VecNT* G,MatNT* H)->bool {
      if(!energyZ(x,E,G,H,i))
        return false;
      if(!SmallScaleNewton<N,MatNT>::template energySoft<Penalty>(*this,x,E,G,H))
        return false;
      return true;
    };
    if(!opt.optimize(_alphaZ[i],z=_z.col(i),E,G,H,energyFunc,tolG))
      succ=false;
#endif
    _z.col(i)=z;
  }
  return succ;
}
template <int N,int M,int MO>
void CollisionObstacle<N,M,MO>::reset(int mask) {
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
template <int N,int M,int MO>
void CollisionObstacle<N,M,MO>::save(int id,int mask) {
  OptimizerTerm::save(id,mask);
  auto ptr=std::dynamic_pointer_cast<CollisionObstacle<N,M,MO>>(_saved[id]);
  if(mask&MASK_Y0)
    ptr->_d0=_d0;
  if(mask&MASK_Z)
    ptr->_z=_z;
}
template <int N,int M,int MO>
void CollisionObstacle<N,M,MO>::load(int id,int mask) {
  OptimizerTerm::load(id,mask);
  auto ptr=std::dynamic_pointer_cast<CollisionObstacle<N,M,MO>>(_saved[id]);
  if(mask&MASK_Y0)
    _d0=ptr->_d0;
  if(mask&MASK_Z) {
    _z=ptr->_z;
  }
}
template <int N,int M,int MO>
void CollisionObstacle<N,M,MO>::debugEnergy(int m,T tolG) {
  for(int i=0; i<m; i++) {
    Eigen::Matrix<T,N,M> y;
    Eigen::Matrix<T,N,MO> y2;
    for(int c=0; c<M; c++)
      y.col(c).setRandom();
    for(int c=0; c<MO; c++)
      y2.col(c)=VecNT::Random()+VecNT::Constant(2);
    _y=concatCol<MatT>(_y,y.reshaped());
    _yObs.push_back(y2);
    _terms[i]=i;
  }
  _yLast.setRandom(N*M,m);
  _L.setRandom(N*M,m);
  _Ax.setRandom(N*M,m);

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
      VecNMdT yd;
      yd.template segment<N*M>(0)=_y.col(i);
      yd[N*M]=_d0[i];
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
      CollisionMatrix<N,M> H;
      VecNMdT yd,yd2,dyd,G,G2;
      yd.template segment<N*M>(0)=_y.col(i);
      yd[N*M]=_d0[i];
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
}
template <int N,int M,int MO>
typename CollisionObstacle<N,M,MO>::T CollisionObstacle<N,M,MO>::eps() const {
  return _r+_x0*2;
}
//helper
template <int N,int M,int MO>
void CollisionObstacle<N,M,MO>::initializePlane(int i) {
  GJK::Vec3T pAL,pBL;
  bool intersect;
  Eigen::Matrix<T,N,M> a;
  SmallShape<N,M> A(a=_y.col(i).template segment<N*M>(0).reshaped(N,M));
  SmallShape<N,MO> B(_yObs[i]);
  GJK::runGJK(A,B,
              GJK::Mat3X4T::Identity(),
              GJK::Mat3X4T::Identity(),
              pAL,pBL,&intersect);
  _z.col(i)=(pAL-pBL).template segment<N>(0).template cast<T>().normalized();
  _d0[i]=-_z.col(i).dot((pAL+pBL).template segment<N>(0).template cast<T>())/2+_r/2;
}
template <int N,int M,int MO>
bool CollisionObstacle<N,M,MO>::energyYd(const VecNMdT& yd,T& E,VecNMdT* G,CollisionMatrix<N,M>* H,int i) const {
  E=0;
  if(!_evalgOnly) {
    E+=_beta*(yd.template segment<N*M>(0)-_Ax.col(i)).squaredNorm()/2;
    E+=_L.col(i).dot(_Ax.col(i)-yd.template segment<N*M>(0));
    E+=_betaY*(yd.template segment<N*M>(0)-_yLast.col(i)).squaredNorm()/2;
  }
  if(G) {
    G->setZero();
    if(!_evalgOnly) {
      G->template segment<N*M>(0)+=_beta*(yd.template segment<N*M>(0)-_Ax.col(i));
      G->template segment<N*M>(0)-=_L.col(i);
      G->template segment<N*M>(0)+=_betaY*(yd.template segment<N*M>(0)-_yLast.col(i));
    }
  }
  if(H) {
    H->setZero();
    H->setN(_z.col(i));
    if(!_evalgOnly)
      for(int i2=0; i2<M; i2++)
        H->addIdentity(i2,_beta+_betaY);
  }
  VecNT pos;
  int off=0,d0Off=N*M;
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
  for(int r=0; r<MO; r++) {
    pos=_yObs[i].col(r);
    E+=Penalty::eval<FLOAT>(-_z.col(i).dot(pos)-d0,G?&D:NULL,H?&DD:NULL,0,_coef);
    if(!isfinite(E))
      return false;
    if(G)
      G->coeffRef(d0Off)-=D;
    if(H)
      H->addIdentity(M,DD);
  }
  return true;
}
template <int N,int M,int MO>
bool CollisionObstacle<N,M,MO>::energyYDirect(const VecNMT& y,T& E,VecNMT* G,MatNMT* H,int i,bool projPSD) const {
  /*E = 0;
  if(G)
    G->setZero();
  if(H) {
    H->setZero();
  VecNT pos;
  T En;
  VecNdT Gn;
  MatNdT Hn;
  Gn.setZero();
  Hn.setZero();
  int off=0,d0Off=N*M;
  T D=0,DD=0,d0=_d0[i];
  //positive shape
  for(int r=0; r<M; r++,off+=N) {
    pos=yd.template segment<N>(off);
    E+=Penalty::eval<FLOAT>(_z.col(i).dot(pos)+d0-_r,G?&D:NULL,H?&DD:NULL,0,_coef);
    if(!isfinite(E))
      return false;
    if(G)
      G->template segment<N>(off)+=D*_z.col(i);
    if(H)
      H->template block<N,N>(off,off)+=DD*_z.col(i)*_z.col(i).transpose();
    Hn+=;
  }*/
  return true;
}
template <int N,int M,int MO>
bool CollisionObstacle<N,M,MO>::energyZ(const VecNT& z,T& E,VecNT* G,MatNT* H,int i) const {
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
  for(int r=0; r<MO; r++) {
    pos=_yObs[i].col(r);
    E+=Penalty::eval<FLOAT>(-z.dot(pos)-_d0[i],G?&D:NULL,H?&DD:NULL,0,_coef);
    if(!isfinite(E))
      return false;
    if(G)
      *G-=D*pos;
    if(H)
      *H+=DD*pos*pos.transpose();
  }
  return true;
}
//instance 2D
template class CollisionObstacle<2,1,2>;
template class CollisionObstacle<2,2,2>;
template class CollisionObstacle<2,2,4>;
template class CollisionObstacle<2,6,2>;
//instance 3D
template class CollisionObstacle<3,1,3>;
template class CollisionObstacle<3,3,3>;
template class CollisionObstacle<3,3,6>;
template class CollisionObstacle<3,6,3>;
}
