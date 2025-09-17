#include "Optimizer.h"
#include "ARAP.h"
#include "SmoothL1.h"
#include "MassSpring.h"
#include <Utils/SparseUtils.h>

namespace PHYSICSMOTION {
template <int N>
Optimizer<N>::Optimizer() {
  _gss.push_back(std::shared_ptr<OptimizerTerm>(new ARAP<N>()));
  _gss.push_back(std::shared_ptr<OptimizerTerm>(new MassSpring<N>()));
  _gss.push_back(std::shared_ptr<OptimizerTerm>(new SmoothL1<N>()));
}
template <int N>
void Optimizer<N>::setCB(std::function<bool()> cb) {
  _cb=cb;
}
template <int N>
void Optimizer<N>::addARAPElement(std::array<int,N+1> id,T k,const MatNT& F0) {
  findTerm<ARAP<N>>()->addElement(id,k,F0);
}
template <int N>
void Optimizer<N>::addSpring(std::array<int,2> id,T k,T l,T lL,T lH) {
  findTerm<MassSpring<N>>()->addSpring(id,k,l,lL,lH);
}
template <int N>
void Optimizer<N>::addSmoothL1(std::array<int,2> id,T k,T eps) {
  findTerm<SmoothL1<N>>()->addL1(id,k,eps);
}
template <int N>
void Optimizer<N>::addQuadratic(const SMatT* H,const Vec* g) {
  if(H)
    _H=*H;
  if(g)
    _g=*g;
}
template <int N>
void Optimizer<N>::assemble(int nX,const LinearConstraint& cons) {
  //assemble A
  for(const auto& g:Optimizer<N>::_gss)
    g->assembleA(&nX);
  //handle fix
  STrips id;
  id.clear();
  int nXF=0;
  for(const auto& f:cons._fixes)
    if(f.second.first>0)
      continue; //not a hard constraint
    else {
      addBlockId<T>(id,nXF,f.first*N,N,1.0f);
      nXF+=N;
    }
  //handle general
  if(cons._general.size()>0) {
    addBlock<T,0,int>(id,nXF,0,cons._general);
    nXF+=cons._general.rows();
  }
  //build KKT
  _Cons.resize(nXF,nX);
  _Cons.setFromTriplets(id.begin(),id.end());
  if(nXF>0)
    _CCTInv.compute(_Cons*_Cons.transpose());
}
template <int N>
void Optimizer<N>::insertCollisionDetector(std::shared_ptr<CollisionChecker<N>> coll) {
  _coll=coll;
}
template <int N>
void Optimizer<N>::save(int id,int mask) {
  if(mask&OptimizerTerm::MASK_X) {
    if(id>=(int)_xSaved.size())
      _xSaved.resize(id+1);
    _xSaved[id]=_x;
  }
  for(const auto& g:_gss)
    g->save(id,mask);
}
template <int N>
void Optimizer<N>::load(int id,int mask) {
  if(mask&OptimizerTerm::MASK_X) {
    ASSERT(id<(int)_xSaved.size())
    _x=_xSaved[id];
  }
  for(const auto& g:_gss)
    g->load(id,mask);
}
//helper
template <int N>
void Optimizer<N>::init(T tolG,bool directMode) {
  for(const auto& g:Optimizer<N>::_gss) {
    g->setDirectMode(false);
    g->reset();
    g->y()=g->Ax()=g->A()*Optimizer<N>::_x;
    g->updateZ(tolG);
    g->evalG(true,true);
    g->setDirectMode(directMode);
  }
}
template <int N>
void Optimizer<N>::project(Vec& G) {
  if(_Cons.rows()==0)
    return;
  _d=_CCTInv.solve(_Cons*G.segment(0,_Cons.cols()));
  G.segment(0,_Cons.cols())-=_Cons.transpose()*_d;
}
template <int N>
bool Optimizer<N>::updateZ(T tolG) {
  bool succ=true;
  for(const auto& g:_gss)
    if(!g->updateZ(tolG))
      succ=false;
  return succ;
}
//instance
template class Optimizer<2>;
template class Optimizer<3>;
}
