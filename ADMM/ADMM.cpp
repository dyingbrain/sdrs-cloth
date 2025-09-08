#include "ADMM.h"
#include <Utils/SparseUtils.h>
#include <Utils/Timing.h>

namespace PHYSICSMOTION {
template <int N>
ADMM<N>::ADMM():_beta(1),_betaX(0),_betaY(0) {}
template <int N>
void ADMM<N>::assemble(int nX,const LinearConstraint& cons) {
  Optimizer<N>::assemble(nX,cons);
  assembleLHS(nX);
  //handle fix
  int nXF=0;
  for(const auto& f:cons._fixes)
    if(f.second.first>0)
      continue; //not a hard constraint
    else nXF+=N;
  //handle general
  if(cons._general.size()>0)
    nXF+=cons._general.rows();
  //setup RHS
  _RHS.setZero(nX+nXF);
  nXF=0;
  for(const auto& f:cons._fixes)
    if(f.second.first>0)
      continue; //not a hard constraint
    else {
      _RHS.template segment<N>(nX+nXF)=f.second.second;
      nXF+=N;
    }
}
template <int N>
void ADMM<N>::optimize(const OptimizerParam& param) {
  static const int SLOT_BEFORE_Z_UPDATE=0;
  static const int SLOT_BEST_COLLISION_FREE=1;
  Vec G;
  std::string collString;
  _beta=param._initBeta;
  _betaX=param._initBetaX;
  _betaY=param._initBetaY;
  assembleLHS(Optimizer<N>::_x.size());
  Optimizer<N>::init(param._tolG);
  //timing
  double collTime=0,cbTime=0;
  while(!Optimizer<N>::_cb());
  TBEG();
  //track best, collision-free solution
  T EBest=eval(NULL),E;
  Optimizer<N>::save(SLOT_BEST_COLLISION_FREE,OptimizerTerm::MASK_X|OptimizerTerm::MASK_Y|OptimizerTerm::MASK_Y0|OptimizerTerm::MASK_L|OptimizerTerm::MASK_Z);
  //collision information
  if(Optimizer<N>::_coll)
    collString=Optimizer<N>::_coll->info(*this);
  //main loop
  int nrYUpdateFail=0,nrZUpdateFail=0,nrZUpdate=0,nCollCheck=0;
  for(int iter=1; iter<param._maxIter; iter++) {
    //updateX
    updateX();
    //updateY
    if(!updateY(param._tolGYZ))
      nrYUpdateFail++;
    //updateZ
    Optimizer<N>::save(SLOT_BEFORE_Z_UPDATE,OptimizerTerm::MASK_Z);
    if(!Optimizer<N>::updateZ(param._tolGYZ))
      nrZUpdateFail++;
    bool doCheck=(iter%param._checkI)==0;
    bool doPrint=(iter%param._printI)==0;
    E=eval((doCheck || doPrint)?&G:NULL);
    //track best, collision-free solution
    if(!param._ensureMonotonic || E<=EBest) {
      //collision check
      if(Optimizer<N>::_coll) {
        TBEG();
        const Vec& xLastCollFree=Optimizer<N>::_xSaved[SLOT_BEST_COLLISION_FREE];
        const Vec& xCurrCollFree=Optimizer<N>::_x;
        Eigen::Matrix<int,2,1> numTerms=Optimizer<N>::_coll->generate(xLastCollFree,xCurrCollFree,*this,true,true,SLOT_BEST_COLLISION_FREE);
        if(numTerms[0]>0 || numTerms[1]>0) {
          collString=Optimizer<N>::_coll->info(*this);
          //re-assemble LHS
          assembleLHS(Optimizer<N>::_x.size());
          //and re-evaluate EBest reflecting the added collision functions
          EBest=eval(NULL);
          collTime+=TENDV();
          continue;
        } else {
          nCollCheck++;
          if(param._collisionRemoveI>0 && nCollCheck%param._collisionRemoveI==0) {
            numTerms=Optimizer<N>::_coll->remove(xCurrCollFree,*this,param._collisionRemoveMargin);
            if(numTerms[0]>0 || numTerms[1]>0) {
              collString=Optimizer<N>::_coll->info(*this);
              //re-assemble LHS
              assembleLHS(Optimizer<N>::_x.size());
              //and re-evaluate EBest reflecting the added collision functions
              E=eval(NULL);
            }
          }
        }
        collTime+=TENDV();
      }
      Optimizer<N>::save(SLOT_BEST_COLLISION_FREE,OptimizerTerm::MASK_X|OptimizerTerm::MASK_Y|OptimizerTerm::MASK_Y0|OptimizerTerm::MASK_L|OptimizerTerm::MASK_Z);
      EBest=E;
      nrZUpdate++;
    } else Optimizer<N>::load(SLOT_BEFORE_Z_UPDATE,OptimizerTerm::MASK_Z);
    //termination
    if(doCheck) {
      //condition 1: xy too small
      if(param._tolXY>0) {
        T dx=(Optimizer<N>::_x-_xLast).cwiseAbs().maxCoeff(),dy=0;
        for(const auto& g:Optimizer<N>::_gss)
          if(g->y().size()>0)
            dy+=(g->y()-g->yLast()).cwiseAbs().maxCoeff();
        if(dx<param._tolXY && dy<param._tolXY)
          break;
      }
      //condition 2: g too small
      Optimizer<N>::project(G);
      if(G.cwiseAbs().maxCoeff()<param._tolG)
        break;
    }
    //print
    if(doPrint) {
      TBEG();
      if(EBest==E)
        while(!Optimizer<N>::_cb());
      cbTime+=TENDV();
      T gNorm=G.cwiseAbs().maxCoeff();
      std::cout << "Iter=" << iter
                << " E=" << E
                << " EBest=" << EBest
                << " gNorm=" << gNorm
                << " nrYFail=" << nrYUpdateFail
                << " nrZFail=" << nrZUpdateFail
                << " nrZUpdate=" << nrZUpdate
                << " " << collString
                << " TotalTime=" << TQUERYV() << " CollTime=" << collTime << " CBTime=" << cbTime << std::endl;
    }
    //update Lagrangian
    updateL();
  }
  //we need to return _xBest instead of _x because _x might be invalid
  Optimizer<N>::load(SLOT_BEST_COLLISION_FREE,OptimizerTerm::MASK_X);
  //timing
  TENDV();
}
//helper
template <int N>
void ADMM<N>::assembleLHS(int nX) {
  //assemble LHS
  STrips id;
  _LHS.resize(nX,nX);
  addBlockId(id,0,0,nX,_betaX);
  _LHS.setFromTriplets(id.begin(),id.end());
  //term by term
  _LHS+=Optimizer<N>::_H;
  for(const auto& g:Optimizer<N>::_gss)
    _LHS+=g->A().transpose()*g->A()*_beta;
  //build KKT
  _LHS=buildKKT<T,0,int>(_LHS,Optimizer<N>::_Cons,0);
  _LHSInv.compute(_LHS);
}
template <int N>
typename ADMM<N>::T ADMM<N>::eval(Vec* G) {
  T E=0;
  if(G)
    G->setZero(Optimizer<N>::_x.size());
  //term by term: f
  E+=Optimizer<N>::_x.dot(Optimizer<N>::_H*Optimizer<N>::_x)/2+Optimizer<N>::_g.dot(Optimizer<N>::_x);
  if(G)
    *G+=Optimizer<N>::_H*Optimizer<N>::_x+Optimizer<N>::_g;
  //term by term: g
  for(const auto& g:Optimizer<N>::_gss) {
    g->Ax()=g->A()*Optimizer<N>::_x;
    E+=g->evalG(G!=NULL,false);
    if(G)
      *G+=g->A().transpose()*g->G();
  }
  return E;
}
template <int N>
void ADMM<N>::updateX() {
  _xLast=Optimizer<N>::_x;
  _RHS.segment(0,Optimizer<N>::_x.size())=_betaX*_xLast;
  //term by term: f
  _RHS.segment(0,Optimizer<N>::_x.size())-=Optimizer<N>::_g;
  //term by term: g
  for(const auto& g:Optimizer<N>::_gss)
    _RHS.segment(0,Optimizer<N>::_x.size())+=g->A().transpose()*(_beta*g->y()-g->L());
  Optimizer<N>::_x=_LHSInv.solve(_RHS).segment(0,Optimizer<N>::_x.size());
}
template <int N>
bool ADMM<N>::updateY(T tolG) {
  bool succ=true;
  //term by term
  for(const auto& g:Optimizer<N>::_gss) {
    g->Ax()=g->A()*Optimizer<N>::_x;
    if(!g->updateY(_betaY,_beta,tolG))
      succ=false;
  }
  return succ;
}
template <int N>
void ADMM<N>::updateL() {
  //term by term
  for(const auto& g:Optimizer<N>::_gss)
    g->updateL(_beta);
}
//instance
template class ADMM<2>;
template class ADMM<3>;
}
