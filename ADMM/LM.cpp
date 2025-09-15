#include "LM.h"
#include <Utils/SparseUtils.h>
#include <Utils/Timing.h>

namespace PHYSICSMOTION {
template <int N>
void LM<N>::optimize(const OptimizerParam& param) {
  Optimizer<N>::init(param._tolG);
  T alpha=param._initAlpha,E,E2;
  SMatT H,HKKT;
  std::string collString;
  Vec x=LM<N>::assembleX(),x2,G,GKKT;
  int nrZUpdateFail=0,nCollCheck=0;
  //timing
  double collTime=0,cbTime=0;
  while(!Optimizer<N>::_cb());
  TBEG();
  //collision information
  if(Optimizer<N>::_coll)
    collString=Optimizer<N>::_coll->info(*this);
  //main loop
  for(int iter=1; iter<param._maxIter; iter++) {
    //evaluate gradient
    E=LM<N>::evalGD(x,&G,&H);
    //collision remove
    nCollCheck++;
    if(Optimizer<N>::_coll && param._collisionRemoveI>0 && nCollCheck%param._collisionRemoveI==0) {
      TBEG();
      auto xCurrCollFree=x.segment(0,Optimizer<N>::_x.size());
      Eigen::Matrix<int,2,1> numTerms=Optimizer<N>::_coll->remove(xCurrCollFree,*this,param._collisionRemoveMargin);
      if(numTerms[0]>0 || numTerms[1]>0) {
        collString=Optimizer<N>::_coll->info(*this);
        //revert to last x
        Optimizer<N>::_x=x.segment(0,Optimizer<N>::_x.size());
        x=LM<N>::assembleX();
        collTime+=TENDV();
        continue;
      }
      collTime+=TENDV();
    }
    //update solution using LM
    solveSchur(x2,x,G,H,alpha);
    //collision check
    if(Optimizer<N>::_coll) {
      TBEG();
      auto xLastCollFree=x.segment(0,Optimizer<N>::_x.size());
      auto xCurrCollFree=x2.segment(0,Optimizer<N>::_x.size());
      Eigen::Matrix<int,2,1> numTerms=Optimizer<N>::_coll->generate(xLastCollFree,xCurrCollFree,*this,true,true);
      if(numTerms[0]>0 || numTerms[1]>0) {
        collString=Optimizer<N>::_coll->info(*this);
        //revert to last x
        Optimizer<N>::_x=x.segment(0,Optimizer<N>::_x.size());
        x=LM<N>::assembleX();
        collTime+=TENDV();
        continue;
      }
      collTime+=TENDV();
    }
    //termination
    if(isfinite(E) && G.cwiseAbs().maxCoeff()<param._tolG)
      break;
    if(alpha>1/param._tolAlpha)
      break;
    //update
    E2=LM<N>::evalGD(x2,NULL,NULL);
    if(isfinite(E2) && E2<E) {
      x=x2;
      Optimizer<N>::_x=x.segment(0,Optimizer<N>::_x.size());
      alpha=std::max<T>(alpha*param._decAlpha,param._tolAlpha);
    } else {
      alpha=alpha*param._incAlpha;
      continue;
    }
    //print
    bool doPrint=(iter%param._printI)==0;
    bool doDebugGradient=param._debugGradientI>0 && (iter%param._debugGradientI)==0;
    if(doPrint) {
      TBEG();
      while(!Optimizer<N>::_cb());
      cbTime+=TENDV();
      project(G);
      T gNorm=G.cwiseAbs().maxCoeff();
      std::cout << "Iter=" << iter
                << " E=" << E
                << " gNorm=" << gNorm
                << " alpha=" << alpha
                << " nrZFail=" << nrZUpdateFail
                << " " << collString
                << " TotalTime=" << TQUERYV() << " CollTime=" << collTime << " CBTime=" << cbTime << std::endl;
    }
    if(doDebugGradient)
      LM<N>::debugGradient(x);
  }
  //timing
  TENDV();
}
template <int N>
typename LM<N>::Vec LM<N>::assembleX() {
  Vec x;
  return x;
}
template <int N>
typename LM<N>::T LM<N>::evalGD(const Vec& x,Vec* G,SMatT* H) {
  return 0;
}
template <int N>
void LM<N>::solveSchur(Vec& x2,const Vec& x,const Vec& G,const Vec& H,T alpha) {
}
template <int N>
void LM<N>::debugGradient(const Vec& x) {
}
//instance
template class LM<2>;
template class LM<3>;
}
