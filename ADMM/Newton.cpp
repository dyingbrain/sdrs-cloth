#include "Newton.h"
#include <Utils/SparseUtils.h>
#include <Utils/Timing.h>

namespace PHYSICSMOTION {
template <int N>
void Newton<N>::optimize(const OptimizerParam& param) {
  Optimizer<N>::init(param._tolG);
  T alpha=1/param._initAlpha,E,E2;
  SMatT H;
  std::string collString;
  Vec x=GradientDescend<N>::assembleX(),x2,d,G;
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
    E=GradientDescend<N>::evalGD(x,&G,&H);
    //collision remove
    nCollCheck++;
    if(Optimizer<N>::_coll && param._collisionRemoveI>0 && nCollCheck%param._collisionRemoveI==0) {
      TBEG();
      auto xCurrCollFree=x.segment(0,Optimizer<N>::_x.size());
      Eigen::Matrix<int,2,1> numTerms=Optimizer<N>::_coll->removeByDistance(xCurrCollFree,*this,param._collisionRemoveMargin);
      if(numTerms[0]>0 || numTerms[1]>0) {
        collString=Optimizer<N>::_coll->info(*this);
        //revert to last x
        Optimizer<N>::_x=x.segment(0,Optimizer<N>::_x.size());
        x=GradientDescend<N>::assembleX();
        collTime+=TENDV();
        continue;
      }
      collTime+=TENDV();
    }
    //solve Newton's direction
    solve(d,x,G,H,alpha);
    //We are use LM, so just update d
    x2=x+d;
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
        x=GradientDescend<N>::assembleX();
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
    E2=GradientDescend<N>::evalGD(x2,NULL,NULL);
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
      GradientDescend<N>::project(G);
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
      GradientDescend<N>::debugGradient(x);
    //updateZ
    if(!GradientDescend<N>::updateZ(x,param._tolGYZ))
      nrZUpdateFail++;
  }
  //timing
  TENDV();
}
template <int N>
void Newton<N>::solve(Vec& d,const Vec& x,const Vec& G,const SMatT& H,T alpha){
  _GKKT=concat<Vec>(G,Vec::Zero(Optimizer<N>::_Cons.rows()));
  _HKKT=buildKKT<T,0,int>(H,Optimizer<N>::_Cons,alpha);
  //Select solver based on PSD property of matrix
  if(Optimizer<N>::_Cons.size()==0) {
    //std::cout << "BEG-LDLT-Solve size=" << _HKKT.rows() << std::endl;
    _invHSym.compute(_HKKT);
    d=-_invHSym.solve(_GKKT);
    //std::cout << "END-LDLT-Solve" << std::endl;
  } else {
    //std::cout << "BEG-LU-Solve size=" << _HKKT.rows() << std::endl;
    _invH.compute(_HKKT);
    d=-_invH.solve(_GKKT).segment(0,x.size());
    //std::cout << "END-LU-Solve" << std::endl;
  }
}
//instance
template class Newton<2>;
template class Newton<3>;
}
