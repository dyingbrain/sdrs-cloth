#include "GradientDescend.h"
#include <Utils/SparseUtils.h>
#include <Utils/Timing.h>

namespace PHYSICSMOTION {
template <int N>
void GradientDescend<N>::optimize(const OptimizerParam& param) {
  Optimizer<N>::init(param._tolG);
  T alpha=param._initAlpha,E;
  Vec x=assembleX(),x2,G;
  std::string collString;
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
    E=evalGD(x,&G);
    Optimizer<N>::project(G);
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
        x=assembleX();
        collTime+=TENDV();
        continue;
      }
      collTime+=TENDV();
    }
    //line search
    lineSearch(x2=x,E,-G,alpha,param);
    evalGD(x,NULL); //we need to reset y0 for collision detection
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
        x=assembleX();
        collTime+=TENDV();
        continue;
      }
      collTime+=TENDV();
    }
    x=x2;
    //termination
    if(G.cwiseAbs().maxCoeff()<param._tolG)
      break;
    if(alpha<=param._tolAlpha)
      break;
    //print
    bool doPrint=(iter%param._printI)==0;
    bool doDebugGradient=param._debugGradientI>0 && (iter%param._debugGradientI)==0;
    if(doPrint) {
      TBEG();
      while(!Optimizer<N>::_cb());
      cbTime+=TENDV();
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
      debugGradient(x);
    //updateZ
    if(!updateZ(x,param._tolGYZ))
      nrZUpdateFail++;
  }
  //timing
  TENDV();
}
template <int N>
typename GradientDescend<N>::Vec GradientDescend<N>::assembleX() {
  Vec x=Optimizer<N>::_x;
  for(const auto& g:Optimizer<N>::_gss) {
    int nY0=g->y0().size();
    if(nY0>0)
      x=concat<Vec>(x,g->y0());
  }
  return x;
}
template <int N>
typename GradientDescend<N>::T GradientDescend<N>::evalGD(const Vec& x,Vec* G,SMatT* H) {
  T E=0;
  if(G)
    G->setZero(x.size());
  //term by term: f
  int off=GradientDescend<N>::_g.size();
  E+=x.segment(0,off).dot(Optimizer<N>::_H*x.segment(0,off))/2+Optimizer<N>::_g.dot(x.segment(0,off));
  if(G)
    G->segment(0,off)+=Optimizer<N>::_H*x.segment(0,off)+Optimizer<N>::_g;
  if(H) {
    STrips HTrips;
    H->resize(x.size(),x.size());
    addBlock(HTrips,0,0,Optimizer<N>::_H);
    H->setFromTriplets(HTrips.begin(),HTrips.end());
  }
  //term by term: g
  for(const auto& g:Optimizer<N>::_gss) {
    g->y()=g->Ax()=g->A()*x.segment(0,g->A().cols());
    int nY0=g->y0().size();
    if(nY0>0)
      g->y0()=x.segment(off,nY0);
    E+=g->evalG(G!=NULL,false,H,off);
    if(G) {
      G->segment(0,g->A().cols())+=g->A().transpose()*g->G();
      if(nY0>0)
        G->segment(off,nY0)+=g->G0();
    }
    off+=g->y0().size();
  }
  return E;
}
template <int N>
bool GradientDescend<N>::updateZ(const Vec& x,T tolG) {
  int off=GradientDescend<N>::_g.size();
  for(const auto& g:Optimizer<N>::_gss) {
    g->y()=g->Ax()=g->A()*x.segment(0,g->A().cols());
    int nY0=g->y0().size();
    if(nY0>0)
      g->y0()=x.segment(off,nY0);
    off+=g->y0().size();
  }
  return Optimizer<N>::updateZ(tolG);
}
template <int N>
void GradientDescend<N>::debugGradient(const Vec& x) {
  SMatT H;
  Vec dx,G,G2;
  T E=evalGD(x,&G,&H);
  dx=Vec::Random(G.size());

  DEFINE_NUMERIC_DELTA_T(T)
  T E2=evalGD(x+dx*DELTA,&G2);
  DEBUG_GRADIENT("G",G.dot(dx),G.dot(dx)-(E2-E)/DELTA)
  DEBUG_GRADIENT("H",(H*dx).norm(),(H*dx-(G2-G)/DELTA).norm())
}
template <int N>
void GradientDescend<N>::lineSearch(Vec& x,T& E,const Vec& D,T& alpha,const OptimizerParam& param) {
  Vec x2=x+D*alpha;
  T E2=evalGD(x2,NULL);
  while(E2>=E && alpha>param._tolAlpha) {
    alpha*=param._decAlpha;
    x2=x+D*alpha;
    E2=evalGD(x2,NULL);
  }
  if(alpha<=param._tolAlpha)
    return;
  x=x2;
  E=E2;
  alpha*=param._incAlpha;
}
//instance
template class GradientDescend<2>;
template class GradientDescend<3>;
}
