#include "DirectNewton.h"
#include <Utils/SparseUtils.h>
#include <Utils/Timing.h>

namespace PHYSICSMOTION {
static const int SLOT_BEFORE_Z_UPDATE = 0;
template <int N>
void DirectNewton<N>::optimize(const OptimizerParam& param) {
  Optimizer<N>::init(param._tolG,true);
  T alpha=param._initAlpha,E,E2;
  SMatT H;
  std::string collString;
  Vec x=GradientDescend<N>::assembleX(),x2,G,d;
  int nCollCheck=0;
  //timing
  double collTime=0,cbTime=0;
  while(!Optimizer<N>::_cb());
  TBEG();
  //collision information
  if(Optimizer<N>::_coll)
    collString=Optimizer<N>::_coll->info(*this);
  //main loop
  Optimizer<N>::save(SLOT_BEFORE_Z_UPDATE,OptimizerTerm::MASK_Z|OptimizerTerm::MASK_Y0);
  for(int iter=1; iter<param._maxIter; iter++) {
    //evaluate gradient
    E=DirectNewton<N>::evalGD(x,&G,&H);
    //find search direction
    Newton<N>::solve(d,x,G,H,param._psdEps);
    //clamp step size to avoid excessively large number of collisions
    /*if (iter == 1 && Optimizer<N>::_coll) {
      Vec sz=d.segment(0,Optimizer<N>::_x.size()).reshaped(N,Optimizer<N>::_x.size()/N).colwise().norm();
      alpha=std::min<T>(alpha,Optimizer<N>::_coll->getEps(*this)/std::max<T>(sz.maxCoeff(),Epsilon<T>::finiteDifferenceEps()));
	  std::cout << "Alpha clamped to " << alpha << " to avoid excessive collisions!" << std::endl;
    }*/
    //line search
    DirectNewton<N>::lineSearch(x2=x,E,d,alpha,param);
    //collision check
    if(Optimizer<N>::_coll) {
      TBEG();
      auto xLastCollFree=x.segment(0,Optimizer<N>::_x.size());
      auto xCurrCollFree=x2.segment(0,Optimizer<N>::_x.size());
      Eigen::Matrix<int,2,1> numTerms=Optimizer<N>::_coll->generate(xLastCollFree,xCurrCollFree,*this,true,true);
      if(numTerms[0]>0 || numTerms[1]>0) {
        collString=Optimizer<N>::_coll->info(*this);
        //revert to last x and reduce alpha
        alpha*=param._decAlpha*param._decAlpha;
        collTime+=TENDV();
        continue;
      }
      collTime+=TENDV();
    }
    //step-size too small
    if(alpha<param._tolAlpha)
      break;
    Optimizer<N>::_x=x2.segment(0,Optimizer<N>::_x.size());
    E2=DirectNewton<N>::evalGD(x2,NULL,NULL);
    //collision remove by energy value
    nCollCheck++;
    if(Optimizer<N>::_coll && param._collisionRemoveI>0 && nCollCheck%param._collisionRemoveI==0) {
      TBEG();
      auto xCurrCollFree=x.segment(0,Optimizer<N>::_x.size());
      Eigen::Matrix<int,2,1> numTerms=Optimizer<N>::_coll->removeByEnergy(xCurrCollFree,*this,Epsilon<T>::finiteDifferenceEps());
      if(numTerms[0]>0 || numTerms[1]>0)
        collString=Optimizer<N>::_coll->info(*this);
      collTime+=TENDV();
    }
    //Accept solution
    x=x2;
    Optimizer<N>::save(SLOT_BEFORE_Z_UPDATE,OptimizerTerm::MASK_Z|OptimizerTerm::MASK_Y0);
    //termination
    Optimizer<N>::project(G);
    if (isfinite(E) && G.cwiseAbs().maxCoeff() < param._tolG)
        break;
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
                << " " << collString
                << " TotalTime=" << TQUERYV() << " CollTime=" << collTime << " CBTime=" << cbTime << std::endl;
    }
    if(doDebugGradient)
      DirectNewton<N>::debugGradient(x);
  }
  //timing
  TENDV();
}
template <int N>
typename DirectNewton<N>::T DirectNewton<N>::evalGD(const Vec& x,Vec* G,SMatT* H,bool projPSD) {
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
  Optimizer<N>::load(SLOT_BEFORE_Z_UPDATE,OptimizerTerm::MASK_Z|OptimizerTerm::MASK_Y0);
  for(const auto& g:Optimizer<N>::_gss) {
    g->y()=g->Ax()=g->A()*x.segment(0,g->A().cols());
    int nY0=g->y0().size();
    if(nY0>0)
      g->y0()=x.segment(off,nY0);
    E+=g->evalGDirect(G!=NULL,H,off,projPSD);
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
void DirectNewton<N>::debugGradient(const Vec& x) {
  SMatT H;
  Vec dx,G,G2;
  T E=DirectNewton<N>::evalGD(x,&G,&H,false);
  dx=Vec::Random(G.size());

  DEFINE_NUMERIC_DELTA_T(T)
  T E2=DirectNewton<N>::evalGD(x+dx*DELTA,&G2,NULL,false);
  DEBUG_GRADIENT("G",G.dot(dx),G.dot(dx)-(E2-E)/DELTA)
  DEBUG_GRADIENT("H",(H*dx).norm(),(H*dx-(G2-G)/DELTA).norm())
}
template <int N>
void DirectNewton<N>::lineSearch(Vec& x,T& E,const Vec& D,T& alpha,const OptimizerParam& param) {
  Vec x2=x+D*alpha;
  T E2=DirectNewton<N>::evalGD(x2,NULL);
  while(E2>=E && alpha>param._tolAlpha) {
    alpha*=param._decAlpha;
    x2=x+D*alpha;
    E2=DirectNewton<N>::evalGD(x2,NULL);
  }
  if(alpha<=param._tolAlpha)
    return;
  x=x2;
  E=E2;
  alpha*=param._incAlpha;
}
//instance
template class DirectNewton<2>;
template class DirectNewton<3>;
}
