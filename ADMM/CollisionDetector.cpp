#include "CollisionDetector.h"
#include "CollisionObstacle.h"
#include "CollisionSelf.h"
#include "Narrow.h"
#include <Utils/Utils.h>

namespace PHYSICSMOTION {
//CollisionDetector
template <int N,int M>
CollisionDetector<N,M>::CollisionDetector(int nCP,const MeshExact& obs):_obs(obs) {
  int order=M-1;
  int numSegment=(nCP-1)/order;
  //build segment indices
  std::deque<int> roots;
  for(int i=0; i<numSegment; i++) {
    Eigen::Matrix<int,M,1> seg;
    for(int d=0; d<=order; d++)
      seg[d]=i*order+d;
    _meshIss.push_back(seg);
    //build leaves
    Node<int,BBoxExact> n;
    n._l=n._r=n._parent=-1;
    n._nrCell=1;
    n._cell=i;
    roots.push_back((int)_meshBVH.size());
    _meshBVH.push_back(n);
  }
  //build internal nodes
  while(roots.size()>1) {
    //merge
    Node<int,BBoxExact> n;
    //connect l
    n._l=roots.front();
    _meshBVH[n._l]._parent=(int)_meshBVH.size();
    roots.pop_front();
    //connect r
    n._r=roots.front();
    _meshBVH[n._r]._parent=(int)_meshBVH.size();
    roots.pop_front();
    //other parts
    n._parent=-1;
    n._nrCell=_meshBVH[n._l]._nrCell+_meshBVH[n._r]._nrCell;
    n._cell=-1;
    roots.push_back((int)_meshBVH.size());
    _meshBVH.push_back(n);
  }

  //layer reordering for parallel computation
  Node<int,BBoxExact>::layerReorder(_meshBVH,_layerOffsets);
  ASSERT_MSG(_layerOffsets[1]==(int)_meshIss.size(),"LayerOffset error!");

  //parallel computation
  _nProc=omp_get_num_procs();
  _sss.resize(_nProc);
}
template <int N,int M>
CollisionDetector<N,M>::CollisionDetector(const MeshExact& mesh,const MeshExact& obs):_obs(obs) {
  _meshBVH=mesh.getBVH();
  _meshIss.resize(mesh.iss().size());
  int sz=std::min<int>(M,3);
  for(int i=0; i<(int)mesh.iss().size(); i++) {
    _meshIss[i].setZero();
    _meshIss[i].segment(0,sz)=mesh.iss()[i].segment(0,sz);
  }

  //layer reordering for parallel computation
  Node<int,BBoxExact>::layerReorder(_meshBVH,_layerOffsets);
  ASSERT_MSG(_layerOffsets[1]==(int)_meshIss.size(),"LayerOffset error!");
  //Node<int,BBoxExact>::parityCheck(_meshBVH);
  //for(int i=0; i<(int)_layerOffsets.size()-1; i++)
  //  std::cout << "Layer " << i << " has " << _layerOffsets[i+1]-_layerOffsets[i] << " nodes!" << std::endl;

  //parallel computation
  _nProc=omp_get_num_procs();
  _sss.resize(_nProc);
}
template <int N,int M>
void CollisionDetector<N,M>::generateBroadBF(const Vec& x,const Vec& x2,T eps,T epsObs,bool noNeighbor,bool self,bool obs) {
  _pssSelf.clear();
  _pssObs.clear();
  //update BVH leaves
  OMP_PARALLEL_FOR_
  for(int i=0; i<(int)_meshBVH.size(); i++)
    if(_meshBVH[i]._cell>=0)
      updateBVH(x,x2,eps,i);
  //self
  if(self && !_meshIss.empty())
    for(int l=0; l<(int)_meshIss.size(); l++)
      for(int r=l+1; r<(int)_meshIss.size(); r++) {
        //intersection test
        const BBoxExact& bbl=_meshBVH[l]._bb;
        const BBoxExact& bbr=_meshBVH[r]._bb;
        if(!bbl.intersect(bbr))
          continue;
        //add entry
        bool isNeighbor=false;
        const Eigen::Matrix<int,M,1>& tl=_meshIss[l];
        const Eigen::Matrix<int,M,1>& tr=_meshIss[r];
        if(noNeighbor)
          for(int i=0; i<M; i++)
            for(int j=0; j<M; j++)
              if(tl[i]==tr[j])
                isNeighbor=true;
        if(!isNeighbor)
          _pssSelf.push_back(combine(l,r));
      }
  //obs
  if(obs && !_meshIss.empty() && !_obs.iss().empty()) {
    BBoxExact::Vec3T pos=BBoxExact::Vec3T::Zero();
    pos.template segment<N>(0).setConstant((BBoxExact::T)epsObs);
    for(int l=0; l<(int)_meshIss.size(); l++) {
      BBoxExact bbl=_meshBVH[l]._bb.enlarged(pos);
      for(int r=0; r<(int)_obs.iss().size(); r++) {
        //intersection test
        const BBoxExact& bbr=_obs.getBVH()[r]._bb;
        if(!bbl.intersect(bbr))
          continue;
        //add entry
        _pssObs.push_back(combine(l,r));
      }
    }
  }
}
template <int N,int M>
void CollisionDetector<N,M>::generateBroad(const Vec& x,const Vec& x2,T eps,T epsObs,bool noNeighbor,bool parallel,bool self,bool obs) {
  _pssSelf.clear();
  _pssObs.clear();
  //update BVH nodes
  updateBVH(x,x2,eps,parallel);
  //self
  if(self && !_meshIss.empty())
    generateSelf(combine(_meshBVH.size()-1,_meshBVH.size()-1),x,noNeighbor,parallel);
  //obs
  if(obs && !_meshIss.empty() && !_obs.iss().empty())
    generateObs(combine(_meshBVH.size()-1,_obs.getBVH().size()-1),x,epsObs,parallel);
}
template <int N,int M>
Eigen::Matrix<int,2,1> CollisionDetector<N,M>::generate(const Vec& x,const Vec& x2,Optimizer<N>& opt,bool noNeighbhor,bool parallel,int restoreSlot) {
  std::shared_ptr<CollisionSelf<N,M>> self=opt.template findTerm<CollisionSelf<N,M>>();
  std::shared_ptr<CollisionObstacle<N,M>> obs=opt.template findTerm<CollisionObstacle<N,M>>();
  //broad phase
  T eps=self?self->eps()/2:0,epsObs=std::max<T>(obs->eps()-eps,0);
  generateBroad(x,x2,eps,epsObs,noNeighbhor,parallel,self!=NULL,obs!=NULL);
  //narrow phase
  Eigen::Matrix<int,2,1> ret=Eigen::Matrix<int,2,1>::Zero();
  Narrow<N,M> narrow(_meshIss,_obs);
  if(self)
    ret[0]=narrow.self(x,x2,_pssSelf,self->terms(),self->eps());
  if(obs)
    ret[1]=narrow.obs(x,x2,_pssObs,obs->terms(),obs->eps());
  //restore ADMM state
  if(ret.sum()>0 && restoreSlot>=0)
    opt.load(restoreSlot,OptimizerTerm::MASK_X|OptimizerTerm::MASK_Y|OptimizerTerm::MASK_Y0|OptimizerTerm::MASK_L|OptimizerTerm::MASK_Z);
  //insert energy
  if(ret[0]>0) {
    int ARows=self->A().rows();
    int ACols=self->A().cols();
    self->insertCollisions(*this);
    self->assembleA(NULL);
    int ARowsNew=self->A().rows()-ARows;
    //updateZ at position x
    self->y().segment(ARows,ARowsNew)=
      self->Ax().segment(ARows,ARowsNew)=
        self->A().block(ARows,0,ARowsNew,ACols)*x.segment(0,ACols);
    self->updateZ(-1);  //notify that we are only initializing Z
    self->evalG(true,true,NULL,0);
    //update ADMM state
    if(restoreSlot>=0)
      self->save(restoreSlot,OptimizerTerm::MASK_Y|OptimizerTerm::MASK_Y0|OptimizerTerm::MASK_L|OptimizerTerm::MASK_Z);
  }
  if(ret[1]>0) {
    int ARows=obs->A().rows();
    int ACols=obs->A().cols();
    obs->insertCollisions(*this);
    obs->assembleA(NULL);
    int ARowsNew=obs->A().rows()-ARows;
    //updateZ at position x
    obs->y().segment(ARows,ARowsNew)=
      obs->Ax().segment(ARows,ARowsNew)=
        obs->A().block(ARows,0,ARowsNew,ACols)*x.segment(0,ACols);
    obs->updateZ(-1);  //notify that we are only initializing Z
    obs->evalG(true,true,NULL,0);
    //update ADMM state
    if(restoreSlot>=0)
      obs->save(restoreSlot,OptimizerTerm::MASK_Y|OptimizerTerm::MASK_Y0|OptimizerTerm::MASK_L|OptimizerTerm::MASK_Z);
  }
  return ret;
}
template <int N,int M>
Eigen::Matrix<int,2,1> CollisionDetector<N,M>::removeByDistance(const Vec& x,Optimizer<N>& opt,T margin) {
  std::shared_ptr<CollisionSelf<N,M>> self=opt.template findTerm<CollisionSelf<N,M>>();
  std::shared_ptr<CollisionObstacle<N,M>> obs=opt.template findTerm<CollisionObstacle<N,M>>();
  Eigen::Matrix<int,2,1> ret=Eigen::Matrix<int,2,1>::Zero();
  if(self) {
    self->Ax()=self->A()*x.segment(0,self->A().cols());
    ret[0]=self->removeCollisionsByDistance(self->eps()*margin);
  }
  if(obs) {
    obs->Ax()=obs->A()*x.segment(0,obs->A().cols());
    ret[1]=obs->removeCollisionsByDistance(obs->eps()*margin);
  }
  return ret;
}
template <int N,int M>
Eigen::Matrix<int,2,1> CollisionDetector<N,M>::removeByEnergy(const Vec& x,Optimizer<N>& opt,T thres) {
  std::shared_ptr<CollisionSelf<N,M>> self=opt.template findTerm<CollisionSelf<N,M>>();
  std::shared_ptr<CollisionObstacle<N,M>> obs=opt.template findTerm<CollisionObstacle<N,M>>();
  Eigen::Matrix<int,2,1> ret=Eigen::Matrix<int,2,1>::Zero();
  if(self) {
    self->Ax()=self->A()*x.segment(0,self->A().cols());
    ret[0]=self->removeCollisionsByEnergy(thres);
  }
  if(obs) {
    obs->Ax()=obs->A()*x.segment(0,obs->A().cols());
    ret[1]=obs->removeCollisionsByEnergy(thres);
  }
  return ret;
}
template <int N,int M>
std::string CollisionDetector<N,M>::info(const Optimizer<N>& opt) const {
  std::shared_ptr<CollisionSelf<N,M>> self=opt.template findTerm<CollisionSelf<N,M>>();
  std::shared_ptr<CollisionObstacle<N,M>> obs=opt.template findTerm<CollisionObstacle<N,M>>();
  std::ostringstream oss;
  if(self)
    oss << "#SelfColl=" << self->n();
  if(self && obs)
    oss << " ";
  if(obs)
    oss << "#ObsColl=" << obs->n();
  return oss.str();
}
template <int N,int M>
typename CollisionDetector<N,M>::T CollisionDetector<N,M>::getEps(Optimizer<N>& opt) const {
  std::shared_ptr<CollisionSelf<N,M>> self=opt.template findTerm<CollisionSelf<N,M>>();
  std::shared_ptr<CollisionObstacle<N,M>> obs=opt.template findTerm<CollisionObstacle<N,M>>();
  T eps=self?self->eps()/2:0,epsObs=std::max<T>(obs->eps()-eps,0);
  return std::max<T>(eps,epsObs);
}
//getter
template <int N,int M>
void CollisionDetector<N,M>::extractPos(const MeshExact& m,Vec& x) const {
  x.resize((int)m.vss().size()*N);
  for(int i=0; i<(int)m.vss().size(); i++)
    x.template segment<N>(i*N)=m.vss()[i].template segment<N>(0).template cast<T>();
}
template <int N,int M>
const std::vector<typename CollisionDetector<N,M>::ID>& CollisionDetector<N,M>::selfCollisions() const {
  return _pssSelf.getVector();
}
template <int N,int M>
const std::vector<typename CollisionDetector<N,M>::ID>& CollisionDetector<N,M>::obsCollisions() const {
  return _pssObs.getVector();
}
template <int N,int M>
Eigen::Matrix<int,M*2,1> CollisionDetector<N,M>::selfCollisionId(ID c) const {
  int l,r;
  separate(c,l,r);

  Eigen::Matrix<int,M*2,1> ret;
  ret.template segment<M>(0)=_meshIss[l].template segment<M>(0);
  ret.template segment<M>(M)=_meshIss[r].template segment<M>(0);
  return ret;
}
template <int N,int M>
Eigen::Matrix<int,M,1> CollisionDetector<N,M>::obsCollisionId(ID c) const {
  int l,r;
  separate(c,l,r);

  Eigen::Matrix<int,M,1> ret;
  ret.template segment<M>(0)=_meshIss[l].template segment<M>(0);
  return ret;
}
template <int N,int M>
typename CollisionDetector<N,M>::VecNT CollisionDetector<N,M>::obsCollisionPos(ID c,int j) const {
  int l,r;
  separate(c,l,r);
  return _obs.vss()[_obs.iss()[r][j]].template segment<N>(0).template cast<T>();
}
//ID
template <int N,int M>
void CollisionDetector<N,M>::separate(ID pair,int& l,int& r) {
  static const ID mask=((ID)1<<(ID)32)-1;
  l=(int)(pair>>32);
  r=(int)(pair&mask);
}
template <int N,int M>
typename CollisionDetector<N,M>::ID CollisionDetector<N,M>::combine(ID l,ID r) {
  return (l<<32)+r;
}
//helper
template <int N,int M>
void CollisionDetector<N,M>::generateSelf(ID nodePair,const Vec& x,bool noNeighbor,bool parallel) {
  typedef std::chrono::high_resolution_clock Clock;
  //serial computation
  int l,r;
  STACK& ss=_sss[omp_get_thread_num()];
  ss.assign(1,nodePair);
  while(!ss.empty()) {
    separate(ss.front(),l,r);
    ss.pop_front();
    const auto& nl=_meshBVH[l];
    const auto& nr=_meshBVH[r];
    if(!nl._bb.intersect(nr._bb))
      continue;
    else if(nl._cell>=0 && nr._cell>=0) {
      if(nl._cell<nr._cell) {
        bool isNeighbor=false;
        const Eigen::Matrix<int,M,1>& tl=_meshIss[nl._cell];
        const Eigen::Matrix<int,M,1>& tr=_meshIss[nr._cell];
        if(noNeighbor)
          for(int i=0; i<M; i++)
            for(int j=0; j<M; j++)
              if(tl[i]==tr[j])
                isNeighbor=true;
        if(!isNeighbor)
          _pssSelf.push_back(combine(nl._cell,nr._cell));
      }
    } else if(nl._cell>=0) {
      ss.push_back(combine(l,nr._l));
      ss.push_back(combine(l,nr._r));
    } else if(nr._cell>=0) {
      ss.push_back(combine(nl._l,r));
      ss.push_back(combine(nl._r,r));
    } else {
      ss.push_back(combine(nl._l,nr._l));
      ss.push_back(combine(nl._l,nr._r));
      ss.push_back(combine(nl._r,nr._l));
      ss.push_back(combine(nl._r,nr._r));
    }
    if(parallel && (int)ss.size()>=_nProc)
      break;
  }
  //parallel computation
  if(!ss.empty()) {
    _nodePairs.assign(ss.begin(),ss.end());
    OMP_PARALLEL_FOR_
    for(int i=0; i<(int)_nodePairs.size(); i++)
      generateSelf(_nodePairs[i],x,noNeighbor,false);
  }
}
template <int N,int M>
void CollisionDetector<N,M>::generateObs(ID nodePair,const Vec& x,T epsObs,bool parallel) {
  BBoxExact::Vec3T pos=BBoxExact::Vec3T::Zero();
  pos.template segment<N>(0).setConstant((BBoxExact::T)epsObs);
  //serial computation
  int l,r;
  STACK& ss=_sss[omp_get_thread_num()];
  ss.assign(1,nodePair);
  while(!ss.empty()) {
    separate(ss.front(),l,r);
    ss.pop_front();
    const auto& nl=_meshBVH[l];
    const auto& nr=_obs.getBVH()[r];
    if(!nl._bb.intersect(nr._bb.enlarged(pos)))
      continue;
    else if(nl._cell>=0 && nr._cell>=0) {
      _pssObs.push_back(combine(nl._cell,nr._cell));
    } else if(nl._cell>=0) {
      ss.push_back(combine(l,nr._l));
      ss.push_back(combine(l,nr._r));
    } else if(nr._cell>=0) {
      ss.push_back(combine(nl._l,r));
      ss.push_back(combine(nl._r,r));
    } else {
      ss.push_back(combine(nl._l,nr._l));
      ss.push_back(combine(nl._l,nr._r));
      ss.push_back(combine(nl._r,nr._l));
      ss.push_back(combine(nl._r,nr._r));
    }
    if(parallel && (int)ss.size()>=_nProc)
      break;
  }
  //parallel computation
  if(!ss.empty()) {
    _nodePairs.assign(ss.begin(),ss.end());
    OMP_PARALLEL_FOR_
    for(int i=0; i<(int)_nodePairs.size(); i++)
      generateObs(_nodePairs[i],x,epsObs,false);
  }
}
template <int N,int M>
void CollisionDetector<N,M>::updateBVH(const Vec& x,const Vec& x2,T eps,bool parallel) {
  if(parallel) {
    for(int l=0; l<(int)_layerOffsets.size()-1; l++) {
      if((int)_meshBVH.size()-l<=(1<<8)) {
        for(int i=l; i<(int)_meshBVH.size(); i++)
          updateBVH(x,x2,eps,i);
      } else {
        int beg=_layerOffsets[l],end=_layerOffsets[l+1];
        OMP_PARALLEL_FOR_
        for(int i=beg; i<end; i++)
          updateBVH(x,x2,eps,i);
      }
    }
  } else {
    for(int i=0; i<(int)_meshBVH.size(); i++)
      updateBVH(x,x2,eps,i);
  }
}
template <int N,int M>
void CollisionDetector<N,M>::updateBVH(const Vec& x,const Vec& x2,T eps,int i) {
  BBoxExact::Vec3T pos=BBoxExact::Vec3T::Zero();
  BBoxExact::Vec3T delta=BBoxExact::Vec3T::Zero();
  delta.template segment<N>(0).setConstant((BBoxExact::T)eps);
  auto& n=_meshBVH[i];
  if(n._cell>=0) {
    n._bb=BBoxExact();
    const auto& t=_meshIss[n._cell];
    for(int tvid=0; tvid<M; tvid++) {
      pos.template segment<N>(0)=x.template segment<N>(t[tvid]*N).template cast<BBoxExact::T>();
      n._bb.setUnion(pos);
      pos.template segment<N>(0)=x2.template segment<N>(t[tvid]*N).template cast<BBoxExact::T>();
      n._bb.setUnion(pos);
    }
    n._bb=n._bb.enlarged(delta);
  } else {
    auto& nl=_meshBVH[_meshBVH[i]._l];
    auto& nr=_meshBVH[_meshBVH[i]._r];
    n._bb=BBoxExact();
    n._bb.setUnion(nl._bb);
    n._bb.setUnion(nr._bb);
  }
}
//instance
template class CollisionDetector<2,2>;
template class CollisionDetector<3,3>;
template class CollisionDetector<2,6>;
template class CollisionDetector<3,6>;
}
