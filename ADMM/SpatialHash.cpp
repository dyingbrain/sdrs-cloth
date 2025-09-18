#include "SpatialHash.h"
#include "CollisionObstacle.h"
#include "CollisionSelf.h"
#include "Narrow.h"
#ifdef FLOAT
#undef FLOAT
#include "pradsort.hpp"
#endif
#include <Utils/Utils.h>

namespace PHYSICSMOTION {
void radixSort(int* val,int* key,int N) {
  prsort::pradsort(val,key,N,8,NULL);
}
//CollisionDetector
template <int N>
CollisionDetector<N,1>::CollisionDetector(const MeshExact& obs):_obs(obs) {}
template <int N>
void CollisionDetector<N,1>::generateBroadBF(const Vec& x,const Vec& x2,T eps,T epsObs,bool self,bool obs) {
  _pssSelf.clear();
  _pssObs.clear();
  int n=x.size()/N;
  //self
  if(self)
    for(int l=0; l<n; l++) {
      BBoxExact bbl=getBB(x,x2,l,eps);
      for(int r=l+1; r<n; r++) {
        //intersection test
        BBoxExact bbr=getBB(x,x2,r,eps);
        if(!bbl.intersect(bbr))
          continue;
        //add entry
        _pssSelf.push_back(combine(l,r));
      }
    }
  //obs
  if(obs && !_obs.iss().empty())
    for(int l=0; l<n; l++) {
      BBoxExact bbl=getBB(x,x2,l,eps+epsObs);
      for(int r=0; r<(int)_obs.iss().size(); r++) {
        //intersection test
        if(!bbl.intersect(_obs.getBVH()[r]._bb))
          continue;
        //add entry
        _pssObs.push_back(combine(l,r));
      }
    }
}
template <int N>
void CollisionDetector<N,1>::generateBroad(const Vec& x,const Vec& x2,T eps,T epsObs,bool self,bool obs) {
  _pssSelf.clear();
  _pssObs.clear();
  int n=x.size()/N;
  buildHash(x,x2,eps);
  //self
  if(self) {
    OMP_PARALLEL_FOR_
    for(int vid=0; vid<n; vid++) {
      T searchRange=_nodes[vid]._radius+_R;
      Eigen::Matrix<int,N,1> L=hash(_nodes[vid]._ctr-VecNT::Constant(searchRange)).cwiseMax(Eigen::Matrix<int,N,1>::Zero());
      Eigen::Matrix<int,N,1> U=hash(_nodes[vid]._ctr+VecNT::Constant(searchRange)).cwiseMin(_nrCell-Eigen::Matrix<int,N,1>::Ones());
      if(N==2) {
        //2D
        for(int x=L[0],offX=L.dot(_stride); x<=U[0]; x++,offX+=_stride[0])
          for(int y=L[1],offY=offX; y<=U[1]; y++,offY+=_stride[1])
            for(int start=_starts[offY],end=_ends[offY]; start<end; start++) {
              int vid0=vid,vid1=_offsets[start];
              if(vid0<vid1 && _nodes[vid0]._bb.intersect(_nodes[vid1]._bb))
                _pssSelf.push_back(combine(vid0,vid1));
            }
      } else {
        //3D
        for(int x=L[0],offX=L.dot(_stride); x<=U[0]; x++,offX+=_stride[0])
          for(int y=L[1],offY=offX; y<=U[1]; y++,offY+=_stride[1])
            for(int z=L[2],offZ=offY; z<=U[2]; z++,offZ+=_stride[2])
              for(int start=_starts[offZ],end=_ends[offZ]; start<end; start++) {
                int vid0=vid,vid1=_offsets[start];
                if(vid0<vid1 && _nodes[vid0]._bb.intersect(_nodes[vid1]._bb))
                  _pssSelf.push_back(combine(vid0,vid1));
              }
      }
    }
  }
  //obs
  if(obs && !_obs.iss().empty()) {
    BBoxExact::Vec3T pos=BBoxExact::Vec3T::Zero();
    pos.template segment<N>(0).setConstant((BBoxExact::T)epsObs);
    OMP_PARALLEL_FOR_
    for(int vid=0; vid<n; vid++)
      generateObs(vid,_nodes[vid]._bb.enlarged(pos));
  }
}
template <int N>
Eigen::Matrix<int,2,1> CollisionDetector<N,1>::generate(const Vec& x,const Vec& x2,Optimizer<N>& opt,bool noNeighbhor,bool,int restoreSlot) {
  std::shared_ptr<CollisionSelf<N,1>> self=opt.template findTerm<CollisionSelf<N,1>>();
  std::shared_ptr<CollisionObstacle<N,1>> obs=opt.template findTerm<CollisionObstacle<N,1>>();
  //broad phase
  T eps=self?self->eps()/2:0,epsObs=std::max<T>(obs->eps()-eps,0);
  generateBroad(x,x2,eps,epsObs,self!=NULL,obs!=NULL);
  //narrow phase
  Eigen::Matrix<int,2,1> ret=Eigen::Matrix<int,2,1>::Zero();
  Narrow<N,1> narrow(_obs);
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
template <int N>
Eigen::Matrix<int,2,1> CollisionDetector<N,1>::removeByDistance(const Vec& x,Optimizer<N>& opt,T margin) {
  std::shared_ptr<CollisionSelf<N,1>> self=opt.template findTerm<CollisionSelf<N,1>>();
  std::shared_ptr<CollisionObstacle<N,1>> obs=opt.template findTerm<CollisionObstacle<N,1>>();
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
template <int N>
Eigen::Matrix<int,2,1> CollisionDetector<N,1>::removeByEnergy(const Vec& x,Optimizer<N>& opt,T thres) {
  std::shared_ptr<CollisionSelf<N,1>> self=opt.template findTerm<CollisionSelf<N,1>>();
  std::shared_ptr<CollisionObstacle<N,1>> obs=opt.template findTerm<CollisionObstacle<N,1>>();
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
template <int N>
std::string CollisionDetector<N,1>::info(const Optimizer<N>& opt) const {
  std::shared_ptr<CollisionSelf<N,1>> self=opt.template findTerm<CollisionSelf<N,1>>();
  std::shared_ptr<CollisionObstacle<N,1>> obs=opt.template findTerm<CollisionObstacle<N,1>>();
  std::ostringstream oss;
  if(self)
    oss << "#SelfColl=" << self->n();
  if(self && obs)
    oss << " ";
  if(obs)
    oss << "#ObsColl=" << obs->n();
  return oss.str();
}
template <int N>
typename CollisionDetector<N,1>::T CollisionDetector<N,1>::getEps(Optimizer<N>& opt) const {
  std::shared_ptr<CollisionSelf<N,1>> self=opt.template findTerm<CollisionSelf<N,1>>();
  std::shared_ptr<CollisionObstacle<N,1>> obs=opt.template findTerm<CollisionObstacle<N,1>>();
  T eps=self?self->eps()/2:0,epsObs=std::max<T>(obs->eps()-eps,0);
  return std::max<T>(eps,epsObs);
}
//getter
template <int N>
void CollisionDetector<N,1>::extractPos(const MeshExact& m,Vec& x) const {
  x.resize((int)m.vss().size()*N);
  for(int i=0; i<(int)m.vss().size(); i++)
    x.template segment<N>(i*N)=m.vss()[i].template segment<N>(0).template cast<T>();
}
template <int N>
const std::vector<typename CollisionDetector<N,1>::ID>& CollisionDetector<N,1>::selfCollisions() const {
  return _pssSelf.getVector();
}
template <int N>
const std::vector<typename CollisionDetector<N,1>::ID>& CollisionDetector<N,1>::obsCollisions() const {
  return _pssObs.getVector();
}
template <int N>
Eigen::Matrix<int,2,1> CollisionDetector<N,1>::selfCollisionId(ID c) const {
  int l,r;
  separate(c,l,r);
  return Eigen::Matrix<int,2,1>(l,r);
}
template <int N>
Eigen::Matrix<int,1,1> CollisionDetector<N,1>::obsCollisionId(ID c) const {
  int l,r;
  separate(c,l,r);
  return Eigen::Matrix<int,1,1>(l);
}
template <int N>
typename CollisionDetector<N,1>::VecNT CollisionDetector<N,1>::obsCollisionPos(ID c,int j) const {
  int l,r;
  separate(c,l,r);
  return _obs.vss()[_obs.iss()[r][j]].template segment<N>(0).template cast<T>();
}
//ID
template <int N>
void CollisionDetector<N,1>::separate(ID pair,int& l,int& r) {
  static const ID mask=((ID)1<<(ID)32)-1;
  l=(int)(pair>>32);
  r=(int)(pair&mask);
}
template <int N>
typename CollisionDetector<N,1>::ID CollisionDetector<N,1>::combine(ID l,ID r) {
  return (l<<32)+r;
}
//helper
template <int N>
int CollisionDetector<N,1>::hashOff(const VecNT& pt) const {
  return hash(pt).dot(_stride);
}
template <int N>
Eigen::Matrix<int,N,1> CollisionDetector<N,1>::hash(const VecNT& pt) const {
  return ((pt-_bb.minCorner().template segment<N>(0).template cast<T>())*_invR).array().floor().matrix().template cast<int>().cwiseMin(_nrCell-Eigen::Matrix<int,N,1>::Ones());
}
template <int N>
BBoxExact CollisionDetector<N,1>::getBB(const Vec& x,const Vec& x2,int id,T eps) const {
  BBoxExact::Vec3T pos=BBoxExact::Vec3T::Zero();
  BBoxExact bb;
  pos.template segment<N>(0)=x.template segment<N>(id*N).template cast<BBoxExact::T>();
  bb.setUnion(pos);
  pos.template segment<N>(0)=x2.template segment<N>(id*N).template cast<BBoxExact::T>();
  bb.setUnion(pos);
  //eps
  if(eps>0) {
    pos.template segment<N>(0).setConstant((BBoxExact::T)eps);
    return bb.enlarged(pos);
  } else return bb;
}
template <int N>
typename CollisionDetector<N,1>::SpatialHashNode CollisionDetector<N,1>::getNode(const Vec& x,const Vec& x2,int id,T eps) const {
  SpatialHashNode node;
  node._bb=getBB(x,x2,id,0);
  node._ctr=(node._bb.minCorner()+node._bb.maxCorner()).template segment<N>(0).template cast<T>()/2;
  node._radius=(T)(node._bb.maxCorner()-node._bb.minCorner()).norm()/2+eps;
  //enlarge
  BBoxExact::Vec3T pos=BBoxExact::Vec3T::Zero();
  pos.template segment<N>(0).setConstant((BBoxExact::T)eps);
  node._bb=node._bb.enlarged(pos);
  return node;
}
template <int N>
void CollisionDetector<N,1>::reduce(std::function<void(SpatialHashNode&,SpatialHashNode&)> op) {
  for(int off=2,offHalf=1; offHalf<(int)_nodes.size(); off<<=1,offHalf<<=1) {
    OMP_PARALLEL_FOR_
    for(int i=0; i<(int)_nodes.size(); i+=off)
      if(i+offHalf<(int)_nodes.size())
        op(_nodes[i],_nodes[i+offHalf]);
  }
}
template <int N>
void CollisionDetector<N,1>::generateObs(int vid,const BBoxExact& bb) {
  //serial computation
  std::stack<int> ss;
  ss.push((int)_obs.getBVH().size()-1);
  while(!ss.empty()) {
    int r=ss.top();
    ss.pop();
    const auto& n=_obs.getBVH()[r];
    if(!n._bb.intersect(bb))
      continue;
    else if(n._cell>=0)
      _pssObs.push_back(combine(vid,r));
    else {
      ss.push(n._l);
      ss.push(n._r);
    }
  }
}
template <int N>
void CollisionDetector<N,1>::buildHash(const Vec& x,const Vec& x2,T eps) {
  int n=x.size()/N;
  _nodes.resize(n);
  OMP_PARALLEL_FOR_
  for(int id=0; id<n; id++)
    _nodes[id]=getNode(x,x2,id,eps);
  //find range
  _nodesBkg=_nodes;
  reduce([&](SpatialHashNode& a,SpatialHashNode& b) {
    a._radius=std::max(a._radius,b._radius);
    a._bb.setUnion(b._bb);
  });
  _R=_nodes[0]._radius;
  _bb=_nodes[0]._bb;
  _nodes.swap(_nodesBkg);
  //compute parameter
  _invR=1/_R;
  _nrCell.array()=((_bb.maxCorner()-_bb.minCorner()).template segment<N>(0).template cast<T>()*_invR).cwiseMax(VecNT::Ones()).array().ceil().template cast<int>();
  _stride.setOnes();
  for(int i=0; i<N; i++)
    for(int j=i+1; j<N; j++)
      _stride[i]*=_nrCell[j];
  //sort
  _indices.resize(n);
  _offsetsInv.resize(n);
  _offsets.resize(n);
  OMP_PARALLEL_FOR_
  for(int i=0; i<n; i++) {
    _indices[i]=hashOff(_nodes[i]._ctr);
    _offsetsInv[i]=i;
  }
  radixSort(_indices.data(),_offsetsInv.data(),n);
  //mark start/end
  _starts.assign(_nrCell.prod(),-1);
  _ends.assign(_nrCell.prod(),-1);
  OMP_PARALLEL_FOR_
  for(int i=0; i<n; i++) {
    _offsets[_offsetsInv[i]]=i;
    if(i==0)
      _starts[_indices[i]]=i;
    else if(_indices[i]!=_indices[i-1])
      _ends[_indices[i-1]]=_starts[_indices[i]]=i;
    if(i==n-1)
      _ends[_indices[i]]=i+1;
  }
}
//instance
template class CollisionDetector<2,1>;
template class CollisionDetector<3,1>;
}
