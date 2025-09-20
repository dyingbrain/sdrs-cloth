#include "Deformable.h"
#include "ADMM.h"
#include "GradientDescend.h"
#include "Newton.h"
#include "DirectNewton.h"
#include "CollisionSelf.h"
#include "CollisionObstacle.h"
#include "CollisionDetector.h"
#include "SpatialHash.h"
#include <Utils/Utils.h>
#include <Environment/EnvironmentUtils.h>

namespace PHYSICSMOTION {
template <int N>
void Deformable<N>::setAgents(const Vec& x,const Vec& xT,T r) {
  _x=_xL=_xLL=_x0=x;
  _xT=xT;
  _r=r;
}
template <int N>
void Deformable<N>::setARAP3D(const MeshExact& mesh,double feature_angle,double size) {
  int off=0,offV=_x0.size()/N;
  ASSERT_MSG(N==3,"setARAP3D can only be called for Deformable<3>!")
  _solver=NULL;
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  std::vector<Eigen::Matrix<int,4,1>> tss;
  for(const auto& v:mesh.vss())
    vss.push_back(v.template cast<double>());
  generateTetMesh(vss,mesh.iss(),"tet.mesh",feature_angle,size);
  readTetMesh(vss,tss,"tet.mesh");
  writeTetVTK(vss,tss,"tet.vtk");

  //extract surface
  extractSurfaceFromTetMesh(tss,iss);
  makeUniform(iss);
  if(volume(vss,iss)<0)
    makeInsideOut(iss);

  //build surface
  off=(int)_sss.size();
  _sss.resize(off+iss.size());
  for(int j=0; j<(int)iss.size(); j++)
    for(int i=0; i<N; i++)
      _sss[off+j][i]=iss[j][i]+offV;

  //build element
  off=(int)_tss.size();
  _tss.resize(off+tss.size());
  for(int j=0; j<(int)tss.size(); j++) {
    ARAPElement& e=_tss[off+j];
    e._stiffness=_springMaterial._k;
    for(int i=0; i<N+1; i++)
      e._indices[i]=tss[j][i]+offV;
  }

  //build vertex
  Vec x0;
  x0.setZero(vss.size()*N);
  for(int i=0; i<(int)vss.size(); i++)
    x0.template segment<N>(i*N)=vss[i].template segment<N>(0).template cast<T>();
  _x=_xL=_xLL=_x0=concat(_x0,x0);
}
template <int N>
void Deformable<N>::setARAP2D(const MeshExact& mesh) {
  int off=0,offV=_x0.size()/N;
  ASSERT_MSG(N==2,"setARAP2D can only be called for Deformable<2>!")
  _solver=NULL;

  //build element
  off=(int)_tss.size();
  _tss.resize(off+mesh.iss().size());
  for(int j=0; j<(int)mesh.iss().size(); j++) {
    ARAPElement& e=_tss[off+j];
    e._stiffness=_springMaterial._k;
    for(int i=0; i<N+1; i++)
      e._indices[i]=mesh.iss()[j][i]+offV;
  }

  //build vertex
  Vec x0;
  x0.setZero(mesh.vss().size()*N);
  for(int i=0; i<(int)mesh.vss().size(); i++)
    x0.template segment<N>(i*N)=mesh.vss()[i].template segment<N>(0).template cast<T>();
  _x=_xL=_xLL=_x0=concat(_x0,x0);
}
template <int N>
void Deformable<N>::setMassSpring(const MeshExact& mesh) {
  int off=0,offV=_x0.size()/N;
  _solver=NULL;

  //build edge
  for(const auto& t:mesh.iss())
    for(int d=0; d<N; d++) {
      Eigen::Matrix<int,2,1> e(t[d],t[(d+1)%3]);
      if(e[0]==-1 || e[1]==-1)
        continue;
      sort2(e[0],e[1]);
      if(e[0]!=e[1])
        _ess[(e.array()+offV).matrix()]=_springMaterial;
    }

  //build surface
  if(N==3) {
    off=(int)_sss.size();
    _sss.resize(off+mesh.iss().size());
    for(int j=0; j<(int)mesh.iss().size(); j++)
      for(int i=0; i<N; i++)
        _sss[off+j][i]=mesh.iss()[j][i]+offV;
  }

  Vec x0;
  x0.setZero(mesh.vss().size()*N);
  for(int i=0; i<(int)mesh.vss().size(); i++)
    x0.template segment<N>(i*N)=mesh.vss()[i].template segment<N>(0).template cast<T>();
  _x=_xL=_xLL=_x0=concat(_x0,x0);
}
template <int N>
bool Deformable<N>::setUAVTrajectory(const std::vector<VecNT>& waypoints,int numSegmentBeg,int numSegmentEnd) {
  bool succ=_UAVTraj.setNodesExhaustive(waypoints,_obs,_margin*2+_r,6,numSegmentBeg,numSegmentEnd);
  if(!succ) {
    _UAVTraj=ADMMUAVTrajectory<N>();
    return false;
  }
  _UAVTraj.buildSmoothness(2);
  _x=_xL=_xLL=_x0=_UAVTraj.getCP();
  return true;
}
template <int N>
void Deformable<N>::setObstacle(const MeshExact& mesh) {
  _obs=mesh;
}
template <int N>
void Deformable<N>::fix(int vid,T k,const VecNT* xStar) {
  _cons._fixes[vid]=std::make_pair(k,xStar?*xStar:_x0.template segment<N>(vid*N));
  _solver=NULL;
}
template <int N>
void Deformable<N>::fix(std::function<bool(const VecNT& a)> f,T k) {
  for(int i=0; i<_x0.size(); i+=N)
    if(f(_x0.template segment<N>(i)))
      fix(i/N,k);
}
template <int N>
int Deformable<N>::findClosestVid(const VecNT& pos) const {
  int id=-1;
  T minDist=0;
  for(int i=0; i<(int)_x0.size(); i+=N) {
    T dist=(_x0.template segment<N>(i)-pos).norm();
    if(id==-1 || dist<minDist) {
      minDist=dist;
      id=i/N;
    }
  }
  return id;
}
template <int N>
void Deformable<N>::addL1(int id,const VecNT& pos,T k) {
  int n=_x0.size()/N;

  Vec x;
  x.resize(_x0.size()+N);
  x.segment(0,_x0.size())=_x0;
  x.template segment<N>(_x0.size())=pos;
  _x=_xL=_xLL=_x0=x;

  fix(n,0);
  addL1({id,n},k);
}
template <int N>
void Deformable<N>::addL1(std::array<int,2> id,T k) {
  Eigen::Matrix<int,2,1> e(id[0],id[1]);
  sort2(e[0],e[1]);
  _l1ss[e]=k;
}
template <int N>
typename Deformable<N>::T Deformable<N>::dt() const {
  return _dt;
}
template <int N>
typename Deformable<N>::T Deformable<N>::r() const {
  return _r;
}
template <int N>
typename Deformable<N>::T Deformable<N>::actualR() const {
  return _r/(1+_margin);
}
template <int N>
void Deformable<N>::setDt(T dt) {
  _dt=dt;
  _solver=NULL;
}
template <int N>
void Deformable<N>::setR(T r) {
  _r=r;
  _solver=NULL;
}
template <int N>
void Deformable<N>::setK(T k) {
  _springMaterial._k=k;
  _solver=NULL;
}
template <int N>
void Deformable<N>::setB(T b) {
  _springMaterial._b=b;
  _solver=NULL;
}
template <int N>
void Deformable<N>::setCL(T cL) {
  _springMaterial._cL=cL;
  _solver=NULL;
}
template <int N>
void Deformable<N>::setCH(T cH) {
  _springMaterial._cH=cH;
  _solver=NULL;
}
template <int N>
void Deformable<N>::setMargin(T margin) {
  _margin=margin;
}
template <int N>
void Deformable<N>::setDamping(T damping) {
  _damping=damping;
}
template <int N>
void Deformable<N>::setCollCoef(T collCoef) {
  _collCoef=collCoef;
}
template <int N>
void Deformable<N>::setG(const VecNT& g) {
  _g=g;
}
template <int N>
void Deformable<N>::setUAVTrajResolution(int RES) {
  _trajSubd=_UAVTraj.subdivide(RES);
}
template <int N>
void Deformable<N>::setCB(std::function<bool()> cb) {
  _cb=cb;
  _hasCB=true;
}
template <int N>
void Deformable<N>::solve(const OptimizerParam& param) {
  if(param._type==OptimizerParam::ADMM && !std::dynamic_pointer_cast<ADMM<N>>(_solver))
    _solver=NULL;
  else if(param._type==OptimizerParam::GD && !std::dynamic_pointer_cast<GradientDescend<N>>(_solver))
    _solver=NULL;
  else if(param._type==OptimizerParam::NEWTON && !std::dynamic_pointer_cast<Newton<N>>(_solver))
    _solver=NULL;
  else if(param._type==OptimizerParam::DIRECT_NEWTON && !std::dynamic_pointer_cast<DirectNewton<N>>(_solver))
    _solver=NULL;
  if(!_solver) {
    if(param._type==OptimizerParam::ADMM)
      _solver.reset(new ADMM<N>);
    else if(param._type==OptimizerParam::GD)
      _solver.reset(new GradientDescend<N>);
    else if(param._type==OptimizerParam::NEWTON)
      _solver.reset(new Newton<N>);
    else if(param._type==OptimizerParam::DIRECT_NEWTON)
      _solver.reset(new DirectNewton<N>);
    else {
      ASSERT_MSG(false,"Unknown solver type!")
    }
    assemble();
  }
  //update state
  _xLL=_xL;
  _xL=_x;

  //gradient
  Vec g=Vec::Zero(_x0.size());
  MatNXT gV=_g*Vec::Ones(_x0.size()/N).transpose();
  //simulation
  if(_dt>0) {
    if(_r>0)
      g+=_mass*(-_xL-_xT*_dt)/(_dt*_dt);
    else g+=_mass*(_xLL-_xL*2)/(_dt*_dt);
    if(_damping>0)
      g+=_xL*_damping/_dt;
  } else if(_UAVTraj.numCP()>0) {
    //UAV planning
  } else if(_r>0)
    //navigation
    g-=_mass*_xT;
  //gravity
  g-=_mass*gV.reshaped();
  fix(NULL,&g);

  //solve
  _solver->addQuadratic(NULL,&g);
  _solver->x()=_x;
  _solver->setCB([&]() {
    if(_hasCB) {
      _x=_solver->x().segment(0,_x.size());
      return _cb();
    }
    return true;
  });
  _solver->optimize(param);
  _x=_solver->x();
}
template <int N>
const MeshExact& Deformable<N>::buildSelfCollisionMesh() {
  std::vector<Eigen::Matrix<double,N,1>> vss;
  std::vector<Eigen::Matrix<int,N,1>> iss;
  //insert vertices from x0
  vss.resize(_x0.size()/N);
  for(int i=0; i<_x0.size()/N; i++)
    vss[i]=_x0.template segment<N>(i*N).template cast<double>();
  //insert indices of different cases
  if(!_sss.empty()) {
    //if surface is labeled, use it
    iss.resize(_sss.size());
    for(int i=0; i<(int)_sss.size(); i++)
      iss[i]=_sss[i];
  } else if(!_ess.empty()) {
    //if edge exists, use it
    int off=0;
    iss.assign(_ess.size(),Eigen::Matrix<int,N,1>::Zero());
    for(const auto& spr:_ess)
      iss[off++].template segment<2>(0)=spr.first;
  } else if(!_tss.empty()) {
    //if element exists, extract surface and use it
    EdgeMap emap;
    ASSERT(N==2)
    for(const auto& e:_tss) {
      const auto& t=e._indices;
      for(int d=0; d<N+1; d++) {
        Eigen::Matrix<int,2,1> ev(t[d],t[(d+1)%3]);
        sort2(ev[0],ev[1]);
        if(ev[0]!=ev[1]) {
          if(emap.find(ev)==emap.end())
            emap[ev]=1;
          else emap[ev]++;
        }
      }
    }
    Eigen::Matrix<int,N,1> i=Eigen::Matrix<int,N,1>::Zero();
    for(const auto& e:emap)
      if(e.second==1) {
        i.template segment<2>(0)=e.first;
        iss.push_back(i);
      }
  }
  //rebuild
  _self=MeshExact(vss,iss,true);
  return _self;
}
//helper
template <int N>
void Deformable<N>::assemble() {
  ASSERT(_solver);

  //hessian
  buildMassMatrix(_m);

  //assemble hessian
  SMatT H;
  STrips HTrips;
  fix(&HTrips,NULL);
  H.resize(_x0.size(),_x0.size());
  H.setFromTriplets(HTrips.begin(),HTrips.end());
  if(_dt>0) {
    //simulation
    H+=_mass/(_dt*_dt);
    if(_damping>0)
      H+=_mass*_damping/_dt;
  } else if(_UAVTraj.numCP()>0)
    //UAV planning
    H+=_UAVTraj.getSmoothH();
  else if(_r>0)
    //navigation
    H+=_mass;

  //setup quadratic term
  _solver->addQuadratic(&H,NULL);

  //setup mass spring 1-ring term
  Ring ring1;
  buildRing1(ring1);
  for(const auto& ess:ring1) {
    int v1=ess.first;
    const Spring& material=ess.second.second;
    for(const auto& v2:ess.second.first) {
      if(v1>=v2)
        continue;
      T l=(_x0.template segment<N>(v1*N)-_x0.template segment<N>(v2*N)).norm();
      if(material._k>0)
        _solver->addSpring({v1,v2},material._k,l,l*material._cL,l*material._cH);
    }
  }

  //setup mass spring 2-ring term
  Ring ring2;
  buildRing2(ring1,ring2);
  for(const auto& ess:ring2) {
    int v1=ess.first;
    const Spring& material=ess.second.second;
    for(const auto& v2:ess.second.first) {
      if(v1>=v2)
        continue;
      T l=(_x0.template segment<N>(v1*N)-_x0.template segment<N>(v2*N)).norm();
      if(material._b>0)
        _solver->addSpring({v1,v2},material._b,l,l*material._cL,l*material._cH);
    }
  }

  //setup ARAP term
  for(const auto& e:_tss) {
    MatNT F0;
    const auto& t=e._indices;
    for(int d=0; d<N; d++)
      F0.col(d)=_x0.template segment<N>(t[d+1]*N)-_x0.template segment<N>(t[0]*N);
    if(abs(F0.determinant())<1e-4f) {
      std::cout << "Detected degenerate tet " << t.transpose() << " vol=" << F0.determinant();
      ASSERT(false)
    }
    std::array<int,N+1> a;
    for(int i=0; i<=N; i++)
      a[i]=t[i];
    _solver->addARAPElement(a,e._stiffness,F0);
  }

  //setup l1
  for(const auto& l1:_l1ss)
    _solver->addSmoothL1({l1.first[0],l1.first[1]},l1.second);

  //setup collision
  if(_collCoef>0) {
    if(_UAVTraj.numCP()>0) {
      //UAV trajectory planning
      int nCP=_UAVTraj.numCP();
      _solver->insertTerm(std::shared_ptr<CollisionObstacle<N,6,N>>(new CollisionObstacle<N,6,N>(_r,_margin,_collCoef)));
      _solver->insertCollisionDetector(std::shared_ptr<CollisionDetector<N,6>>(new CollisionDetector<N,6>(nCP,_obs)));
      _cons._general=_UAVTraj.getCPCons();
    } else if(_r>0) {
      //agent navigation
      _solver->insertTerm(std::shared_ptr<CollisionSelf<N,1>>(new CollisionSelf<N,1>(_r,_margin,_collCoef)));
      _solver->insertTerm(std::shared_ptr<CollisionObstacle<N,1,N>>(new CollisionObstacle<N,1,N>(_r,_margin,_collCoef)));
      _solver->insertCollisionDetector(std::shared_ptr<CollisionDetector<N,1>>(new CollisionDetector<N,1>(_obs)));
    } else {
      //other
      buildSelfCollisionMesh();
      _solver->insertTerm(std::shared_ptr<CollisionSelf<N,N>>(new CollisionSelf<N,N>(_r,_margin,_collCoef)));
      _solver->insertTerm(std::shared_ptr<CollisionObstacle<N,N,N>>(new CollisionObstacle<N,N,N>(_r,_margin,_collCoef)));
      _solver->insertCollisionDetector(std::shared_ptr<CollisionDetector<N,N>>(new CollisionDetector<N,N>(_self,_obs)));
    }
  }

  //assemble solver
  _solver->assemble(_x0.size(),_cons);
}
template <int N>
void Deformable<N>::fix(STrips* HTrips,Vec* g) const {
  for(const auto& f:_cons._fixes) {
    int vid=f.first;
    T k=f.second.first;
    if(k<=0)
      continue; //use hard constraint
    const VecNT& xStar=f.second.second;
    if(HTrips)
      addBlockId(*HTrips,vid*N,vid*N,N,k);
    if(g)
      g->template segment<N>(vid*N)-=xStar*k;
  }
}
template <int N>
void Deformable<N>::buildRing1(Ring& ring1) const {
  Eigen::Matrix<int,2,1> e;
  for(const auto& spr:_ess) {
    int v0=spr.first[0];
    int v1=spr.first[1];
    ring1[v0].first.insert(v1);
    ring1[v1].first.insert(v0);
    ring1[v0].second=spr.second;
  }
}
template <int N>
void Deformable<N>::buildRing2(const Ring& ring1,Ring& ring2) const {
  for(const auto& ess:ring1) {
    int v1=ess.first;
    std::unordered_set<int> vss;
    //insert all neighbors of 2-ring vertices
    for(int v2:ess.second.first) {
      auto it=ring1.find(v2);
      if(it!=ring1.end())
        vss.insert(it->second.first.begin(),it->second.first.end());
    }
    //erase the center vertex and 1-ring vertices
    vss.erase(v1);
    for(int v2:ess.second.first)
      if(vss.find(v2)!=vss.end())
        vss.erase(v2);
    //insert into ring2
    ring2[v1]=std::make_pair(vss,ess.second.second);
  }
}
template <int N>
void Deformable<N>::buildMassMatrix(T coef) {
  STrips HTrips;
  _mass.resize(_x0.size(),_x0.size());
  if(_tss.empty())
    buildPointMassMatrix(HTrips,coef);
  else buildFEMMassMatrix(HTrips,coef);
  _mass.setFromTriplets(HTrips.begin(),HTrips.end());
}
template <int N>
void Deformable<N>::buildPointMassMatrix(STrips& HTrips,T coef) const {
  for(int i=0; i<_x0.size(); i+=N)
    addBlockId(HTrips,i,i,N,coef);
}
template <int N>
void Deformable<N>::buildFEMMassMatrix(STrips& HTrips,T coef) const {
  Eigen::Matrix<T,N+1,N+1> H;
  H.setOnes();
  H.diagonal()*=2;
  for(const auto& e:_tss) {
    const auto& t=e._indices;
    //compute volume
    MatNT F;
    for(int r=0; r<N; r++)
      F.col(r)=_x0.template segment<N>(t[r+1]*N)-_x0.template segment<N>(t[0]*N);
    //HCoef*H.row(0).sum()=area/(N+1);
    T HCoef=abs(F.determinant())/(N+1)/H.row(0).sum();
    //insert entries
    for(int r=0; r<N+1; r++)
      for(int c=0; c<N+1; c++)
        addBlockId<FLOAT>(HTrips,t[r]*N,t[c]*N,N,coef*HCoef*H(r,c));
  }
}
//instance
template class Deformable<2>;
template class Deformable<3>;
}
