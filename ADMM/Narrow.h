#ifndef NARROW_H
#define NARROW_H

#include "CollisionDetector.h"
#include "SmallShape.h"
#include <Environment/GJK.h>
#include <Environment/DistanceFunction.h>
#include <Utils/Utils.h>
//use CCD library
//#define USE_CCD
#define USE_SIMPLE_TRIANGLE_CCD
//#define USE_EXACT_CCD
#ifdef USE_EXACT_CCD
#include "CCD/ccd.hpp"
#else
#include "CTCD/CTCD.h"
#endif

namespace PHYSICSMOTION {
struct EdgeEdgeHash {
  size_t operator()(const Eigen::Matrix<int,4,1>& key) const {
    size_t seed=0;
    std::hash<int> h;
    hash_combine(seed,h(key[0]));
    hash_combine(seed,h(key[1]));
    hash_combine(seed,h(key[2]));
    hash_combine(seed,h(key[3]));
    return seed;
  }
  bool operator()(const Eigen::Matrix<int,4,1>& a, const Eigen::Matrix<int,4,1>& b) const {
    for(int i=0; i<4; i++)
      if(a[i]<b[i])
        return true;
      else if(a[i]>b[i])
        return false;
    return false;
  }
};
static bool edgeEdgeCCD(const Eigen::Vector3d& q0s,const Eigen::Vector3d& p0s,
                        const Eigen::Vector3d& q1s,const Eigen::Vector3d& p1s,
                        const Eigen::Vector3d& q0e,const Eigen::Vector3d& p0e,
                        const Eigen::Vector3d& q1e,const Eigen::Vector3d& p1e,double eps=1e-4f) {
#ifdef USE_EXACT_CCD
  return ccd::edgeEdgeCCD(q0s,p0s,q1s,p1s,
                          q0e,p0e,q1e,p1e);
#else
  double t;
  return CTCD::edgeEdgeCTCD(q0s,p0s,q1s,p1s,
                            q0e,p0e,q1e,p1e,(double)eps,t);
#endif
}
static bool vertexFaceCCD(const Eigen::Vector3d& q0s,const Eigen::Vector3d& q1s,
                          const Eigen::Vector3d& q2s,const Eigen::Vector3d& q3s,
                          const Eigen::Vector3d& q0e,const Eigen::Vector3d& q1e,
                          const Eigen::Vector3d& q2e,const Eigen::Vector3d& q3e,double eps=1e-4f) {
#ifdef USE_EXACT_CCD
  return ccd::vertexFaceCCD(q0s,q1s,q2s,q3s,
                            q0e,q1e,q2e,q3e);
#else
  double t;
  return CTCD::vertexFaceCTCD(q0s,q1s,q2s,q3s,
                              q0e,q1e,q2e,q3e,(double)eps,t);
#endif
}
template <int N,int M>
class NarrowBase {
 public:
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  typedef Eigen::Matrix<T,N,1> VecNT;
  typedef unsigned long long ID;
  virtual ~NarrowBase() = default;
  int self(const Vec& x,const Vec& x2,ParallelVector<ID>& IDs,const std::unordered_map<ID,int>& existing,T eps) const {
    int ret=0;
    auto& v=IDs.getVector();
    OMP_PARALLEL_FOR_
    for(int i=0; i<(int)v.size(); i++) {
      if(existing.find(v[i])!=existing.end())
        v[i]=CollisionDetector<N,M>::INVALID_ID;
      else self(x,x2,v[i],eps);
      if(v[i]!=CollisionDetector<N,M>::INVALID_ID)
        OMP_ATOMIC_
        ret++;
    }
    return ret;
  }
  int obs(const Vec& x,const Vec& x2,ParallelVector<ID>& IDs,const std::unordered_map<ID,int>& existing,T eps) const {
    int ret=0;
    auto& v=IDs.getVector();
    OMP_PARALLEL_FOR_
    for(int i=0; i<(int)v.size(); i++) {
      if(existing.find(v[i])!=existing.end())
        v[i]=CollisionDetector<N,M>::INVALID_ID;
      else obs(x,x2,v[i],eps);
      if(v[i]!=CollisionDetector<N,M>::INVALID_ID)
        ret++;
    }
    return ret;
  }
  virtual void self(const Vec& x,const Vec& x2,ID& nodePair,T eps) const=0;
  virtual void obs(const Vec& x,const Vec& x2,ID& nodePair,T eps) const=0;
};
template <int N,int M>
class Narrow : public NarrowBase<N,M> {
 public:
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  typedef Eigen::Matrix<T,N,1> VecNT;
  typedef std::vector<Eigen::Matrix<int,M,1>> Iss;
  typedef unsigned long long ID;
  using NarrowBase<N,M>::self;
  using NarrowBase<N,M>::obs;
  Narrow(const Iss& iss,const MeshExact& obs):_meshIss(iss),_obs(obs) {}
  void self(const Vec& x,const Vec& x2,ID& nodePair,T eps) const override {
    //getID
    int l,r;
    CollisionDetector<N,M>::separate(nodePair,l,r);
    const Eigen::Matrix<int,M,1>& ml=_meshIss[l];
    const Eigen::Matrix<int,M,1>& mr=_meshIss[r];
    GJK::Vec3T pAL,pBL;
    bool intersect;
    GJK::T dist;
    //CCD
    Eigen::Matrix<T,N,M*2> a;
    Eigen::Matrix<T,N,M*2> b;
    for(int i=0; i<M; i++) {
      a.col(i)=x.template segment<N>(ml[i]*N);
      b.col(i)=x.template segment<N>(mr[i]*N);
      a.col(i+M)=x2.template segment<N>(ml[i]*N);
      b.col(i+M)=x2.template segment<N>(mr[i]*N);
    }
    //GJK
    SmallShape<N,M*2> A(a);
    SmallShape<N,M*2> B(b);
    dist=GJK::runGJK(A,B,
                     GJK::Mat3X4T::Identity(),
                     GJK::Mat3X4T::Identity(),
                     pAL,pBL,&intersect);
    if((T)dist>eps)
      nodePair=CollisionDetector<N,M>::INVALID_ID;
  }
  void obs(const Vec& x,const Vec& x2,ID& nodePair,T eps) const override {
    //getID
    int l,r;
    CollisionDetector<N,M>::separate(nodePair,l,r);
    const Eigen::Matrix<int,M,1>& ml=_meshIss[l];
    const Eigen::Matrix<int,3,1>& mr=_obs.iss()[r];
    GJK::Vec3T pAL,pBL;
    bool intersect;
    GJK::T dist;
    //CCD
    Eigen::Matrix<T,N,M*2> a;
    Eigen::Matrix<T,N,N> b;
    for(int i=0; i<M; i++) {
      a.col(i)=x.template segment<N>(ml[i]*N);
      a.col(i+M)=x2.template segment<N>(ml[i]*N);
    }
    for(int i=0; i<N; i++)
      b.col(i)=_obs.vss()[mr[i]].template segment<N>(0).template cast<T>();
    //GJK
    SmallShape<N,M*2> A(a);
    SmallShape<N,N> B(b);
    dist=GJK::runGJK(A,B,
                     GJK::Mat3X4T::Identity(),
                     GJK::Mat3X4T::Identity(),
                     pAL,pBL,&intersect);
    if((T)dist>eps)
      nodePair=CollisionDetector<N,M>::INVALID_ID;
  }
 protected:
  const Iss& _meshIss;
  const MeshExact& _obs;
};
template <int N>
class Narrow<N,1> : public NarrowBase<N,1> {
 public:
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  typedef Eigen::Matrix<T,N,1> VecNT;
  typedef unsigned long long ID;
  using NarrowBase<N,1>::self;
  using NarrowBase<N,1>::obs;
  Narrow(const MeshExact& obs):_obs(obs) {}
  void self(const Vec& x,const Vec& x2,ID& nodePair,T eps) const override {
    //getID
    int l,r;
    CollisionDetector<N,1>::separate(nodePair,l,r);
    GJK::Vec3T pAL,pBL;
    bool intersect;
    GJK::T dist;
    //CCD
    Eigen::Matrix<T,N,2> a;
    Eigen::Matrix<T,N,2> b;
    a.col(0)=x.template segment<N>(l*N);
    b.col(0)=x.template segment<N>(r*N);
    a.col(1)=x2.template segment<N>(l*N);
    b.col(1)=x2.template segment<N>(r*N);
    //GJK
    SmallShape<N,2> A(a);
    SmallShape<N,2> B(b);
    dist=GJK::runGJK(A,B,
                     GJK::Mat3X4T::Identity(),
                     GJK::Mat3X4T::Identity(),
                     pAL,pBL,&intersect);
    if((T)dist>eps)
      nodePair=CollisionDetector<N,1>::INVALID_ID;
  }
  void obs(const Vec& x,const Vec& x2,ID& nodePair,T eps) const override {
    //getID
    int l,r;
    CollisionDetector<N,1>::separate(nodePair,l,r);
    const Eigen::Matrix<int,3,1>& mr=_obs.iss()[r];
    GJK::Vec3T pAL,pBL;
    bool intersect;
    GJK::T dist;
    //CCD
    Eigen::Matrix<T,N,2> a;
    Eigen::Matrix<T,N,N> b;
    a.col(0)=x.template segment<N>(l*N);
    a.col(1)=x2.template segment<N>(l*N);
    for(int i=0; i<N; i++)
      b.col(i)=_obs.vss()[mr[i]].template segment<N>(0).template cast<T>();
    //GJK
    SmallShape<N,2> A(a);
    SmallShape<N,N> B(b);
    dist=GJK::runGJK(A,B,
                     GJK::Mat3X4T::Identity(),
                     GJK::Mat3X4T::Identity(),
                     pAL,pBL,&intersect);
    if((T)dist>eps)
      nodePair=CollisionDetector<N,1>::INVALID_ID;
  }
 protected:
  const MeshExact& _obs;
};
#ifdef USE_CCD
template <>
class Narrow<2,2> : public NarrowBase<2,2> {
 public:
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  typedef Eigen::Matrix<T,2,1> VecNT;
  typedef std::vector<Eigen::Matrix<int,2,1>> Iss;
  typedef unsigned long long ID;
  using NarrowBase<2,2>::self;
  using NarrowBase<2,2>::obs;
  Narrow(const Iss& iss,const MeshExact& obs):_meshIss(iss),_obs(obs) {}
  void self(const Vec& x,const Vec& x2,ID& nodePair,T eps) const override {
    //getID
    int l,r;
    CollisionDetector<2,2>::separate(nodePair,l,r);
    const Eigen::Matrix<int,2,1>& ml=_meshIss[l];
    const Eigen::Matrix<int,2,1>& mr=_meshIss[r];
    Eigen::Matrix<T,2,2> a;
    Eigen::Matrix<T,2,2> b;
    GJK::Vec3T pAL,pBL;
    bool intersect;
    GJK::T dist;
    //CCD
    if(edgeEdgeCCD(extract(x,ml[0]),extract(x,ml[1]),
                   extract(x,mr[0]),extract(x,mr[1]),
                   extract(x2,ml[0]),extract(x2,ml[1]),
                   extract(x2,mr[0]),extract(x2,mr[1]))) {
      //need to insert anyway
    } else {
      //only check x2
      for(int i=0; i<2; i++) {
        a.col(i)=x2.template segment<2>(ml[i]*2);
        b.col(i)=x2.template segment<2>(mr[i]*2);
      }
      //GJK
      SmallShape<2,2> A(a);
      SmallShape<2,2> B(b);
      dist=GJK::runGJK(A,B,
                       GJK::Mat3X4T::Identity(),
                       GJK::Mat3X4T::Identity(),
                       pAL,pBL,&intersect);
      if((T)dist>eps)
        nodePair=CollisionDetector<2,2>::INVALID_ID;
    }
  }
  void obs(const Vec& x,const Vec& x2,ID& nodePair,T eps) const override {
    //getID
    int l,r;
    CollisionDetector<2,2>::separate(nodePair,l,r);
    const Eigen::Matrix<int,2,1>& ml=_meshIss[l];
    const Eigen::Matrix<int,3,1>& mr=_obs.iss()[r];
    Eigen::Matrix<T,2,2> a;
    Eigen::Matrix<T,2,2> b;
    GJK::Vec3T pAL,pBL;
    bool intersect;
    GJK::T dist;
    //CCD
    if(edgeEdgeCCD(extract(x,ml[0]),extract(x,ml[1]),
                   extract(_obs,mr[0]),extract(_obs,mr[1]),
                   extract(x2,ml[0]),extract(x2,ml[1]),
                   extract(_obs,mr[0]),extract(_obs,mr[1]))) {
      //need to insert anyway
    } else {
      //only check x2
      for(int i=0; i<2; i++) {
        a.col(i)=x2.template segment<2>(ml[i]*2);
        b.col(i)=_obs.vss()[mr[i]].template segment<2>(0).template cast<T>();
      }
      //GJK
      SmallShape<2,2> A(a);
      SmallShape<2,2> B(b);
      dist=GJK::runGJK(A,B,
                       GJK::Mat3X4T::Identity(),
                       GJK::Mat3X4T::Identity(),
                       pAL,pBL,&intersect);
      if((T)dist>eps)
        nodePair=CollisionDetector<2,2>::INVALID_ID;
    }
  }
 private:
  static Eigen::Vector3d extract(const Vec& x,int off) {
    off*=2;
    return Eigen::Vector3d((double)x[off],(double)x[off+1],0);
  }
  static Eigen::Vector3d extract(const MeshExact& obs,int off) {
    return obs.vss()[off].template cast<double>();
  }
  const Iss& _meshIss;
  const MeshExact& _obs;
};
#ifdef USE_SIMPLE_TRIANGLE_CCD
template <>
class Narrow<3,3> : public NarrowBase<3,3> {
 public:
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  typedef Eigen::Matrix<T,3,1> VecNT;
  typedef std::vector<Eigen::Matrix<int,3,1>> Iss;
  typedef unsigned long long ID;
  using NarrowBase<3,3>::self;
  using NarrowBase<3,3>::obs;
  Narrow(const Iss& iss,const MeshExact& obs):_meshIss(iss),_obs(obs) {}
  void self(const Vec& x,const Vec& x2,ID& nodePair,T eps) const override {
    //getID
    int l,r;
    CollisionDetector<3,3>::separate(nodePair,l,r);
    const Eigen::Matrix<int,3,1>& ml=_meshIss[l];
    const Eigen::Matrix<int,3,1>& mr=_meshIss[r];
    Eigen::Matrix<T,3,3> a;
    Eigen::Matrix<T,3,3> b;
    GJK::Vec3T pAL,pBL;
    bool intersect;
    GJK::T dist;
    //CCD
    if(CCDSelf(x,x2,ml,mr)) {
      //need to insert anyway
    } else {
      //only check x2
      for(int i=0; i<3; i++) {
        a.col(i)=x2.template segment<3>(ml[i]*3);
        b.col(i)=x2.template segment<3>(mr[i]*3);
      }
      //GJK
      SmallShape<3,3> A(a);
      SmallShape<3,3> B(b);
      dist=GJK::runGJK(A,B,
                       GJK::Mat3X4T::Identity(),
                       GJK::Mat3X4T::Identity(),
                       pAL,pBL,&intersect);
      if((T)dist>eps)
        nodePair=CollisionDetector<3,3>::INVALID_ID;
    }
  }
  void obs(const Vec& x,const Vec& x2,ID& nodePair,T eps) const override {
    //getID
    int l,r;
    CollisionDetector<3,3>::separate(nodePair,l,r);
    const Eigen::Matrix<int,3,1>& ml=_meshIss[l];
    const Eigen::Matrix<int,3,1>& mr=_obs.iss()[r];
    Eigen::Matrix<T,3,3> a;
    Eigen::Matrix<T,3,3> b;
    GJK::Vec3T pAL,pBL;
    bool intersect;
    GJK::T dist;
    //CCD
    if(CCDObs(x,x2,ml,mr)) {
      //need to insert anyway
    } else {
      //only check x2
      for(int i=0; i<3; i++) {
        a.col(i)=x2.template segment<3>(ml[i]*3);
        b.col(i)=_obs.vss()[mr[i]].template segment<3>(0).template cast<T>();
      }
      //GJK
      SmallShape<3,3> A(a);
      SmallShape<3,3> B(b);
      dist=GJK::runGJK(A,B,
                       GJK::Mat3X4T::Identity(),
                       GJK::Mat3X4T::Identity(),
                       pAL,pBL,&intersect);
      if((T)dist>eps)
        nodePair=CollisionDetector<3,3>::INVALID_ID;
    }
  }
 private:
  bool CCDSelf(const Vec& x,const Vec& x2,const Eigen::Matrix<int,3,1>& ml,const Eigen::Matrix<int,3,1>& mr) const {
    for(int d=0; d<3; d++) {
      if(CCDSelfVT(x,x2,ml[d],mr))
        return true;
      if(CCDSelfVT(x,x2,mr[d],ml))
        return true;
      for(int d2=0; d2<3; d2++)
        if(CCDSelfEE(x,x2,Eigen::Matrix<int,2,1>(ml[d],ml[(d+1)%3]),Eigen::Matrix<int,2,1>(mr[d2],mr[(d2+1)%3])))
          return true;
    }
    return false;
  }
  bool CCDSelfVT(const Vec& x,const Vec& x2,int ml,const Eigen::Matrix<int,3,1>& mr) const {
    return vertexFaceCCD(extract(x,ml),extract(x,mr[0]),extract(x,mr[1]),extract(x,mr[2]),
                         extract(x2,ml),extract(x2,mr[0]),extract(x2,mr[1]),extract(x2,mr[2]));
  }
  bool CCDSelfEE(const Vec& x,const Vec& x2,const Eigen::Matrix<int,2,1>& ml,const Eigen::Matrix<int,2,1>& mr) const {
    return edgeEdgeCCD(extract(x,ml[0]),extract(x,ml[1]),extract(x,mr[0]),extract(x,mr[1]),
                       extract(x2,ml[0]),extract(x2,ml[1]),extract(x2,mr[0]),extract(x2,mr[1]));
  }
  bool CCDObs(const Vec& x,const Vec& x2,const Eigen::Matrix<int,3,1>& ml,const Eigen::Matrix<int,3,1>& mr) const {
    for(int d=0; d<3; d++) {
      if(CCDObsVT(x,x2,ml[d],mr))
        return true;
      if(CCDObsTV(x,x2,ml,mr[d]))
        return true;
      for(int d2=0; d2<3; d2++)
        if(CCDObsEE(x,x2,Eigen::Matrix<int,2,1>(ml[d],ml[(d+1)%3]),Eigen::Matrix<int,2,1>(mr[d2],mr[(d2+1)%3])))
          return true;
    }
    return false;
  }
  bool CCDObsTV(const Vec& x,const Vec& x2,const Eigen::Matrix<int,3,1>& ml,int mr) const {
    return vertexFaceCCD(extract(x,ml[0]),extract(x,ml[1]),extract(x,ml[2]),extract(_obs,mr),
                         extract(x2,ml[0]),extract(x2,ml[1]),extract(x2,ml[2]),extract(_obs,mr));
  }
  bool CCDObsVT(const Vec& x,const Vec& x2,int ml,const Eigen::Matrix<int,3,1>& mr) const {
    return vertexFaceCCD(extract(x,ml),extract(_obs,mr[0]),extract(_obs,mr[1]),extract(_obs,mr[2]),
                         extract(x2,ml),extract(_obs,mr[0]),extract(_obs,mr[1]),extract(_obs,mr[2]));
  }
  bool CCDObsEE(const Vec& x,const Vec& x2,const Eigen::Matrix<int,2,1>& ml,const Eigen::Matrix<int,2,1>& mr) const {
    return edgeEdgeCCD(extract(x,ml[0]),extract(x,ml[1]),extract(_obs,mr[0]),extract(_obs,mr[1]),
                       extract(x2,ml[0]),extract(x2,ml[1]),extract(_obs,mr[0]),extract(_obs,mr[1]));
  }
  //extract
  static Eigen::Vector3d extract(const Vec& x,int off) {
    off*=3;
    return Eigen::Vector3d(x[off],x[off+1],x[off+2]);
  }
  static Eigen::Vector3d extract(const MeshExact& obs,int off) {
    return obs.vss()[off].template cast<double>();
  }
  const Iss& _meshIss;
  const MeshExact& _obs;
};
#else
template <>
class Narrow<3,3> {
 public:
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  typedef Eigen::Matrix<T,3,1> VecNT;
  typedef std::vector<Eigen::Matrix<int,3,1>> Iss;
  typedef unsigned long long ID;
  Narrow(const Iss& iss,const MeshExact& obs):_meshIss(iss),_obs(obs) {}
  int self(const Vec& x,const Vec& x2,ParallelVector<ID>& IDs,const std::unordered_map<ID,int>& existing,T eps) {
    int ret=0;
    auto& v=IDs.getVector();
    //build hash
    _vt.clear();
    _vtHash.clear();
    _ee.clear();
    _eeHash.clear();
    for(int i=0; i<(int)v.size(); i++)
      if(existing.find(v[i])!=existing.end())
        v[i]=CollisionDetector<3,3>::INVALID_ID;
      else if(v[i]!=CollisionDetector<3,3>::INVALID_ID) {
        int l,r;
        CollisionDetector<3,3>::separate(v[i],l,r);
        const Eigen::Matrix<int,3,1>& ml=_meshIss[l];
        const Eigen::Matrix<int,3,1>& mr=_meshIss[r];
        for(int d=0; d<3; d++) {
          //vt
          Eigen::Matrix<int,2,1> vt(ml[d],r);
          if(_vtHash.find(vt)==_vtHash.end()) {
            _vtHash[vt]=false;
            _vt.push_back(vt);
          }
          //tv
          Eigen::Matrix<int,2,1> tv(mr[d],l);
          if(_vtHash.find(tv)==_vtHash.end()) {
            _vtHash[tv]=false;
            _vt.push_back(tv);
          }
          //ee
          for(int d2=0; d2<3; d2++) {
            Eigen::Matrix<int,4,1> ee(ml[d],ml[(d+1)%3],mr[d2],mr[(d2+1)%3]);
            makeUnique(ee,true);
            if(_eeHash.find(ee)==_eeHash.end()) {
              _eeHash[ee]=false;
              _ee.push_back(ee);
            }
          }
        }
      }
    //parallel CCD
    OMP_PARALLEL_FOR_
    for(int i=0; i<(int)_vt.size(); i++) {
      const Eigen::Matrix<int,2,1>& vt=_vt[i];
      _vtHash[vt]=CCDSelfVT(x,x2,vt[0],_meshIss[vt[1]]);
    }
    OMP_PARALLEL_FOR_
    for(int i=0; i<(int)_ee.size(); i++) {
      const Eigen::Matrix<int,4,1>& ee=_ee[i];
      _eeHash[ee]=CCDSelfEE(x,x2,ee.template segment<2>(0),ee.template segment<2>(2));
    }
    //narrow
    OMP_PARALLEL_FOR_
    for(int i=0; i<(int)v.size(); i++)
      if(v[i]!=CollisionDetector<3,3>::INVALID_ID)
        self(x,x2,v[i],eps);
    return ret;
  }
  int obs(const Vec& x,const Vec& x2,ParallelVector<ID>& IDs,const std::unordered_map<ID,int>& existing,T eps) {
    int ret=0;
    auto& v=IDs.getVector();
    //build hash
    _vt.clear();
    _vtHash.clear();
    _tv.clear();
    _tvHash.clear();
    _ee.clear();
    _eeHash.clear();
    for(int i=0; i<(int)v.size(); i++)
      if(existing.find(v[i])!=existing.end())
        v[i]=CollisionDetector<3,3>::INVALID_ID;
      else if(v[i]!=CollisionDetector<3,3>::INVALID_ID) {
        int l,r;
        CollisionDetector<3,3>::separate(v[i],l,r);
        const Eigen::Matrix<int,3,1>& ml=_meshIss[l];
        const Eigen::Matrix<int,3,1>& mr=_obs.iss()[r];
        for(int d=0; d<3; d++) {
          //vt
          Eigen::Matrix<int,2,1> vt(ml[d],r);
          if(_vtHash.find(vt)==_vtHash.end()) {
            _vtHash[vt]=false;
            _vt.push_back(vt);
          }
          //tv
          Eigen::Matrix<int,2,1> tv(l,mr[d]);
          if(_tvHash.find(tv)==_tvHash.end()) {
            _tvHash[tv]=false;
            _tv.push_back(tv);
          }
          //ee
          for(int d2=0; d2<3; d2++) {
            Eigen::Matrix<int,4,1> ee(ml[d],ml[(d+1)%3],mr[d2],mr[(d2+1)%3]);
            makeUnique(ee,false);
            if(_eeHash.find(ee)==_eeHash.end()) {
              _eeHash[ee]=false;
              _ee.push_back(ee);
            }
          }
        }
      }
    //parallel CCD
    OMP_PARALLEL_FOR_
    for(int i=0; i<(int)_vt.size(); i++) {
      const Eigen::Matrix<int,2,1>& vt=_vt[i];
      _vtHash[vt]=CCDObsVT(x,x2,vt[0],_meshIss[vt[1]]);
    }
    OMP_PARALLEL_FOR_
    for(int i=0; i<(int)_tv.size(); i++) {
      const Eigen::Matrix<int,2,1>& tv=_tv[i];
      _tvHash[tv]=CCDObsTV(x,x2,_meshIss[tv[0]],tv[1]);
    }
    OMP_PARALLEL_FOR_
    for(int i=0; i<(int)_ee.size(); i++) {
      const Eigen::Matrix<int,4,1>& ee=_ee[i];
      _eeHash[ee]=CCDObsEE(x,x2,ee.template segment<2>(0),ee.template segment<2>(2));
    }
    //narrow
    OMP_PARALLEL_FOR_
    for(int i=0; i<(int)v.size(); i++)
      if(v[i]!=CollisionDetector<3,3>::INVALID_ID)
        obs(x,x2,v[i],eps);
    return ret;
  }
  void self(const Vec& x,const Vec& x2,ID& nodePair,T eps) const {
    //getID
    int l,r;
    CollisionDetector<3,3>::separate(nodePair,l,r);
    const Eigen::Matrix<int,3,1>& ml=_meshIss[l];
    const Eigen::Matrix<int,3,1>& mr=_meshIss[r];
    Eigen::Matrix<T,3,3> a;
    Eigen::Matrix<T,3,3> b;
    GJK::Vec3T pAL,pBL;
    bool intersect;
    GJK::T dist;
    //CCD
    if(CCDSelf(x,x2,l,r)) {
      //need to insert anyway
    } else {
      //only check x2
      for(int i=0; i<3; i++) {
        a.col(i)=x2.template segment<3>(ml[i]*3);
        b.col(i)=x2.template segment<3>(mr[i]*3);
      }
      //GJK
      SmallShape<3,3> A(a);
      SmallShape<3,3> B(b);
      dist=GJK::runGJK(A,B,
                       GJK::Mat3X4T::Identity(),
                       GJK::Mat3X4T::Identity(),
                       pAL,pBL,&intersect);
      if((T)dist>eps)
        nodePair=CollisionDetector<3,3>::INVALID_ID;
    }
  }
  void obs(const Vec& x,const Vec& x2,ID& nodePair,T eps) const {
    //getID
    int l,r;
    CollisionDetector<3,3>::separate(nodePair,l,r);
    const Eigen::Matrix<int,3,1>& ml=_meshIss[l];
    const Eigen::Matrix<int,3,1>& mr=_obs.iss()[r];
    Eigen::Matrix<T,3,3> a;
    Eigen::Matrix<T,3,3> b;
    GJK::Vec3T pAL,pBL;
    bool intersect;
    GJK::T dist;
    //CCD
    if(CCDObs(x,x2,l,r)) {
      //need to insert anyway
    } else {
      //only check x2
      for(int i=0; i<3; i++) {
        a.col(i)=x2.template segment<3>(ml[i]*3);
        b.col(i)=_obs.vss()[mr[i]].template segment<3>(0).template cast<T>();
      }
      //GJK
      SmallShape<3,3> A(a);
      SmallShape<3,3> B(b);
      dist=GJK::runGJK(A,B,
                       GJK::Mat3X4T::Identity(),
                       GJK::Mat3X4T::Identity(),
                       pAL,pBL,&intersect);
      if((T)dist>eps)
        nodePair=CollisionDetector<3,3>::INVALID_ID;
    }
  }
 protected:
  bool CCDSelf(const Vec& x,const Vec& x2,int l,int r) const {
    const Eigen::Matrix<int,3,1>& ml=_meshIss[l];
    const Eigen::Matrix<int,3,1>& mr=_meshIss[r];
    for(int d=0; d<3; d++) {
      Eigen::Matrix<int,2,1> vt(ml[d],r);
      if(_vtHash.find(vt)->second)
        return true;
      Eigen::Matrix<int,2,1> tv(mr[d],l);
      if(_vtHash.find(tv)->second)
        return true;
      for(int d2=0; d2<3; d2++) {
        Eigen::Matrix<int,4,1> ee(ml[d],ml[(d+1)%3],mr[d2],mr[(d2+1)%3]);
        makeUnique(ee,true);
        if(_eeHash.find(ee)->second)
          return true;
      }
    }
    return false;
  }
  bool CCDSelfVT(const Vec& x,const Vec& x2,int ml,const Eigen::Matrix<int,3,1>& mr) const {
    return vertexFaceCCD(extract(x,ml),extract(x,mr[0]),extract(x,mr[1]),extract(x,mr[2]),
                         extract(x2,ml),extract(x2,mr[0]),extract(x2,mr[1]),extract(x2,mr[2]));
  }
  bool CCDSelfEE(const Vec& x,const Vec& x2,const Eigen::Matrix<int,2,1>& ml,const Eigen::Matrix<int,2,1>& mr) const {
    return edgeEdgeCCD(extract(x,ml[0]),extract(x,ml[1]),extract(x,mr[0]),extract(x,mr[1]),
                       extract(x2,ml[0]),extract(x2,ml[1]),extract(x2,mr[0]),extract(x2,mr[1]));
  }
  bool CCDObs(const Vec& x,const Vec& x2,int l,int r) const {
    const Eigen::Matrix<int,3,1>& ml=_meshIss[l];
    const Eigen::Matrix<int,3,1>& mr=_obs.iss()[r];
    for(int d=0; d<3; d++) {
      Eigen::Matrix<int,2,1> vt(ml[d],r);
      if(_vtHash.find(vt)->second)
        return true;
      Eigen::Matrix<int,2,1> tv(l,mr[d]);
      if(_tvHash.find(tv)->second)
        return true;
      for(int d2=0; d2<3; d2++) {
        Eigen::Matrix<int,4,1> ee(ml[d],ml[(d+1)%3],mr[d2],mr[(d2+1)%3]);
        makeUnique(ee,false);
        if(_eeHash.find(ee)->second)
          return true;
      }
    }
    return false;
  }
  bool CCDObsTV(const Vec& x,const Vec& x2,const Eigen::Matrix<int,3,1>& ml,int mr) const {
    return vertexFaceCCD(extract(x,ml[0]),extract(x,ml[1]),extract(x,ml[2]),extract(_obs,mr),
                         extract(x2,ml[0]),extract(x2,ml[1]),extract(x2,ml[2]),extract(_obs,mr));
  }
  bool CCDObsVT(const Vec& x,const Vec& x2,int ml,const Eigen::Matrix<int,3,1>& mr) const {
    return vertexFaceCCD(extract(x,ml),extract(_obs,mr[0]),extract(_obs,mr[1]),extract(_obs,mr[2]),
                         extract(x2,ml),extract(_obs,mr[0]),extract(_obs,mr[1]),extract(_obs,mr[2]));
  }
  bool CCDObsEE(const Vec& x,const Vec& x2,const Eigen::Matrix<int,2,1>& ml,const Eigen::Matrix<int,2,1>& mr) const {
    return edgeEdgeCCD(extract(x,ml[0]),extract(x,ml[1]),extract(_obs,mr[0]),extract(_obs,mr[1]),
                       extract(x2,ml[0]),extract(x2,ml[1]),extract(_obs,mr[0]),extract(_obs,mr[1]));
  }
  //extract
  static Eigen::Vector3d extract(const Vec& x,int off) {
    off*=3;
    return Eigen::Vector3d((double)x[off],(double)x[off+1],(double)x[off+2]);
  }
  static Eigen::Vector3d extract(const MeshExact& obs,int off) {
    return obs.vss()[off].template cast<double>();
  }
  static void makeUnique(Eigen::Matrix<int,4,1>& ee,bool reorder) {
    sort2(ee[0],ee[1]);
    sort2(ee[2],ee[3]);
    if(reorder)
      if(ee[0]<ee[2] || (ee[0]==ee[2] && ee[1]<ee[3])) {
        std::swap(ee[0],ee[2]);
        std::swap(ee[1],ee[3]);
      }
  }
  //data
  const Iss& _meshIss;
  const MeshExact& _obs;
  std::vector<Eigen::Matrix<int,2,1>> _vt;
  std::unordered_map<Eigen::Matrix<int,2,1>,bool,EdgeHash> _vtHash;
  std::vector<Eigen::Matrix<int,2,1>> _tv;
  std::unordered_map<Eigen::Matrix<int,2,1>,bool,EdgeHash> _tvHash;
  std::vector<Eigen::Matrix<int,4,1>> _ee;
  std::unordered_map<Eigen::Matrix<int,4,1>,bool,EdgeEdgeHash> _eeHash;
};
#endif
#endif
}

#endif
