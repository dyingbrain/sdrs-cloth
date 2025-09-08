#ifndef SMALL_SHAPE_H
#define SMALL_SHAPE_H

#include <Environment/ShapeExact.h>

namespace PHYSICSMOTION {
template <int N,int M>
class SmallShape : public ShapeExact {
 public:
  typedef Eigen::Matrix<FLOAT,N,M> DataType;
  SmallShape(const DataType& pts):_pts(pts) {}
  std::string type() const {
    return typeid(SmallShape<N,M>).name();
  }
  const BBoxExact& getBB() const override {
    FUNCTION_NOT_IMPLEMENTED
    static BBoxExact bb;
    return bb;
  }
  Vec3T support(const Vec3T& D,int& id) const override {
    id=0;
    Vec3T pt=Vec3T::Zero();
    pt.template segment<N>(0)=_pts.col(0).template cast<T>();
    T dist,maxDist=pt.dot(D);
    for(int i=1; i<M; i++) {
      pt.template segment<N>(0)=_pts.col(i).template cast<T>();
      dist=pt.dot(D);
      if(dist>maxDist) {
        maxDist=dist;
        id=i;
      }
    }
    pt.template segment<N>(0)=_pts.col(id).template cast<T>();
    return pt;
  }
 private:
  const DataType& _pts;
};
template <int N>
class SmallShapeDynamic : public ShapeExact {
 public:
  typedef Eigen::Matrix<FLOAT,N,-1> DataType;
  SmallShapeDynamic(const DataType& pts):_pts(pts) {}
  std::string type() const {
    return typeid(SmallShapeDynamic<N>).name();
  }
  const BBoxExact& getBB() const override {
    FUNCTION_NOT_IMPLEMENTED
    static BBoxExact bb;
    return bb;
  }
  Vec3T support(const Vec3T& D,int& id) const override {
    id=0;
    Vec3T pt=Vec3T::Zero();
    pt.template segment<N>(0)=_pts.col(0).template cast<T>();
    T dist,maxDist=pt.dot(D);
    for(int i=1; i<_pts.cols(); i++) {
      pt.template segment<N>(0)=_pts.col(i).template cast<T>();
      dist=pt.dot(D);
      if(dist>maxDist) {
        maxDist=dist;
        id=i;
      }
    }
    pt.template segment<N>(0)=_pts.col(id).template cast<T>();
    return pt;
  }
 private:
  const DataType _pts;
};
}

#endif
