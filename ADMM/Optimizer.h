#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "OptimizerTerm.h"

namespace PHYSICSMOTION {
struct OptimizerParam {
  typedef FLOAT T;
  enum TYPE {
    ADMM,
    GD,
    NEWTON,
    UNKNOWN,
  };
  T _tolXY=0;
  T _tolG=1e-4f;
  T _tolGYZ=1e-6f;
  T _initAlpha=1;
  T _tolAlpha=1e-10f;
  T _decAlpha=0.5f;
  T _incAlpha=1/0.8f;
  T _initBeta=1;
  T _initBetaX=1;
  T _initBetaY=0;
  T _collisionRemoveMargin=3;
  int _maxIter=1e4;
  int _checkI=1;
  int _printI=50;
  int _debugGradientI=0;
  int _collisionRemoveI=200;
  bool _ensureMonotonic=true;   //turn this off will make BC-ADMM fail to converge, for debug/comparison only
  TYPE _type=ADMM;
};
template <int N>
class Optimizer;
template <int N>
class CollisionChecker {
 public:
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  virtual ~CollisionChecker() = default;
  virtual Eigen::Matrix<int,2,1> generate(const Vec& x,const Vec& x2,Optimizer<N>& opt,bool noNeighbhor,bool parallel,int restoreSlot=-1)=0;
  virtual Eigen::Matrix<int,2,1> remove(const Vec& x,Optimizer<N>& opt,T margin)=0;
  virtual std::string info(const Optimizer<N>& opt) const=0;
};
template <int N>
class Optimizer {
 public:
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  typedef Eigen::Matrix<T,N,1> VecNT;
  typedef Eigen::Matrix<T,N,N> MatNT;
  typedef Eigen::Triplet<T,int> STrip;
  typedef ParallelVector<STrip> STrips;
  typedef Eigen::SparseMatrix<T,0,int> SMatT;
  struct LinearConstraint {
    std::unordered_map<int,std::pair<T,VecNT>> _fixes;
    SMatT _general;
  };
  Optimizer();
  virtual ~Optimizer() = default;
  Vec& x() {
    return _x;
  }
  const Vec& x() const {
    return _x;
  }
  void setCB(std::function<bool()> cb);
  void addARAPElement(std::array<int,N+1> id,T k,const MatNT& F0);
  void addSpring(std::array<int,2> id,T k,T l,T lL,T lH);
  void addSmoothL1(std::array<int,2> id,T k,T eps=1e-2f);
  void addQuadratic(const SMatT* H,const Vec* g);
  virtual void assemble(int nX,const LinearConstraint& cons= {});
  virtual void optimize(const OptimizerParam& param)=0;
  template <typename TERM>
  std::shared_ptr<TERM> findTerm() const {
    for(const auto& term:_gss)
      if(std::dynamic_pointer_cast<TERM>(term))
        return std::dynamic_pointer_cast<TERM>(term);
    return NULL;
  }
  template <typename TERM>
  void insertTerm(std::shared_ptr<TERM> term) {
    _gss.push_back(term);
  }
  void insertCollisionDetector(std::shared_ptr<CollisionChecker<N>> coll);
  void save(int id,int mask);
  void load(int id,int mask);
 protected:
  void init(T tolG);
  void project(Vec& G);
  bool updateZ(T tolG);
  //data
  Vec _x,_d;
  std::vector<Vec> _xSaved;
  //matrix A
  SMatT _Cons;
  Eigen::SparseLU<SMatT> _CCTInv;
  //term 1: quadratic term: f
  SMatT _H;
  Vec _g;
  //term 2: g
  std::vector<std::shared_ptr<OptimizerTerm>> _gss;
  //collision checker
  std::shared_ptr<CollisionChecker<N>> _coll;
  //callback
  std::function<bool()> _cb=[]() {
    return true;
  };
};
}

#endif
