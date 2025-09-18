#ifndef DEFORMABLE_H
#define DEFORMABLE_H

#include "Optimizer.h"
#include "ADMMUAVTrajectory.h"
#include <Environment/MeshExact.h>
#include <unordered_map>
#include <unordered_set>

namespace PHYSICSMOTION {
//(M/(dt*dt*2))*||x-2*xL+xLL||^2+P(x)-g^T*M*x+P(x)
template <int N>
class Deformable {
 public:
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  struct ARAPElement {
    Eigen::Matrix<int,N+1,1> _indices;
    T _stiffness;
  };
  struct Spring {
    T _k=1e3f,_b=1e2f,_cL=0,_cH=0;
  };
  typedef Eigen::Matrix<T,N,1> VecNT;
  typedef Eigen::Matrix<T,N,N> MatNT;
  typedef Eigen::Matrix<T,N,-1> MatNXT;
  typedef Eigen::Triplet<T,int> STrip;
  typedef ParallelVector<STrip> STrips;
  typedef Eigen::SparseMatrix<T,0,int> SMatT;
  typedef std::unordered_map<int,std::pair<std::unordered_set<int>,Spring>> Ring;
  typedef std::unordered_map<Eigen::Matrix<int,2,1>,Spring,EdgeHash> Springs;
  typedef std::unordered_map<Eigen::Matrix<int,2,1>,int,EdgeHash> EdgeMap;
  typedef std::unordered_map<Eigen::Matrix<int,2,1>,T,EdgeHash> L1Constraint;
  typedef typename Optimizer<N>::LinearConstraint LinearConstraint;
 public:
  const auto& getObs() const {
    return _obs;
  }
  const auto& getEss() const {
    return _ess;
  }
  const auto& getTss() const {
    return _tss;
  }
  const auto& getSss() const {
    return _sss;
  }
  const auto& getCons() const {
    return _cons;
  }
  const auto& getL1ss() const {
    return _l1ss;
  }
  const auto& x() const {
    return _x;
  }
  const auto& getTrajSubd() const {
    return _trajSubd;
  }
  const auto& getUAVTraj() const {
    return _UAVTraj;
  }
  void setAgents(const Vec& x,const Vec& xT,T r);
  void setARAP3D(const MeshExact& mesh,double feature_angle=91,double size=1);      //we generate tet mesh and replace the given mesh with surface of tet mesh
  void setARAP2D(const MeshExact& mesh);
  void setMassSpring(const MeshExact& mesh);
  bool setUAVTrajectory(const std::vector<VecNT>& waypoints,int numSegmentBeg,int numSegmentEnd);
  void setObstacle(const MeshExact& mesh);
  void fix(int vid,T k,const VecNT* xStar=NULL);
  void fix(std::function<bool(const VecNT& a)> f,T k);
  auto& getFix() {return _cons._fixes;}
  int findClosestVid(const VecNT& pos) const;
  void addL1(int id,const VecNT& pos,T k);
  void addL1(std::array<int,2> id,T k);
  T dt() const;
  T r() const;
  T actualR() const;
  void setDt(T dt);
  void setR(T r);
  void setK(T k);
  void setB(T b);
  void setCL(T cL);
  void setCH(T cH);
  void setMargin(T margin);
  void setCollCoef(T collCoef);
  void setG(const VecNT& g);
  void setUAVTrajResolution(int RES);
  void setCB(std::function<bool()> cb);
  void solve(const OptimizerParam& param=OptimizerParam());
  const MeshExact& buildSelfCollisionMesh();
 private:
  void assemble();
  void fix(STrips* HTrips,Vec* g) const;
  void buildRing1(Ring& ring1) const;
  void buildRing2(const Ring& ring1,Ring& ring2) const;
  void buildMassMatrix(T coef);
  void buildPointMassMatrix(STrips& HTrips,T coef) const;
  void buildFEMMassMatrix(STrips& HTrips,T coef) const;
  //param
  SMatT _mass,_trajSubd;
  Vec _x,_xL,_xLL,_x0,_xT;
  VecNT _g=VecNT::Zero();
  Spring _springMaterial;
  T _dt=1e-2f,_r=0,_m=1,_margin=.2f,_collCoef=.1f;
  //topology
  Springs _ess;                                 //mass springs
  std::vector<ARAPElement> _tss;                //ARAP elements
  std::vector<Eigen::Matrix<int,N,1>> _sss;     //surface of ARAP elements
  LinearConstraint _cons;                       //linear constraints
  L1Constraint _l1ss;                           //l1 terms
  ADMMUAVTrajectory<N> _UAVTraj;                //UAV trajectory
  MeshExact _self,_obs;                         //obstacle
  //solver
  std::shared_ptr<Optimizer<N>> _solver;
  //callback
  std::function<bool()> _cb=[]() {
    return true;
  };
  bool _hasCB=false;
};
}

#endif
