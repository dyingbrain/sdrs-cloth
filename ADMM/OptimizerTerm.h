#ifndef OPTIMIZER_TERM_H
#define OPTIMIZER_TERM_H

#include <Eigen/Dense>
#include <ConvexHull/Barrier.h>
#include <Utils/SparseUtils.h>
#include <Utils/Pragma.h>
#include <unordered_set>

//#define OPTIMIZE_ON_SPHERE
namespace PHYSICSMOTION {
class OptimizerTerm {
 public:
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  typedef Eigen::Triplet<T,int> STrip;
  typedef ParallelVector<STrip> STrips;
  typedef Eigen::SparseMatrix<T,0,int> SMatT;
  struct HessianBlock {
    MatT _blk;
    int _nY;
  };
  enum SaveMask {
    MASK_X      =1,
    MASK_Y      =2,
    MASK_Y0     =4,
    MASK_L      =8,
    MASK_Z      =16,
    MASK_ALPHA  =32,
  };
  virtual ~OptimizerTerm() = default;
  VecM y();
  VecCM y() const;
  VecCM yLast() const;
  VecM Ax();
  VecCM Ax() const;
  VecM L();
  VecCM L() const;
  VecCM G() const;
  const SMatT& A() const;
  void assembleA(int* nX);
  int m() const;
  //virtual function
  virtual VecM y0();
  virtual VecCM y0() const;
  virtual VecCM G0() const;
  virtual int n() const=0;
  virtual std::shared_ptr<OptimizerTerm> copy() const=0;
  virtual T evalG(bool calcG,bool initL,SMatT* H=NULL,int y0Off=-1)=0;
  virtual bool updateY(T betaY,T beta,T tolG)=0;
  virtual bool updateZ(T tolG);
  virtual void updateL(T beta);
  virtual void reset(int mask=MASK_Y|MASK_Y0|MASK_L|MASK_Z|MASK_ALPHA);
  virtual void save(int id,int mask);
  virtual void load(int id,int mask);
 protected:
  void initializeL();
  void assembleHessian(SMatT& H,int y0Off);
  void removeColumns(const std::unordered_set<int>& cols,int rows);
  //data
  std::vector<std::shared_ptr<OptimizerTerm>> _saved;
  std::vector<HessianBlock> _HBlks;
  MatT _y,_yLast,_L,_Ax,_G;
  STrips _E;
  SMatT _A;
  T _beta,_betaY;
  bool _evalgOnly;
};
}

#endif
