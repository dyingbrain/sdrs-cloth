#include "OptimizerTerm.h"

namespace PHYSICSMOTION {
OptimizerTerm::VecM OptimizerTerm::y() {
  if(_y.cols()!=n())
    _y.conservativeResize(m(),n());
  return VecM(_y.data(),_y.size());
}
OptimizerTerm::VecCM OptimizerTerm::y() const {
  return VecCM(_y.data(),_y.size());
}
OptimizerTerm::VecCM OptimizerTerm::yLast() const {
  return VecCM(_yLast.data(),_yLast.size());
}
OptimizerTerm::VecM OptimizerTerm::Ax() {
  if(_Ax.cols()!=n())
    _Ax.conservativeResize(m(),n());
  return VecM(_Ax.data(),_Ax.size());
}
OptimizerTerm::VecCM OptimizerTerm::Ax() const {
  return VecCM(_Ax.data(),_Ax.size());
}
OptimizerTerm::VecM OptimizerTerm::L() {
  return VecM(_L.data(),_L.size());
}
OptimizerTerm::VecCM OptimizerTerm::L() const {
  return VecCM(_L.data(),_L.size());
}
OptimizerTerm::VecCM OptimizerTerm::G() const {
  return VecCM(_G.data(),_G.size());
}
const OptimizerTerm::SMatT& OptimizerTerm::A() const {
  return _A;
}
void OptimizerTerm::assembleA(int* nX) {
  int nRows=0,nCols=nX?*nX:_A.cols();
  for(auto beg=_E.begin(),end=_E.end(); beg!=end; beg++)
    nRows=std::max(nRows,beg->row());
  _A.resize(_E.getVector().empty()?0:nRows+1,nCols);
  _A.setFromTriplets(_E.begin(),_E.end());
}
int OptimizerTerm::m() const {
  if(n()==0)
    return 0;
  return _A.rows()/n();
}
OptimizerTerm::VecM OptimizerTerm::y0() {
  static Vec y0V=Vec::Zero(0);
  return VecM(y0V.data(),0);
}
OptimizerTerm::VecCM OptimizerTerm::y0() const {
  static Vec y0V=Vec::Zero(0);
  return VecCM(y0V.data(),0);
}
OptimizerTerm::VecCM OptimizerTerm::G0() const {
  static Vec G0V=Vec::Zero(0);
  return VecCM(G0V.data(),0);
}
OptimizerTerm::T OptimizerTerm::evalGDirect(bool calcG,SMatT* H,int y0Off,bool projPSD) {
  //Direct evaluation is not implemented by default, it is more involved than evalG
  //FUNCTION_NOT_IMPLEMENTED
  return 0;
}
bool OptimizerTerm::updateZ(T tolG) {
  return true;
}
void OptimizerTerm::updateL(T beta) {
  _beta=beta;
  OMP_PARALLEL_FOR_
  for(int i=0; i<(int)_y.cols(); i++)
    _L.col(i)+=beta*(_Ax.col(i)-_y.col(i));
}
void OptimizerTerm::reset(int mask) {
  if(mask&MASK_Y)
    _y.resize(0,0);
  if(mask&MASK_L)
    _L.resize(0,0);
}
void OptimizerTerm::save(int id,int mask) {
  if(id>=(int)_saved.size())
    _saved.resize(id+1);
  if(!_saved[id])
    _saved[id]=copy();
  if(mask&MASK_Y)
    _saved[id]->_y=_y;
  if(mask&MASK_L)
    _saved[id]->_L=_L;
}
void OptimizerTerm::load(int id,int mask) {
  ASSERT(id<(int)_saved.size() && _saved[id])
  if(mask&MASK_Y)
    _y=_saved[id]->_y;
  if(mask&MASK_L)
    _L=_saved[id]->_L;
}
void OptimizerTerm::setDirectMode(bool direct) {
  _directMode=direct;
}
//helper
void OptimizerTerm::initializeL() {
  int cols=n();
  if(cols==0)
    return;
  int rows=m();
  int colsL=_L.cols();
  _L.conservativeResize(rows,cols);
  OMP_PARALLEL_FOR_
  for(int i=colsL; i<cols; i++)
    _L.col(i)=_G.col(i);
}
void OptimizerTerm::assembleHessian(SMatT& H,int y0Off) {
  //build HY
  SMatT HY;
  int offY=0,offY0=0,HYSize=0;
  {
    STrips HYTrips;
    for(int i=0; i<n(); i++) {
      const HessianBlock& HBlk=_HBlks[i];
      HYSize+=HBlk._blk.rows();
      offY0+=HBlk._nY;
    }
    ASSERT(offY0==_A.rows())
    HY.resize(HYSize,HYSize);
    for(int i=0; i<n(); i++) {
      const HessianBlock& HBlk=_HBlks[i];
      const auto mappedId=[&](int i) {
        return i<HBlk._nY?offY+i:offY0+i-HBlk._nY;
      };
      for(int r=0; r<HBlk._blk.rows(); r++)
        for(int c=0; c<HBlk._blk.cols(); c++)
          HYTrips.push_back(STrip(mappedId(r),mappedId(c),HBlk._blk(r,c)));
      offY+=HBlk._nY;
      offY0+=HBlk._blk.rows()-HBlk._nY;
    }
    HY.setFromTriplets(HYTrips.begin(),HYTrips.end());
    ASSERT(offY==_A.rows())
    ASSERT(offY0==HYSize)
  }
  //build AI
  SMatT AI;
  {
    STrips AITrips=_E;
    AI.resize(HYSize,H.cols());
    addBlockId<T>(AITrips,offY,y0Off,HYSize-offY,1);
    AI.setFromTriplets(AITrips.begin(),AITrips.end());
  }
  H+=AI.transpose()*HY*AI;
}
void OptimizerTerm::removeColumns(const std::unordered_set<int>& cols,int rows) {
  //compact
  int newSize=0;
  for(int col=0; col<n(); col++)
    if(cols.find(col)==cols.end()) {
      _y.col(newSize)=_y.col(col);
      if(_yLast.cols()==n())
        _yLast.col(newSize)=_yLast.col(col);
      if(_L.cols()==n())
        _L.col(newSize)=_L.col(col);
      if(_Ax.cols()==n())
        _Ax.col(newSize)=_Ax.col(col);
      if(_G.cols()==n())
        _G.col(newSize)=_G.col(col);
      newSize++;
    }
  //resize
  _y.conservativeResize(_y.rows(),newSize);
  if(_yLast.cols()==n())
    _yLast.conservativeResize(_yLast.rows(),newSize);
  if(_L.cols()==n())
    _L.conservativeResize(_L.rows(),newSize);
  if(_Ax.cols()==n())
    _Ax.conservativeResize(_Ax.rows(),newSize);
  if(_G.cols()==n())
    _G.conservativeResize(_G.rows(),newSize);
  //update _A
  int row=0;
  STrips trips;
  SMatT ASelect;
  for(int i=0; i<_A.rows(); i++) {
    int col=i/rows;
    if(cols.find(col)==cols.end())
      trips.push_back(STrip(row++,i,1));
  }
  ASelect.resize(row,_A.rows());
  ASelect.setFromTriplets(trips.begin(),trips.end());
  _A=ASelect*_A;
  //update _E
  _E.clear();
  addBlock<T,0,int>(_E,0,0,_A);
}
}
