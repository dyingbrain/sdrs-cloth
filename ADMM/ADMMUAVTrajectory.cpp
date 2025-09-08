#include "ADMMUAVTrajectory.h"
#include <Environment/GJK.h>
#include <Utils/VTKWriter.h>
#include "SmallShape.h"
#include <stack>

namespace PHYSICSMOTION {
//construction
template <int N>
void ADMMUAVTrajectory<N>::setNodes(const std::vector<VecNT>& waypoints,int order,int numSegment,bool debug) {
  setNodes(order,numSegment);
  int numNodes=_traj.getNumNode();
  //evenly divide the trajectory into nodes
  _nodes.resize(numNodes,N);
  //calculate total length
  T totalLength=0;
  for(int i=1; i<(int)waypoints.size(); i++)
    totalLength+=(waypoints[i]-waypoints[i-1]).norm();
  //sample waypoints
  int currId=0;
  T currLen=0,currFrac=0,step=totalLength/(numNodes-1);
  _nodes.row(0)=waypoints[0];
  for(int i=1; i<numNodes; i++) {
    //sample
    currLen=0;
    while(currLen<step) {
      if(currId==(int)waypoints.size()-1) {
        currId=(int)waypoints.size()-1;
        currFrac=0;
        break;
      }
      T segLen=(waypoints[currId+1]-waypoints[currId]).norm();
      if(currLen+segLen-currFrac>step) {
        currFrac+=step-currLen;
        currLen=step;
      } else {
        currLen+=segLen-currFrac;
        currFrac=0;
        currId++;
      }
    }
    //assign node
    if(currId==(int)waypoints.size()-1)
      _nodes.row(i)=waypoints.back();
    else {
      T segLen=(waypoints[currId+1]-waypoints[currId]).norm();
      _nodes.row(i)=waypoints[currId]+(waypoints[currId+1]-waypoints[currId])*currFrac/segLen;
    }
  }
  if(debug)
    writeNodesVTK("nodes.vtk");
  //convert x to CP
  _CPs=_nodeToCP*_nodes;
  //debug to see if constraints are satisfied
  if(debug)
    std::cout << "CP constraint error: " << (_CPCons*_CPs).norm() << std::endl;
}
template <int N>
bool ADMMUAVTrajectory<N>::setNodesExhaustive(const std::vector<VecNT>& waypoints,const MeshExact& obs,T margin,int order,int numSegmentBeg,int numSegmentEnd) {
  for(int numSegment=numSegmentBeg; numSegment<numSegmentEnd; numSegment++) {
    setNodes(waypoints,order,numSegment);
    if(!intersect(obs,margin))
      return true;
  }
  return false;
}
template <int N>
bool ADMMUAVTrajectory<N>::intersect(const MeshExact& obs,T margin) const {
  bool ret=false;
  int numSegment=_traj.getNumSegment();
  OMP_PARALLEL_FOR_
  for(int i=0; i<numSegment; i++)
    if(intersect(i,obs,margin))
      ret=true;
  return ret;
}
template <int N>
void ADMMUAVTrajectory<N>::setNodes(int order,int numSegment) {
  _traj=ThetaTrajectory<T>(1,order,numSegment);
  _nodeToCP=_traj.getControlPointCoeffAll(0);
  //build CPCons
  int rowOff=0;
  STrips trips;
  for(int i=1; i<numSegment; i++) {
    int ctrNode=i*order;
    //first constraint
    trips.push_back(STrip(rowOff,ctrNode,1));
    trips.push_back(STrip(rowOff,ctrNode-2, 0.25));
    trips.push_back(STrip(rowOff,ctrNode+2,-0.25));
    trips.push_back(STrip(rowOff,ctrNode-1,-1));
    rowOff++;
    //second constraint
    trips.push_back(STrip(rowOff,ctrNode,1));
    trips.push_back(STrip(rowOff,ctrNode-2,-0.25));
    trips.push_back(STrip(rowOff,ctrNode+2, 0.25));
    trips.push_back(STrip(rowOff,ctrNode+1,-1));
    rowOff++;
  }
  _CPCons.resize(rowOff,_nodeToCP.rows());
  _CPCons.setFromTriplets(trips.begin(),trips.end());
}
template <int N>
void ADMMUAVTrajectory<N>::buildSmoothness(int d) {
  int order=_traj.getOrder();
  _smoothH.resize(_nodeToCP.rows(),_nodeToCP.rows());
  MatT h=_traj.getCurve().getBasisCrossIntegral(d);
  for(int i=0; i<_traj.getNumSegment(); i++) {
    //select ith segment
    STrips trips;
    for(int j=0; j<=order; j++)
      trips.push_back(STrip(j,i*order+j,1));
    //assemble Hessian term
    SMatT seg;
    seg.resize(order+1,_nodeToCP.rows());
    seg.setFromTriplets(trips.begin(),trips.end());
    _smoothH+=seg.transpose()*h.sparseView()*seg;
  }
}
//getter
template <int N>
typename ADMMUAVTrajectory<N>::Vec ADMMUAVTrajectory<N>::getCP() const {
  MatNXT CPs=_CPs.transpose();
  return CPs.reshaped();
}
template <int N>
int ADMMUAVTrajectory<N>::numCP() const {
  return _CPs.rows();
}
template <int N>
typename ADMMUAVTrajectory<N>::SMatT ADMMUAVTrajectory<N>::getNodeToCP() const {
  return kronecker(_nodeToCP,N);
}
template <int N>
typename ADMMUAVTrajectory<N>::SMatT ADMMUAVTrajectory<N>::getCPCons() const {
  return kronecker(_CPCons,N);
}
template <int N>
typename ADMMUAVTrajectory<N>::SMatT ADMMUAVTrajectory<N>::getSmoothH() const {
  return kronecker(_smoothH,N);
}
template <int N>
typename ADMMUAVTrajectory<N>::SMatT ADMMUAVTrajectory<N>::subdivide(int RES) const {
  int rowOff=0;
  int order=_traj.getOrder();
  int numSegment=_traj.getNumSegment();
  //build triplets
  Vec coeff;
  STrips trips;
  //add first point
  coeff=_traj.getCurve().getCoefficientMatrix()*_traj.getCurve().getBasis(0);
  for(int c=0; c<order+1; c++)
    trips.push_back(STrip(rowOff,c,coeff[c]));
  rowOff++;
  //add subsequent points
  for(int i=0; i<numSegment; i++) {
    for(int t=1; t<=RES; t++) {
      coeff=_traj.getCurve().getCoefficientMatrix()*_traj.getCurve().getBasis(t/(T)RES);
      for(int c=0; c<order+1; c++)
        trips.push_back(STrip(rowOff,i*order+c,coeff[c]));
      rowOff++;
    }
  }
  //build matrix
  SMatT ret;
  ret.resize(rowOff,_nodeToCP.rows());
  ret.setFromTriplets(trips.begin(),trips.end());
  return kronecker(ret,N);
}
template <int N>
const ThetaTrajectory<typename ADMMUAVTrajectory<N>::T>& ADMMUAVTrajectory<N>::getTraj() const {
  return _traj;
}
template <int N>
void ADMMUAVTrajectory<N>::writeTrajVTK(const std::string& path,const Vec* pssIn,int RES) const {
  Vec pss;
  if(!pssIn) {
    MatNXT CPs=_CPs.transpose();
    pss=subdivide(RES)*CPs.reshaped();
  } else pss=*pssIn;
  std::vector<Eigen::Matrix<double,3,1>> vss;
  for(int i=0; i<pss.size(); i+=N) {
    Eigen::Matrix<double,3,1> pos=Eigen::Matrix<double,3,1>::Zero();
    pos.template segment<N>(0)=pss.template segment<N>(i).template cast<double>();
    vss.push_back(pos);
  }

  VTKWriter<double> os("points",path,true);
  os.appendPoints(vss.begin(),vss.end());
  VTKWriter<double>::IteratorIndex<Eigen::Matrix<int,2,1>> beg(0,0,1),end((int)vss.size()-1,0,1);
  os.appendCells(beg,end,VTKWriter<double>::LINE);
}
template <int N>
void ADMMUAVTrajectory<N>::writeNodesVTK(const std::string& path) const {
  std::vector<Eigen::Matrix<double,3,1>> vss;
  for(int i=0; i<_nodes.rows(); i++) {
    Eigen::Matrix<double,3,1> pos=Eigen::Matrix<double,3,1>::Zero();
    pos.template segment<N>(0)=_nodes.row(i).transpose().template segment<N>(0).template cast<double>();
    vss.push_back(pos);
  }

  VTKWriter<double> os("nodes",path,true);
  os.appendPoints(vss.begin(),vss.end());
  VTKWriter<double>::IteratorIndex<Eigen::Matrix<int,1,1>> beg(0,0,0),end((int)vss.size(),0,0);
  os.appendCells(beg,end,VTKWriter<double>::POINT);
}
//helper
template <int N>
bool ADMMUAVTrajectory<N>::intersect(int segId,const MeshExact& obs,T margin) const {
  if(obs.getBVH().empty())
    return false;
  //build bounding box
  BBoxExact bb;
  int order=_traj.getOrder();
  MatNXT pts=MatNXT::Zero(N,order+1);
  for(int i=0; i<order+1; i++) {
    BBoxExact::Vec3T pos=BBoxExact::Vec3T::Zero();
    pos.template segment<N>(0)=_CPs.row(segId*order+i).transpose().template cast<BBoxExact::T>();
    pts.col(i)=_CPs.row(segId*order+i).transpose();
    bb.setUnion(pos);
  }
  //dynamic small shape
  SmallShapeDynamic<N> A(pts);
  //depth first search
  std::stack<int> ss;
  ss.push((int)obs.getBVH().size()-1);
  while(!ss.empty()) {
    int r=ss.top();
    ss.pop();
    const auto& n=obs.getBVH()[r];
    if(!n._bb.intersect(bb))
      continue;
    else if(n._cell>=0) {
      MatNT ptsObs;
      for(int c=0; c<N; c++)
        ptsObs.col(c)=obs.vss()[obs.iss()[n._cell][c]].template segment<N>(0).template cast<T>();
      SmallShape<N,N> B(ptsObs);
      //calc narrow phase distance
      GJK::Vec3T pAL,pBL;
      bool intersect;
      auto dist=(T)GJK::runGJK(A,B,
                               GJK::Mat3X4T::Identity(),
                               GJK::Mat3X4T::Identity(),
                               pAL,pBL,&intersect);
      if(dist<margin)
        return true;
    } else {
      ss.push(n._l);
      ss.push(n._r);
    }
  }
  return false;
}
//instance
template class ADMMUAVTrajectory<2>;
template class ADMMUAVTrajectory<3>;
}
