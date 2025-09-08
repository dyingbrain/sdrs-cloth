#include <ADMM/SpatialHash.h>
#include <Utils/Epsilon.h>
#include <unordered_set>
#include <chrono>

using namespace PHYSICSMOTION;
typedef FLOAT T;
DECL_MAT_VEC_MAP_TYPES_T

void createGrid(int N,MeshExact& m,T rand=0) {
#define ID(I,J) (I)*(N+1)+(J)
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,3,1>> iss;
  for(int i=0; i<=N; i++)
    for(int j=0; j<=N; j++) {
      vss.push_back(Eigen::Matrix<double,3,1>(i,j,0));
      if(rand>0)
        vss.back()+=Eigen::Matrix<double,3,1>::Random()*(double)rand;
    }
  for(int i=0; i<N; i++)
    for(int j=0; j<N; j++) {
      iss.push_back(Eigen::Matrix<int,3,1>(ID(i,j),ID(i+1,j),ID(i+1,j+1)));
      iss.push_back(Eigen::Matrix<int,3,1>(ID(i,j),ID(i+1,j+1),ID(i,j+1)));
    }
  m.init<double>(vss,iss);
#undef ID
}
template <int N,int RES=32>
void debugCollision(int m,T eps,T rand) {
  typedef std::chrono::high_resolution_clock Clock;
  //build mesh grid
  MeshExact mesh,obs;
  createGrid(RES,mesh);
  createGrid(RES,obs,rand);
  typedef CollisionDetector<N,1> Detector;
  typedef typename Detector::ID ID;
  std::cout << "DebugSpatialHash" << " eps=" << eps << std::endl;
  Vec x0,x,x2;
  Detector detector(obs);
  detector.extractPos(mesh,x0);
  for(int i=0; i<m; i++) {
    x.setRandom(x0.size());
    x*=rand;
    x+=x0;
    x2.setRandom(x0.size());
    x2*=rand;
    x2+=x0;
    {
      //serial
      {
        auto beg=Clock::now();
        detector.generateBroadBF(x,x2,eps,eps/2,true,true);
        std::cout << "BroadBF " << std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now()-beg).count();
      }
      std::unordered_set<ID> idSelfA,idObsA;
      for(ID d:detector.selfCollisions())
        if(d!=Detector::INVALID_ID)
          idSelfA.insert(d);
      for(ID d:detector.obsCollisions())
        if(d!=Detector::INVALID_ID)
          idObsA.insert(d);
      //parallel
      {
        auto beg=Clock::now();
        detector.generateBroad(x,x2,eps,eps/2,true,true);
        std::cout << " Broad " << std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now()-beg).count() << std::endl;
      }
      std::unordered_set<ID> idSelfB,idObsB;
      for(ID d:detector.selfCollisions())
        if(d!=Detector::INVALID_ID)
          idSelfB.insert(d);
      for(ID d:detector.obsCollisions())
        if(d!=Detector::INVALID_ID)
          idObsB.insert(d);
      std::cout << "Broadphase: nSelf=" << idSelfA.size() << " nObs=" << idObsA.size() << std::endl;
      //check self
      for(ID d:idSelfA) {
        ASSERT_MSG(idSelfB.find(d)!=idSelfB.end(),"Error detected!")
      }
      for(ID d:idSelfB) {
        ASSERT_MSG(idSelfA.find(d)!=idSelfA.end(),"Error detected!")
      }
      //check obs
      for(ID d:idObsA) {
        ASSERT_MSG(idObsB.find(d)!=idObsB.end(),"Error detected!")
      }
      for(ID d:idObsB) {
        ASSERT_MSG(idObsA.find(d)!=idObsA.end(),"Error detected!")
      }
    }
  }
}
int main(int argc,char** argv) {
  debugCollision<2>(10,0.1f,0.5f);
  debugCollision<3>(10,0.1f,0.5f);
  return 0;
}
