#include <Environment/GJK.h>
#include <Utils/VTKWriter.h>
#include <Utils/IO.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  Mat3T a,b;
  {
    std::ifstream is("gjk.dat",std::ios::binary);
    readBinaryData(a,is);
    readBinaryData(b,is);
  }

  //create two meshes
  std::shared_ptr<MeshExact> A(new MeshExact);
  A->init<double>(std::vector<Eigen::Matrix<double,3,1>>({
      Eigen::Matrix<double,3,1>(a.col(0)),
      Eigen::Matrix<double,3,1>(a.col(1)),
      Eigen::Matrix<double,3,1>(a.col(2))}), {Eigen::Matrix<int,3,1>(0,1,2)},false);
  std::shared_ptr<MeshExact> B(new MeshExact);
  B->init<double>(std::vector<Eigen::Matrix<double,3,1>>({
      Eigen::Matrix<double,3,1>(b.col(0)),
      Eigen::Matrix<double,3,1>(b.col(1)),
      Eigen::Matrix<double,3,1>(b.col(2))}), {Eigen::Matrix<int,3,1>(0,1,2)},false);

  //run GJK
  GJK::Mat3X4T transA,transB;
  transA.setIdentity();
  transB.setIdentity();
  GJK::Vec3T ptAL,ptBL;
  bool intersect;
  GJK::T dist=GJK::runGJK(A,B,transA,transB,ptAL,ptBL,&intersect);
  GJK::Vec3T ptA=transA.template block<3,3>(0,0)*ptAL+transA.col(3);
  GJK::Vec3T ptB=transB.template block<3,3>(0,0)*ptBL+transB.col(3);
  std::cout << "dist=" << dist << " intersect=" <<intersect << std::endl;
  {
    GJK::Vec3T n=(ptA-ptB).normalized();
    GJK::Vec2T a=A->project(transA.template block<3,3>(0,0).transpose()*n)+GJK::Vec2T::Constant(n.dot(transA.col(3)));
    GJK::Vec2T b=B->project(transB.template block<3,3>(0,0).transpose()*n)+GJK::Vec2T::Constant(n.dot(transB.col(3)));
    bool intersect=a[1]>=b[0] && b[1]>=a[0];
    std::cout << "IntervalA=" << a.transpose() << " IntervalB=" << b.transpose() << " IntervalIntersect=" << intersect << std::endl;
  }

  //write result:mesh
  {
    VTKWriter<double> os("GJK","GJKResult.vtk",true);
    A->writeVTK(os,transA);
    B->writeVTK(os,transB);
  }
  //write result:line
  {
    VTKWriter<double> os("GJK","GJKVector.vtk",true);
    std::vector<Eigen::Matrix<double,3,1>> pss;
    std::vector<Eigen::Matrix<int,2,1>> iss({Eigen::Matrix<int,2,1>(0,1)});
    pss.push_back(ptA.template cast<double>());
    pss.push_back(ptB.template cast<double>());
    os.appendPoints(pss.begin(),pss.end());
    os.appendCells(iss.begin(),iss.end(),VTKWriter<double>::LINE);
  }
  return 0;
}
