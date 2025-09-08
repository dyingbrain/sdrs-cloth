#ifndef ENVIRONMENT_UTILS_H
#define ENVIRONMENT_UTILS_H

#include <Eigen/Dense>
#include <unordered_map>
#include <vector>
#include <memory>

namespace PHYSICSMOTION {
struct EdgeHash {
  size_t operator()(const Eigen::Matrix<int,2,1>& key) const;
  bool operator()(const Eigen::Matrix<int,2,1>& a, const Eigen::Matrix<int,2,1>& b) const;
};
struct TriangleHash {
  size_t operator()(const Eigen::Matrix<int,3,1>& key) const;
  bool operator()(const Eigen::Matrix<int,3,1>& a, const Eigen::Matrix<int,3,1>& b) const;
};
struct ConvexHullExact;
extern void buildVertex(const std::vector<Eigen::Matrix<int,2,1>>& iss,
                        std::unordered_map<Eigen::Matrix<int,2,1>,std::pair<int,int>,EdgeHash>& edgeMap);
extern void buildEdge(const std::vector<Eigen::Matrix<int,3,1>>& iss,
                      std::unordered_map<Eigen::Matrix<int,2,1>,std::pair<int,int>,EdgeHash>& edgeMap);
extern void makeUniform(std::vector<Eigen::Matrix<int,3,1>>& iss);
extern void makeInsideOut(std::vector<Eigen::Matrix<int,3,1>>& iss);
extern std::shared_ptr<ConvexHullExact> makeConvexPolygon(const std::vector<Eigen::Matrix<double,3,1>>& vss,const Eigen::Matrix<double,2,1>& height);
extern std::shared_ptr<ConvexHullExact> makeRegularConvexPolygon(int n,double radius,const Eigen::Matrix<double,2,1>& pos,const Eigen::Matrix<double,2,1>& height);
extern double volume(const std::vector<Eigen::Matrix<double,3,1>>& vss,
                     const std::vector<Eigen::Matrix<int,3,1>>& iss);
extern void readTetMesh(std::vector<Eigen::Matrix<double,3,1>>& vss,
                        std::vector<Eigen::Matrix<int,4,1>>& iss,
                        const std::string& input);
extern void extractSurfaceFromTetMesh(const std::vector<Eigen::Matrix<int,4,1>>& ess,
                                      std::vector<Eigen::Matrix<int,3,1>>& iss);
extern void writeTetVTK(const std::vector<Eigen::Matrix<double,3,1>>& vss,
                        const std::vector<Eigen::Matrix<int,4,1>>& iss,
                        const std::string& output);
extern void makeConvex(std::vector<Eigen::Matrix<double,3,1>>& vss,
                       std::vector<Eigen::Matrix<int,3,1>>& iss);
extern void makeConvexProject(std::vector<Eigen::Matrix<double,3,1>>& vss);
extern void generateTetMesh(const std::vector<Eigen::Matrix<double,3,1>>& vss,
                            const std::vector<Eigen::Matrix<int,3,1>>& iss,
                            const std::string& output,double feature_angle=91,
                            double size=1,double dist=1,
                            double angle=25,double ratio=3);
extern void generateTriMesh(std::vector<Eigen::Matrix<double,3,1>>& vss,
                            std::vector<Eigen::Matrix<int,3,1>>& iss,
                            const std::vector<Eigen::Matrix<double,2,1>>& sss,
                            double size=1.,double aspect=.125);
extern std::pair<int,int> addBox(std::vector<Eigen::Matrix<double,3,1>>& vss,
                                 std::vector<Eigen::Matrix<int,3,1>>& iss,
                                 Eigen::Matrix<double,3,1> minC,Eigen::Matrix<double,3,1> maxC);
extern std::pair<int,int> addSphere(std::vector<Eigen::Matrix<double,3,1>>& vss,
                                    std::vector<Eigen::Matrix<int,3,1>>& iss,
                                    Eigen::Matrix<double,3,1> pos,double radius,int res=2);
extern std::pair<int,int> addCircle(std::vector<Eigen::Matrix<double,3,1>>& vss,
                                    std::vector<Eigen::Matrix<int,3,1>>& iss,
                                    Eigen::Matrix<double,3,1> pos,double radius,int res=32);
extern std::pair<int,int> addCapsule(std::vector<Eigen::Matrix<double,3,1>>& vss,
                                     std::vector<Eigen::Matrix<int,3,1>>& iss,
                                     Eigen::Matrix<double,3,1> minC,Eigen::Matrix<double,3,1> maxC,double radius,int res=2);
extern void extrude(std::vector<Eigen::Matrix<double,3,1>>& vss,
                    std::vector<Eigen::Matrix<int,3,1>>& iss,double halfThickness);
}

#endif
