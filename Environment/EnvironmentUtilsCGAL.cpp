#define _USE_MATH_DEFINES
#include <cmath>
#include "EnvironmentUtils.h"
#include <Utils/Pragma.h>
#include <list>

//make convex
#ifdef CGAL_SUPPORT
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
//tet generation
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>
//convex hull generation
#include <CGAL/convex_hull_2.h>
#include <CGAL/convex_hull_3.h>
//polyhedron
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
//constrained delaunay triangulation
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/draw_constrained_triangulation_2.h>
namespace PHYSICSMOTION {
struct PointHash {
  size_t operator()(const Eigen::Matrix<double,3,1>& key) const {
    size_t seed=0;
    std::hash<double> h;
    boost::hash_combine(seed,h(key[0]));
    boost::hash_combine(seed,h(key[1]));
    boost::hash_combine(seed,h(key[2]));
    return seed;
  }
};
template <class HDS>
class BuildMesh : public CGAL::Modifier_base<HDS> {
 public:
  BuildMesh(const std::vector<Eigen::Matrix<double,3,1>>& vss,
            const std::vector<Eigen::Matrix<int,3,1>>& iss):_vss(vss),_iss(iss) {}
  void operator()(HDS& hds) {
    // Postcondition: hds is a valid polyhedral surface.
    CGAL::Polyhedron_incremental_builder_3<HDS> B(hds,true);
    B.begin_surface(_vss.size(),_iss.size());
    typedef typename HDS::Vertex Vertex;
    typedef typename Vertex::Point Point;
    for(auto& v:_vss)
      B.add_vertex(Point(v[0],v[1],v[2]));
    for(auto& t:_iss) {
      B.begin_facet();
      B.add_vertex_to_facet(t[0]);
      B.add_vertex_to_facet(t[1]);
      B.add_vertex_to_facet(t[2]);
      B.end_facet();
    }
    B.end_surface();
  }
  const std::vector<Eigen::Matrix<double,3,1>>& _vss;
  const std::vector<Eigen::Matrix<int,3,1>>& _iss;
};
void makeConvex(std::vector<Eigen::Matrix<double,3,1>>& vss,
                std::vector<Eigen::Matrix<int,3,1>>& iss) {
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Polyhedron_3<K> Polyhedron_3;
  typedef Polyhedron_3::Vertex_iterator Vertex_iterator;
  typedef Polyhedron_3::Facet_iterator Facet_iterator;
  typedef Polyhedron_3::Halfedge_around_facet_circulator Halfedge_facet_circulator;
  typedef K::Point_3 Point_3;

  std::vector<Point_3> points;
  for(int i=0; i<(int)vss.size(); i++)
    points.push_back(Point_3(vss[i][0],vss[i][1],vss[i][2]));

  vss.clear();
  iss.clear();
  Polyhedron_3 poly;
  CGAL::convex_hull_3(points.begin(),points.end(),poly);
  for(Vertex_iterator v=poly.vertices_begin(); v!=poly.vertices_end(); ++v)
    vss.push_back(Eigen::Matrix<double,3,1>(v->point()[0],v->point()[1],v->point()[2]));
  for(Facet_iterator i=poly.facets_begin(); i!=poly.facets_end(); ++i) {
    Halfedge_facet_circulator j=i->facet_begin();
    //CGAL_assertion(CGAL::circulator_size(j)>=3);
    std::vector<int> facet;
    do {
      facet.push_back(std::distance(poly.vertices_begin(),j->vertex()));
    } while(++j!=i->facet_begin());
    for(auto t=0; t<(int)facet.size()-2; t++)
      iss.push_back(Eigen::Matrix<int,3,1>(facet[t],facet[t+1],facet[t+2]));
  }
}
void makeConvexProject(std::vector<Eigen::Matrix<double,3,1>>& vss) {
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef K::Point_2 Point_2;

  std::vector<Point_2> points,pointsOut;
  for(int i=0; i<(int)vss.size(); i++)
    points.push_back(Point_2(vss[i][0],vss[i][1]));

  vss.clear();
  CGAL::convex_hull_2(points.begin(),points.end(),std::back_inserter(pointsOut));
  for(int i=0; i<(int)pointsOut.size(); i++)
    vss.push_back(Eigen::Matrix<double,3,1>(pointsOut[i].x(),pointsOut[i].y(),0));
}
void generateTetMesh(const std::vector<Eigen::Matrix<double,3,1>>& vss,
                     const std::vector<Eigen::Matrix<int,3,1>>& iss,
                     const std::string& output,double feature_angle,
                     double size,double dist,double angle,double ratio) {
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron_3;
  typedef Polyhedron_3::HalfedgeDS HalfedgeDS;

  typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;
  typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
  typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

  typedef K::Point_3 Point_3;
  typedef std::vector<Point_3> Polyline_3;
  typedef std::list<Polyline_3> Polylines;

  using namespace CGAL::parameters;

  //build mesh
  Polyhedron_3 poly;
  BuildMesh<HalfedgeDS> tmesh(vss,iss);
  poly.delegate(tmesh);
  Mesh_domain domain(poly);

  //detect feature
  Polylines polylines;
  std::unordered_map<Eigen::Matrix<int,2,1>,std::pair<int,int>,EdgeHash> edgeMap;
  auto dihedralAngle=[&](Eigen::Matrix<int,2,1> v,std::pair<int,int> t) {
    Eigen::Matrix<double,3,1> v0=vss[v[0]],v1=vss[v[1]],d0=(v1-v0).normalized(),d1,d2;
    //first vertex
    const auto& t1=iss[t.first];
    for(int d=0; d<3; d++)
      if(t1[d]!=v[0] && t1[d]!=v[1])
        d1=vss[t1[d]];
    d1-=v0;
    d1-=d1.dot(d0)*d0;
    //second vertex
    const auto& t2=iss[t.second];
    for(int d=0; d<3; d++)
      if(t2[d]!=v[0] && t2[d]!=v[1])
        d2=vss[t2[d]];
    d2-=v0;
    d2-=d2.dot(d0)*d0;
    //return angle
    double cosVal=d1.dot(d2)/d1.norm()/d2.norm();
    double angle=std::acos(std::min<double>(std::max<double>(cosVal,-1),1))*180/M_PI;
    return angle;
  };
  buildEdge(iss,edgeMap);
  for(const auto& e:edgeMap)
    if(e.second.second==-1 || dihedralAngle(e.first,e.second)<feature_angle) {
      Polyline_3 p;
      const auto& v1=vss[e.first[0]];
      const auto& v2=vss[e.first[1]];
      p.push_back(Point_3(v1[0],v1[1],v1[2]));
      p.push_back(Point_3(v2[0],v2[1],v2[2]));
      polylines.push_back(p);
    }
  domain.add_features(polylines.begin(),polylines.end());

  //generate tet
  Mesh_criteria criteria(facet_angle=angle,
                         facet_size=size,
                         facet_distance=dist,
                         cell_radius_edge_ratio=ratio);
  C3t3 c3t3=CGAL::make_mesh_3<C3t3>(domain,criteria);//,no_perturb(),no_exude());

  //output
  std::ofstream file(output.c_str());
  c3t3.output_to_medit(file);
}
void generateTriMesh(std::vector<Eigen::Matrix<double,3,1>>& vss,
                     std::vector<Eigen::Matrix<int,3,1>>& iss,
                     const std::vector<Eigen::Matrix<double,2,1>>& sss,
                     double size,double aspect) {
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Triangulation_vertex_base_2<K> Vb;
  typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
  typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;

  typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
  typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher;

  typedef CDT::Vertex_handle Vertex_handle;
  typedef CDT::Point Point;

  using namespace CGAL::parameters;

  //insert edge
  CDT cdt;
  std::vector<Vertex_handle> vhs;
  for(int i=0; i<(int)vss.size(); i++)
    vhs.push_back(cdt.insert(Point(vss[i][0],vss[i][1])));
  for(int i=0; i<(int)iss.size(); i++)
    cdt.insert_constraint(vhs[iss[i][0]],vhs[iss[i][1]]);

  //insert seeds
  std::list<Point> list_of_seeds;
  for(int i=0; i<(int)sss.size(); i++)
    list_of_seeds.push_back(Point(sss[i][0],sss[i][1]));

  //triangulate
  CGAL::refine_Delaunay_mesh_2(cdt, seeds(list_of_seeds));
  Mesher mesher(cdt);
  mesher.set_criteria(Criteria(aspect, size));
  mesher.refine_mesh();

  //read back
  int offV=0;
  iss.clear();
  Eigen::Matrix<int,3,1> t;
  std::unordered_map<Eigen::Matrix<double,3,1>,int,PointHash> handles;
  for(CDT::Finite_faces_iterator fit=cdt.finite_faces_begin(); fit!=cdt.finite_faces_end(); ++fit) {
    if(!fit->is_in_domain())
      continue;
    for(int d=0; d<3; d++) {
      auto pt=fit->vertex(d)->point();
      Eigen::Matrix<double,3,1> ptd(pt[0],pt[1],0);
      if(handles.find(ptd)==handles.end())
        handles[ptd]=offV++;
      t[d]=handles[ptd];
    }
    iss.push_back(t);
  }
  vss.resize(offV);
  for(const auto& vhandle:handles)
    vss[vhandle.second]=vhandle.first;
}
}
#else
namespace PHYSICSMOTION {
void makeConvex(std::vector<Eigen::Matrix<double,3,1>>& vss,
                std::vector<Eigen::Matrix<int,3,1>>& iss) {
  FUNCTION_NOT_IMPLEMENTED
}
void makeConvexProject(std::vector<Eigen::Matrix<double,3,1>>& vss) {
  FUNCTION_NOT_IMPLEMENTED
}
void generateTetMesh(const std::vector<Eigen::Matrix<double,3,1>>& vss,
                     const std::vector<Eigen::Matrix<int,3,1>>& iss,
                     const std::string& output,double feature_angle,
                     double size,double dist,double angle,double ratio) {
  FUNCTION_NOT_IMPLEMENTED
}
void generateTriMesh(std::vector<Eigen::Matrix<double,3,1>>& vss,
                     std::vector<Eigen::Matrix<int,3,1>>& iss,
                     const std::vector<Eigen::Matrix<double,2,1>>& sss,
                     double size,double aspect) {
  FUNCTION_NOT_IMPLEMENTED
}
}
#endif
