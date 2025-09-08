#include "EnvironmentUtils.h"
#include "ConvexHullExact.h"
#include <Utils/RotationUtils.h>
#include <Utils/VTKWriter.h>
#include <Utils/Pragma.h>
#include <Utils/Utils.h>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <stack>

//mesh processing
namespace PHYSICSMOTION {
size_t EdgeHash::operator()(const Eigen::Matrix<int,2,1>& key) const {
  size_t seed=0;
  std::hash<int> h;
  hash_combine(seed,h(key[0]));
  hash_combine(seed,h(key[1]));
  return seed;
}
bool EdgeHash::operator()(const Eigen::Matrix<int,2,1>& a,const Eigen::Matrix<int,2,1>& b) const {
  for(int i=0; i<2; i++)
    if(a[i]<b[i])
      return true;
    else if(a[i]>b[i])
      return false;
  return false;
}
size_t TriangleHash::operator()(const Eigen::Matrix<int,3,1>& key) const {
  size_t seed=0;
  std::hash<int> h;
  hash_combine(seed,h(key[0]));
  hash_combine(seed,h(key[1]));
  hash_combine(seed,h(key[2]));
  return seed;
}
bool TriangleHash::operator()(const Eigen::Matrix<int,3,1>& a,const Eigen::Matrix<int,3,1>& b) const {
  for(int i=0; i<3; i++)
    if(a[i]<b[i])
      return true;
    else if(a[i]>b[i])
      return false;
  return false;
}
void buildVertex(const std::vector<Eigen::Matrix<int,2,1>>& iss,
                 std::unordered_map<Eigen::Matrix<int,2,1>,std::pair<int,int>,EdgeHash>& edgeMap) {
  for(int i=0; i<(int)iss.size(); i++) {
    for(int d=0; d<2; d++) {
      //insert vertex
      Eigen::Matrix<int,2,1> v(iss[i][d],-1);
      auto it=edgeMap.find(v);
      if(it==edgeMap.end())
        edgeMap[v]=std::make_pair(i,-1);
      else {
        ASSERT_MSGV(it->second.second==-1,
                    "Non-manifold 2D mesh detected, "
                    "Vertex(%d) alreadying bordering edge: "
                    "(%d,%d) and (%d,%d) when trying to border edge: "
                    "(%d,%d)!",v[0],
                    iss[it->second.first][0],iss[it->second.first][1],
                    iss[it->second.second][0],iss[it->second.second][1],
                    iss[i][0],iss[i][1])
        it->second.second=i;
      }
    }
  }
}
void buildEdge(const std::vector<Eigen::Matrix<int,3,1>>& iss,
               std::unordered_map<Eigen::Matrix<int,2,1>,std::pair<int,int>,EdgeHash>& edgeMap) {
  for(int i=0; i<(int)iss.size(); i++) {
    ASSERT_MSG(iss[i][2]!=-1,"You cannot use BottomUp BVH building in 3D Meshes!")
    for(int d=0; d<3; d++) {
      //edge index
      Eigen::Matrix<int,2,1> e(iss[i][d],iss[i][(d+1)%3]);
      if(e[0]>e[1])
        std::swap(e[0],e[1]);
      //insert edge
      auto it=edgeMap.find(e);
      if(it==edgeMap.end())
        edgeMap[e]=std::make_pair(i,-1);
      else {
        ASSERT_MSGV(it->second.second==-1,
                    "Non-manifold mesh detected, "
                    "Edge(%d,%d) alreadying bordering triangle: "
                    "(%d,%d,%d) and (%d,%d,%d) when trying to border triangle: "
                    "(%d,%d,%d)!",e[0],e[1],
                    iss[it->second.first][0],iss[it->second.first][1],iss[it->second.first][2],
                    iss[it->second.second][0],iss[it->second.second][1],iss[it->second.second][2],
                    iss[i][0],iss[i][1],iss[i][2])
        it->second.second=i;
      }
    }
  }
}
void makeUniform(std::vector<Eigen::Matrix<int,3,1>>& iss,
                 int i,int j,int v0,int v1) {
  int v0i,v1i,v0j,v1j;
  for(int d=0; d<3; d++) {
    if(iss[i][d]==v0)v0i=d;
    if(iss[i][d]==v1)v1i=d;
    if(iss[j][d]==v0)v0j=d;
    if(iss[j][d]==v1)v1j=d;
  }
  bool isI=(v0i+1)%3==v1i;
  bool isJ=(v0j+1)%3==v1j;
  if(isI==isJ)
    std::swap(iss[j][1],iss[j][2]);
}
void makeUniform(std::vector<Eigen::Matrix<int,3,1>>& iss) {
  //initialize edge
  std::unordered_map<Eigen::Matrix<int,2,1>,std::pair<int,int>,EdgeHash> eMap;
  buildEdge(iss,eMap);
  //make uniform
  std::vector<bool> visited(iss.size(),false);
  for(int i=0; i<(int)visited.size(); i++)
    if(!visited[i]) {
      std::stack<int> queue;
      queue.push(i);
      visited[i]=true;
      while(!queue.empty()) {
        int ti=queue.top();
        queue.pop();
        const Eigen::Matrix<int,3,1>& I=iss[ti];
        Eigen::Matrix<int,2,1> e;
        for(int eid=0; eid<3; eid++) {
          e[0]=(int)I[eid];
          e[1]=(int)I[(eid+1)%3];
          ASSERT(e[0]!=e[1])
          if(e[0]>e[1])
            std::swap(e[0],e[1]);
          std::pair<int,int> edg=eMap[e];
          if(edg.second>=0) {
            if(edg.second==ti)
              std::swap(edg.first,edg.second);
            if(!visited[edg.second]) {
              makeUniform(iss,edg.first,edg.second,e[0],e[1]);
              queue.push(edg.second);
              visited[edg.second]=true;
            }
          }
        }
      }
    }
}
void makeInsideOut(std::vector<Eigen::Matrix<int,3,1>>& iss) {
  for(int i=0; i<(int)iss.size(); i++) {
    std::swap(iss[i][0],iss[i][1]);
  }
}
std::shared_ptr<ConvexHullExact> makeConvexPolygon(const std::vector<Eigen::Matrix<double,2,1>>& vss,const Eigen::Matrix<double,2,1>& height) {
  std::vector<Eigen::Matrix<double,3,1>> vss3D;
  for(const Eigen::Matrix<double,2,1>& v:vss) {
    vss3D.push_back(Eigen::Matrix<double,3,1>(v[0],v[1],height[0]));
    vss3D.push_back(Eigen::Matrix<double,3,1>(v[0],v[1],height[1]));
  }
  return std::shared_ptr<ConvexHullExact>(new ConvexHullExact(vss3D));
}
std::shared_ptr<ConvexHullExact> makeRegularConvexPolygon(int n,double radius,const Eigen::Matrix<double,2,1>& pos,const Eigen::Matrix<double,2,1>& height) {
  std::vector<Eigen::Matrix<double,2,1>> vss;
  for(int i=0; i<n; i++) {
    double angle=i*M_PI*2/n;
    vss.push_back(Eigen::Matrix<double,2,1>(cos(angle)*radius+pos[0],sin(angle)*radius+pos[1]));
  }
  return makeConvexPolygon(vss,height);
}
double signedVolume(const Eigen::Matrix<double,3,1>& a,
                    const Eigen::Matrix<double,3,1>& b,
                    const Eigen::Matrix<double,3,1>& c) {
  double v321=c.x()*b.y()*a.z();
  double v231=b.x()*c.y()*a.z();
  double v312=c.x()*a.y()*b.z();
  double v132=a.x()*c.y()*b.z();
  double v213=b.x()*a.y()*c.z();
  double v123=a.x()*b.y()*c.z();
  return (1.0f/6.0f)*(-v321+v231+v312-v132-v213+v123);
}
double volume(const std::vector<Eigen::Matrix<double,3,1>>& vss,
              const std::vector<Eigen::Matrix<int,3,1>>& iss) {
  double volume=0.0;
  for(int it=0; it<(int)iss.size(); it++) {
    const Eigen::Matrix<double,3,1>& p1=vss[iss[it].x()];
    const Eigen::Matrix<double,3,1>& p2=vss[iss[it].y()];
    const Eigen::Matrix<double,3,1>& p3=vss[iss[it].z()];
    volume+=signedVolume(p1,p2,p3);
  }
  return volume;
}
void readTetMesh(std::vector<Eigen::Matrix<double,3,1>>& vss,
                 std::vector<Eigen::Matrix<int,4,1>>& ess,
                 const std::string& input) {
  vss.clear();
  ess.clear();
  std::string line;
  std::ifstream infile(input.c_str());
  std::unordered_set<int> usedVss;
  while(std::getline(infile,line))
    if(line=="Vertices") {
      std::getline(infile,line);
      for(int i=0,n=atoi(line.c_str()); i<n; i++) {
        std::getline(infile,line);
        Eigen::Matrix<double,3,1> v;
        std::istringstream(line) >> v[0] >> v[1] >> v[2];
        vss.push_back(v);
      }
    } else if(line=="Tetrahedra") {
      std::getline(infile,line);
      for(int i=0,n=atoi(line.c_str()); i<n; i++) {
        std::getline(infile,line);
        Eigen::Matrix<int,4,1> e;
        std::istringstream(line) >> e[0] >> e[1] >> e[2] >> e[3];
        e-=Eigen::Matrix<int,4,1>::Ones();
        usedVss.insert(e[0]);
        usedVss.insert(e[1]);
        usedVss.insert(e[2]);
        usedVss.insert(e[3]);
        ess.push_back(e);
      }
    }
  ASSERT_MSG(usedVss.size()==vss.size(),"Dangling vertices detected!");
}
void extractSurfaceFromTetMesh(const std::vector<Eigen::Matrix<int,4,1>>& ess,
                               std::vector<Eigen::Matrix<int,3,1>>& iss) {
  std::unordered_map<Eigen::Matrix<int,3,1>,std::pair<int,int>,TriangleHash> edgeMap;
  for(int i=0; i<(int)ess.size(); i++) {
    const auto& e=ess[i];
    for(int d=0; d<4; d++) {
      Eigen::Matrix<int,3,1> t(e[(d+1)%4],e[(d+2)%4],e[(d+3)%4]);
      sort3(t[0],t[1],t[2]);
      if(edgeMap.find(t)!=edgeMap.end())
        edgeMap.find(t)->second.second=i;
      else edgeMap[t]=std::make_pair(i,-1);
    }
  }
  for(const auto& e:edgeMap)
    if(e.second.second==-1)
      iss.push_back(e.first);
}
void writeTetVTK(const std::vector<Eigen::Matrix<double,3,1>>& vss,
                 const std::vector<Eigen::Matrix<int,4,1>>& iss,
                 const std::string& output) {
  VTKWriter<double> os("tet",output,false);
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(iss.begin(),iss.end(),VTKWriter<double>::TETRA);
}
void subdivide(std::vector<Eigen::Matrix<double,3,1>>& vss,
               std::vector<Eigen::Matrix<int,3,1>>& iss,int voff,int ioff) {
  std::unordered_map<Eigen::Matrix<int,2,1>,int,EdgeHash> umap;
  //vss
  for(int i=ioff; i<(int)iss.size(); i++)
    for(int k=0; k<3; k++) {
      Eigen::Matrix<int,2,1> eab(iss[i][(k+1)%3],iss[i][(k+2)%3]);
      sort2(eab[0],eab[1]);
      if(umap.find(eab)==umap.end()) {
        umap[eab]=(int)vss.size();
        vss.push_back((vss[eab[0]]+vss[eab[1]])/2);
      }
    }
  //iss
  std::vector<Eigen::Matrix<int,3,1>> iss0(iss.begin()+ioff,iss.end());
  iss.resize(ioff);
  for(int i=0; i<(int)iss0.size(); i++) {
    Eigen::Matrix<int,3,1> tri=iss0[i];

    Eigen::Matrix<int,2,1> e01(tri[0],tri[1]);
    sort2(e01[0],e01[1]);
    int tri01=umap.find(e01)->second;

    Eigen::Matrix<int,2,1> e02(tri[0],tri[2]);
    sort2(e02[0],e02[1]);
    int tri02=umap.find(e02)->second;

    Eigen::Matrix<int,2,1> e12(tri[1],tri[2]);
    sort2(e12[0],e12[1]);
    int tri12=umap.find(e12)->second;

    iss.push_back(Eigen::Matrix<int,3,1>(tri[0],tri01,tri02));
    iss.push_back(Eigen::Matrix<int,3,1>(tri[1],tri12,tri01));
    iss.push_back(Eigen::Matrix<int,3,1>(tri[2],tri02,tri12));
    iss.push_back(Eigen::Matrix<int,3,1>(tri01,tri12,tri02));
  }
}
std::pair<int,int> addBox(std::vector<Eigen::Matrix<double,3,1>>& vss,
                          std::vector<Eigen::Matrix<int,3,1>>& iss,
                          Eigen::Matrix<double,3,1> minC,Eigen::Matrix<double,3,1> maxC) {
  int voff=(int)vss.size();
  int ioff=(int)iss.size();
  for(int i=0; i<8; i++)
    vss.push_back(Eigen::Matrix<double,3,1>((i&1)?maxC[0]:minC[0],
                                            (i&2)?maxC[1]:minC[1],
                                            (i&4)?maxC[2]:minC[2]));
  iss.push_back(Eigen::Matrix<int,3,1>(0,1,5)+Eigen::Matrix<int,3,1>::Constant(voff));
  iss.push_back(Eigen::Matrix<int,3,1>(0,5,4)+Eigen::Matrix<int,3,1>::Constant(voff));
  iss.push_back(Eigen::Matrix<int,3,1>(3,2,6)+Eigen::Matrix<int,3,1>::Constant(voff));
  iss.push_back(Eigen::Matrix<int,3,1>(3,6,7)+Eigen::Matrix<int,3,1>::Constant(voff));
  iss.push_back(Eigen::Matrix<int,3,1>(1,3,7)+Eigen::Matrix<int,3,1>::Constant(voff));
  iss.push_back(Eigen::Matrix<int,3,1>(1,7,5)+Eigen::Matrix<int,3,1>::Constant(voff));
  iss.push_back(Eigen::Matrix<int,3,1>(0,4,6)+Eigen::Matrix<int,3,1>::Constant(voff));
  iss.push_back(Eigen::Matrix<int,3,1>(0,6,2)+Eigen::Matrix<int,3,1>::Constant(voff));
  iss.push_back(Eigen::Matrix<int,3,1>(0,2,3)+Eigen::Matrix<int,3,1>::Constant(voff));
  iss.push_back(Eigen::Matrix<int,3,1>(0,3,1)+Eigen::Matrix<int,3,1>::Constant(voff));
  iss.push_back(Eigen::Matrix<int,3,1>(4,5,7)+Eigen::Matrix<int,3,1>::Constant(voff));
  iss.push_back(Eigen::Matrix<int,3,1>(4,7,6)+Eigen::Matrix<int,3,1>::Constant(voff));
  return std::make_pair(voff,ioff);
}
std::pair<int,int> addSphere(std::vector<Eigen::Matrix<double,3,1>>& vss,
                             std::vector<Eigen::Matrix<int,3,1>>& iss,
                             Eigen::Matrix<double,3,1> pos,double radius,int res) {
  auto vioff=addBox(vss,iss,Eigen::Matrix<double,3,1>::Constant(-1),Eigen::Matrix<double,3,1>::Constant(1));
  for(int i=0; i<res; i++)
    subdivide(vss,iss,vioff.first,vioff.second);
  for(int i=vioff.first; i<(int)vss.size(); i++)
    vss[i]=pos+vss[i].normalized()*radius;
  return vioff;
}
std::pair<int,int> addCircle(std::vector<Eigen::Matrix<double,3,1>>& vss,
                             std::vector<Eigen::Matrix<int,3,1>>& iss,
                             Eigen::Matrix<double,3,1> pos,double radius,int res) {
  int voff=(int)vss.size();
  int ioff=(int)iss.size();
  for(int i=0; i<res; i++) {
    double angle=(M_PI*2*i)/res;
    vss.push_back(Eigen::Matrix<double,3,1>(cos(angle),sin(angle),0)*radius+pos);
  }
  vss.push_back(pos);
  for(int i=0; i<res; i++)
    iss.push_back(Eigen::Matrix<int,3,1>(i,(i+1)%res,res)+Eigen::Matrix<int,3,1>::Constant(voff));
  return std::make_pair(voff,ioff);
}
std::pair<int,int> addCapsule(std::vector<Eigen::Matrix<double,3,1>>& vss,
                              std::vector<Eigen::Matrix<int,3,1>>& iss,
                              Eigen::Matrix<double,3,1> minC,Eigen::Matrix<double,3,1> maxC,double radius,int res) {
  double dist=(maxC-minC).norm();
  Eigen::Matrix<double,3,1> ctr=(maxC+minC)/2,w;
  Eigen::Matrix<double,3,1> dir=(maxC-minC).normalized();
  if((dir-Eigen::Matrix<double,3,1>::UnitZ()).norm()<1e-6 ||
      (dir+Eigen::Matrix<double,3,1>::UnitZ()).norm()<1e-6)
    w.setZero();
  else {
    w=Eigen::Matrix<double,3,1>::UnitZ().cross(dir).normalized();
    w*=std::acos(Eigen::Matrix<double,3,1>::UnitZ().dot(dir));
  }
  Eigen::Matrix<double,3,3> R=expWGradV<double,Eigen::Matrix<double,3,1>>(w);
  auto vioff=addSphere(vss,iss,Eigen::Matrix<double,3,1>(0,0,0),radius,res);
  for(int i=vioff.first; i<(int)vss.size(); i++) {
    if(vss[i][2]>0)
      vss[i][2]+=dist/2;
    else if(vss[i][2]<0)
      vss[i][2]-=dist/2;
    vss[i]=R*vss[i]+ctr;
  }
  return vioff;
}
void extrude(std::vector<Eigen::Matrix<double,3,1>>& vss,
             std::vector<Eigen::Matrix<int,3,1>>& iss,double halfThickness) {
  std::vector<Eigen::Matrix<double,3,1>> vssE;
  std::vector<Eigen::Matrix<int,3,1>> issE;
  //front and back piece
  for(auto v:vss)
    vssE.push_back(v-Eigen::Matrix<double,3,1>::UnitZ()*halfThickness);
  for(auto v:vss)
    vssE.push_back(v+Eigen::Matrix<double,3,1>::UnitZ()*halfThickness);
  for(auto i:iss)
    issE.push_back(i);
  for(auto i:iss)
    issE.push_back((i.array()+(int)vss.size()).matrix());
  //build edge
  std::unordered_map<Eigen::Matrix<int,2,1>,std::pair<int,int>,EdgeHash> edgeMap;
  buildEdge(iss,edgeMap);
  //vertical piece
  for(auto e:edgeMap)
    if(e.second.second==-1) {
      int v0=e.first[0];
      int v1=e.first[1];
      int v3=e.first[0]+(int)vss.size();
      int v2=e.first[1]+(int)vss.size();
      issE.push_back(Eigen::Matrix<int,3,1>(v0,v1,v2));
      issE.push_back(Eigen::Matrix<int,3,1>(v0,v2,v3));
    }
  vssE.swap(vss);
  issE.swap(iss);
}
}
