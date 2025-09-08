#ifndef DEFORMABLE_VISUALIZER_H
#define DEFORMABLE_VISUALIZER_H

#include "Deformable.h"
#include <Environment/EnvironmentVisualizer.h>
#include <TinyVisualizer/MakeMesh.h>
#include <TinyVisualizer/MeshShape.h>
#include <TinyVisualizer/CompositeShape.h>

namespace PHYSICSMOTION {
void setVertices(const Deformable<3>& def,std::shared_ptr<DRAWER::MeshShape> m,bool wire=true) {
  std::vector<GLfloat> vss;
  vss.resize(def.x().size());
  Eigen::Map<Eigen::Matrix<GLfloat,-1,1>>(vss.data(),vss.size())=def.x().template cast<GLfloat>();
  if(m->nrVertex()==0)
    for(int i=0; i<(int)vss.size(); i+=3)
      m->addVertex(Eigen::Matrix<GLfloat,3,1>(vss[i],vss[i+1],vss[i+2]));
  else m->setVertices(vss);
  if(!wire)
    m->computeNormals();
}
void setVerticesSubd(const Deformable<3>& def,std::shared_ptr<DRAWER::MeshShape> m,bool wire=true) {
  std::vector<GLfloat> vss;
  Deformable<3>::Vec x=def.getTrajSubd()*def.x();
  vss.resize(x.size());
  Eigen::Map<Eigen::Matrix<GLfloat,-1,1>>(vss.data(),vss.size())=x.template cast<GLfloat>();
  if(m->nrVertex()==0)
    for(int i=0; i<(int)vss.size(); i+=3)
      m->addVertex(Eigen::Matrix<GLfloat,3,1>(vss[i],vss[i+1],vss[i+2]));
  else m->setVertices(vss);
}
void setVertices(const Deformable<2>& def,std::shared_ptr<DRAWER::MeshShape> m,bool wire=true) {
  std::vector<GLfloat> vss;
  vss.resize(def.x().size()/2*3);
  Eigen::Map<Eigen::Matrix<GLfloat,3,-1>> vssM(vss.data(),3,def.x().size()/2);
  vssM.setZero();
  vssM.block(0,0,2,def.x().size()/2)=def.x().reshaped(2,def.x().size()/2).template cast<GLfloat>();
  if(m->nrVertex()==0)
    for(int i=0; i<(int)vss.size(); i+=3)
      m->addVertex(Eigen::Matrix<GLfloat,3,1>(vss[i],vss[i+1],vss[i+2]));
  else m->setVertices(vss);
}
void setVerticesSubd(const Deformable<2>& def,std::shared_ptr<DRAWER::MeshShape> m,bool wire=true) {
  std::vector<GLfloat> vss;
  Deformable<2>::Vec x=def.getTrajSubd()*def.x();
  vss.resize(x.size()/2*3);
  Eigen::Map<Eigen::Matrix<GLfloat,3,-1>> vssM(vss.data(),3,x.size()/2);
  vssM.setZero();
  vssM.block(0,0,2,x.size()/2)=x.reshaped(2,x.size()/2).template cast<GLfloat>();
  if(m->nrVertex()==0)
    for(int i=0; i<(int)vss.size(); i+=3)
      m->addVertex(Eigen::Matrix<GLfloat,3,1>(vss[i],vss[i+1],vss[i+2]));
  else m->setVertices(vss);
}
template <int N>
std::shared_ptr<DRAWER::CompositeShape> visualizeDeformable(Deformable<N>& def,GLfloat width=5,int UAVRes=32) {
  std::shared_ptr<DRAWER::CompositeShape> ret(new DRAWER::CompositeShape);
  //traj/agent
  if(def.getUAVTraj().numCP()>0) {
    def.setUAVTrajResolution(UAVRes);
    std::shared_ptr<DRAWER::MeshShape> m(new DRAWER::MeshShape);
    setVerticesSubd(def,m);
    int n=def.getTrajSubd().rows()/N;
    for(int i=0; i<n-1; i++) {
      m->addIndexSingle(i+0);
      m->addIndexSingle(i+1);
    }
    m->setMode(GL_LINES);
    m->setLineWidth(width);
    m->setColorDiffuse(GL_LINES,1,.5f,0);
    ret->addShape(m);
  } else if(def.r()>0) {
    std::shared_ptr<DRAWER::MeshShape> m;
    if(N==2)
      m=DRAWER::makeCircle(16,true,Eigen::Matrix<GLfloat,2,1>::Zero(),def.actualR());
    else m=DRAWER::makeSphere(16,true,def.actualR());
    std::shared_ptr<DRAWER::CompositeShape> s(new DRAWER::CompositeShape);
    Eigen::Matrix<GLfloat,3,1> pos=Eigen::Matrix<GLfloat,3,1>::Zero();
    for(int i=0; i<def.x().size(); i+=N) {
      std::shared_ptr<DRAWER::Bullet3DShape> b(new DRAWER::Bullet3DShape);
      pos.template segment<N>(0)=def.x().template segment<N>(i).template cast<GLfloat>();
      b->setLocalTranslate(pos);
      b->addShape(m);
      s->addShape(b);
    }
    ret->addShape(s);
  }
  //edge
  if(!def.getEss().empty()) {
    std::shared_ptr<DRAWER::MeshShape> m(new DRAWER::MeshShape);
    setVertices(def,m);
    for(const auto& e:def.getEss()) {
      m->addIndexSingle(e.first[0]);
      m->addIndexSingle(e.first[1]);
    }
    m->setMode(GL_LINES);
    m->setLineWidth(width);
    m->setColorDiffuse(GL_LINES,1,.5f,0);
    ret->addShape(m);
  }
  //element
  if((N==2 && !def.getTss().empty()) || (N==3 && !def.getSss().empty())) {
    std::shared_ptr<DRAWER::MeshShape> m(new DRAWER::MeshShape);
    setVertices(def,m);
    if(N==2)
      for(const auto& e:def.getTss())
        m->addIndex(e._indices.template cast<GLuint>());
    else
      for(const auto& t:def.getSss())
        m->addIndex(t.template cast<GLuint>());
    m->setMode(GL_TRIANGLES);
    m->computeNormals();
    ret->addShape(m);
  }
  //fix
  if(!def.getCons()._fixes.empty()) {
    int id=0;
    std::shared_ptr<DRAWER::MeshShape> m(new DRAWER::MeshShape);
    for(int j=0; j<(int)def.getCons()._fixes.size(); j++)
      for(int i=0; i<2; i++) {
        m->addVertex(Eigen::Matrix<GLfloat,3,1>::Zero());
        m->addIndexSingle(id++);
      }
    m->setMode(GL_LINES);
    m->setLineWidth(width);
    m->setColorAmbient(GL_LINES,1,0,0);
    m->setColorDiffuse(GL_LINES,1,0,0);
    ret->addShape(m);
  }
  //l1
  if(!def.getL1ss().empty()) {
    std::shared_ptr<DRAWER::MeshShape> m(new DRAWER::MeshShape);
    setVertices(def,m);
    for(const auto& e:def.getL1ss()) {
      m->addIndexSingle(e.first[0]);
      m->addIndexSingle(e.first[1]);
    }
    m->setMode(GL_LINES);
    m->setLineWidth(width);
    m->setColorDiffuse(GL_LINES,.5f,0,0);
    ret->addShape(m);
  }
  //obs
  if(!def.getObs().iss().empty()) {
    std::shared_ptr<DRAWER::Shape> m=visualizeMeshExact(def.getObs(),N==2);
    m->setLineWidth(width);
    m->setColorDiffuse(GL_LINES,.5f,.5f,.5f);
    ret->addShape(m);
  }
  return ret;
}
template <int N>
void updateDeformable(std::shared_ptr<DRAWER::CompositeShape> ret,const Deformable<N>& def) {
  std::shared_ptr<DRAWER::MeshShape> m;
  int id=0;
  //traj/agent
  if(def.getUAVTraj().numCP()>0) {
    m=std::dynamic_pointer_cast<DRAWER::MeshShape>(ret->getChild(id++));
    setVerticesSubd(def,m);
  } else if(def.r()>0) {
    std::shared_ptr<DRAWER::CompositeShape> s=std::dynamic_pointer_cast<DRAWER::CompositeShape>(ret->getChild(id++));
    Eigen::Matrix<GLfloat,3,1> pos=Eigen::Matrix<GLfloat,3,1>::Zero();
    for(int i=0,j=0; i<def.x().size(); i+=N,j++) {
      pos.template segment<N>(0)=def.x().template segment<N>(i).template cast<GLfloat>();
      std::shared_ptr<DRAWER::Bullet3DShape> b=std::dynamic_pointer_cast<DRAWER::Bullet3DShape>(s->getChild(j));
      b->setLocalTranslate(pos);
    }
  }
  //edge
  if(!def.getEss().empty()) {
    m=std::dynamic_pointer_cast<DRAWER::MeshShape>(ret->getChild(id++));
    setVertices(def,m);
  }
  //element
  if((N==2 && !def.getTss().empty()) || (N==3 && !def.getSss().empty())) {
    m=std::dynamic_pointer_cast<DRAWER::MeshShape>(ret->getChild(id++));
    setVertices(def,m,false);
  }
  //fix
  if(!def.getCons()._fixes.empty()) {
    m=std::dynamic_pointer_cast<DRAWER::MeshShape>(ret->getChild(id++));
    std::vector<GLfloat> vss;
    for(const auto& f:def.getCons()._fixes) {
      //vertex
      for(int i=0; i<N; i++)
        vss.push_back((GLfloat)def.x()[f.first*N+i]);
      for(int i=N; i<3; i++)
        vss.push_back(0);
      //hook point
      for(int i=0; i<N; i++)
        vss.push_back((GLfloat)f.second.second[i]);
      for(int i=N; i<3; i++)
        vss.push_back(0);
    }
    m->setVertices(vss);
  }
  //l1
  if(!def.getL1ss().empty()) {
    m=std::dynamic_pointer_cast<DRAWER::MeshShape>(ret->getChild(id++));
    setVertices(def,m);
  }
}
}

#endif
