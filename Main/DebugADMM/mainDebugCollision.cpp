#include <ADMM/CollisionSelf.h>
#include <ADMM/CollisionObstacle.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  mpfr_set_default_prec(1024);
  CollisionMatrix<2,4>::debug();
  CollisionMatrix<3,6>::debug();
  CollisionSelf<2,2>(1e-3f,5).debugEnergy(10,1e-4f);
  CollisionSelf<3,3>(1e-3f,5).debugEnergy(10,1e-4f);
  CollisionObstacle<2,2,4>(1e-3f,5).debugEnergy(10,1e-4f);
  CollisionObstacle<3,3,6>(1e-3f,5).debugEnergy(10,1e-4f);
  return 0;
}
