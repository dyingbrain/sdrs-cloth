#include <ADMM/MassSpring.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  MassSpring<2>().debugEnergy(10,1e-4f);
  MassSpring<3>().debugEnergy(10,1e-4f);
  return 0;
}
