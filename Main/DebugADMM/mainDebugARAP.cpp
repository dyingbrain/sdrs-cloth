#include <ADMM/ARAP.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  ARAP<2>().debugEnergy(10,1e-4f);
  ARAP<3>().debugEnergy(10,1e-4f);
  return 0;
}
