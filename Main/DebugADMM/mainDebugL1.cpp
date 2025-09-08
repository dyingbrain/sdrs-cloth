#include <ADMM/SmoothL1.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  SmoothL1<2>().debugEnergy(10);
  SmoothL1<3>().debugEnergy(10);
  return 0;
}
