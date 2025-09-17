#include <ADMM/ARAP.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  mpfr_set_default_prec(1024);
  ARAP<2>().debugEnergy(10,1e-4f);
  ARAP<3>().debugEnergy(10,1e-4f);
  return 0;
}
