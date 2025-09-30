#include "SkirtDatasetCommon.h"

int main(int argc,char** argv) {
  //moving directions
  std::vector<typename Deformable<3>::VecNT> dirs;
  dirs.push_back(Deformable<3>::VecNT({1,0,0}).normalized());
  runDataset(argc,argv,true,dirs);
  return 0;
}
