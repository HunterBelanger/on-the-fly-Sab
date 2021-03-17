#include <ENDFtk.hpp>
#include <Panglos.hpp>

using namespace njoy::ENDFtk;

#include <iostream>
#include <string>

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cout << " ERROR: Must provide 2 arguments.\n";
    return 1;
  }

  int nin = std::stoi(argv[1]);
  std::string fname = "tape" + std::to_string(nin);

  int mat = std::stoi(argv[2]);

  std::cout << " Read tape: " << fname << std::endl;
  std::cout << " MAT: " << mat << std::endl;

  tree::Tape<std::string> pendf = tree::fromFile(fname);

  file::Type<7> MF7 = pendf.material(mat).front().file(7).parse<7>();

  Panglos panglos(MF7);

  std::string input = "";
  std::cout << "Press any key to continue...\n";
  std::cin >> input;

  return 0;
}
