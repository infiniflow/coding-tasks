#include <iostream>
#include <fstream>
#include <vector>
#include "wavelet_tree.hpp"

using namespace std;

int ReadArrayFromFile(const std::string& file_name, vector<uint64_t>& array){
  array.clear();
  ifstream ifs(file_name.c_str());
  if (!ifs){
    cerr << "Unable to open [" << file_name << "]" << endl;
    return -1;
  }

  for (uint64_t val; ifs >> val; ){
    array.push_back(val);
  }
  return 0;
}

int BuildIndex(const string& input_file_name,
	       const string& index_name){
  vector<uint64_t> array;
  if (ReadArrayFromFile(input_file_name, array) == -1){
    return -1;
  }

  wavelet::WaveletTree wa;
  wa.Init(array);
  
  ofstream ofs(index_name.c_str());
  if (!ofs){
    cerr << "Unable to open [" << index_name << "]" << endl;
    return -1;
  }
  wa.Save(ofs);

  if (!ofs){
    cerr << "Save failed [" << index_name << "]" << endl;
    return -1;
  }

  return 0;
}


int main(int argc, char* argv[]){
  
  return 0;
}
