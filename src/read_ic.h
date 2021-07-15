#pragma once

#include<fstream>
#include<string>
#include<sstream>
#include<iostream>
#include<fstream>
#include "idb.h"
#include "grid.h"

	
using namespace std;

class read_ic
{
  
 public:
  read_ic(idb* );
  ~read_ic();
  void set_ic(grid* , EoS* );
  
  
 private:
  idb *IDB;

inline int theta(double _a)
{
  if(_a > 0){return 1;}else{return 0;}
}

  istringstream* iss;
  char   buff[200];

};
















