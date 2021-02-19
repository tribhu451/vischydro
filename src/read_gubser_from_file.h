#pragma once

#include<fstream>
#include<string>
#include<sstream>
#include "global.h"
#include<iostream>
#include<fstream>
#include <TRandom3.h>
#include <TF1.h>
#include "TMath.h"
#include "fluid.h"
#include "eos.h"
#include "input_data.h"
	
using namespace std;

class read_gubser_from_file
{
  
 public:
  read_gubser_from_file(InputData *InData1);
  ~read_gubser_from_file();
  void set_ic(fluid* f, EoS* eos, double tau_0);
  
  
 private:
  InputData *InData;
  istringstream* iss;
  char   buff[200];

};
















