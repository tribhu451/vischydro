#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <chrono>
#include <cmath>
#include "idb.h"
#include "read_id.h"
#include "master.h"

using std::cout;
using std::endl;
using namespace std::chrono;


int main(int argc, char **argv)
{
  
  // reading the input data from the file
  string input_file_name;char* event_no_s ;

  if(argc == 3){
      event_no_s = *(argv+1); input_file_name = *(argv+2);}
  else
     {
      cout<<"[Info] Please give 2 arguments\n"
            "1st argument - event no\n    "
            "2nd argument - input filename"<<endl;
      exit(1);
     }

  read_id* reader = new read_id(); 
  idb *IDB = new idb;
  reader->read_id_from_file(IDB, input_file_name); // read input data base and store it.
  int event_no = atof(event_no_s) ;

  cout<<"event no : "<<event_no<<endl;


  master* head = new master(IDB);
  head->init(); // initialize
  head->run_hydro(); 
 

 delete head;
 return 0;
}



