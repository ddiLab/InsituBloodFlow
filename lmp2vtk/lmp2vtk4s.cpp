# include <iostream>
# include <fstream>
# include <sstream>
# include <unistd.h>
# include <string>
# include <stdlib.h>
# include <vector>
# include "lmp2vtk4s.h"

using namespace std;


int main(int argc, char **argv)
{
  int opt;
  int ncell(0),pxyz[3]={0,0,0};
  double lxyz[6]={0,0,0,0,0,0};
  string fp_name(""),ftp_name(""),fw_name(""),ff_name(""),fb_name(""), fv_name(""), fa_name("");
  string style("");
  while ((opt =getopt(argc,argv,"p:t:b:f:s:w:c:x:y:z:v:a:")) != -1){
    switch(opt){
      case 'a':
        fa_name =string(optarg);
        break;
      case 'p'://position file xyz
        fp_name =string(optarg);
        cout<<"file name "<<fp_name<<endl;
        break;
      case 't'://topology, in.cells
        ftp_name = string(optarg);
        break;
      case 'b':
        fb_name = string(optarg);
        break;
      case 'f'://force file
        ff_name = string(optarg);
        break;
      case 'v'://velocity
        fv_name = string(optarg);
        break;
      case 's'://style: cell or chain
        style = string(optarg);
        break;
      case 'w'://output directory and file name prefix
        fw_name = string(optarg);
        break;
      case 'c': //number of cells
        ncell = atoi(optarg);
        cout<<"ncell "<<ncell<<endl;
        break;
      case 'x'://periodic in x
        pxyz[0] = atoi(optarg);
        break;
      case 'y':
        pxyz[1] = atoi(optarg);
        break;
      case 'z':
        pxyz[2] = atoi(optarg);
        break;
      default:
        //cerr<<"error"<<endl;
        cout<<"usage: a.out -p dump.xyz -t in.cells -w vtk_solid/cell -c 31 -x 1 "<<endl;
        return 0;
    }
  }
  
  //readWrite(ftp_name,fp_name, ff_name, fv_name, fw_name, style, ncell, pxyz);
  readWrite(ftp_name,fp_name, ff_name, fa_name, fv_name, fw_name, style, ncell, pxyz);
  return 0;
}


