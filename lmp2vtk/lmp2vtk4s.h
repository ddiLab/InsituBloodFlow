# include <iostream>
# include <fstream>
# include <sstream>
# include <unistd.h>
# include <string>
# include <stdlib.h>
# include <vector>
//# include <direct>

using namespace std;

void readXYZ(string fp_name, vector<int> &time, vector< vector<double> > &pos);
void readForce(string fp_name, int nn, vector< vector<double> > &force);
void readVelocity(string fp_name, int nn, vector< vector<double> > &vel);

/*changing atom numbers every time step*/
void readXYZ(string fp_name, vector<int> &time, vector<int> &n_atoms, vector< vector<double> > &pos);
//void readForce(string fp_name, vector<int> &n_atoms, vector< vector<double> > &force);
/*eof changing atom numbers*/
void readTopo(string fp_name, vector< vector<int> > &face, vector< vector<int> > &bond, double* lxyz);
void readAdBond(string fp_name, int nb, vector< vector<int> > &bond);
void writeVTKsolid(string fw_name, string style,vector<int> &time, vector< vector<double> > &pos, vector< vector<int> > &face, vector< vector<int> > &bond);
/*write force vector*/
void writeVTKsolid(string fw_name, string style,vector<int> &time,vector<int> &n_atoms, vector< vector<double> > &pos, vector< vector<double> > &force,vector< vector<int> > &face, vector< vector<int> > &bond);
/*write force, velocity vector*/
void writeVTKsolid(string fw_name, string style,vector<int> &time,vector<int> &n_atoms, vector< vector<double> > &pos, vector< vector<double> > &force,vector< vector<double> > &vel,vector< vector<int> > &face, vector< vector<int> > &bond);
void writeVTKsolid(string fw_name, string style,vector<int> &time, vector<int> &n_atoms,vector< vector<double> > &pos, vector< vector<double> > &force,vector< vector<int> > &face, vector< vector<int> > &bond, vector< vector<int> > &adBond);
void periodicCorrection(vector<int> &time, vector< vector<double> > &pos, int nc, int pxyz[3], double lxyz[6]); 
void periodicCorrection(vector<int> &time, vector< vector<double> > &pos, int nc, int nOneCell, int pxyz[3], double lxyz[6]); 
void writeXYZ(string fw_name, vector<int> &time, vector< vector<double> > &pos);

void readWrite(string ftp_name, string fp_name, string ff_name, string fv_name, string fw_name, string style, int ncell, int pxyz[3] );

void readWrite(string ftp_name, string fp_name, string ff_name, string fa_name, string fv_name, string fw_name, string style, int ncell, int pxyz[3] );
//-------- Read and Write single time step, Less memory comsumption---------------------//
void readWrite(string ftp_name, string fp_name, string ff_name, string fa_name, string fv_name, string fw_name, string style, int ncell, int pxyz[3] ){
  
  ifstream fp;
  ifstream ff;
  vector <double> row(3,0.0);
  vector <double> row_f(4,0.0);
  vector <double> row_strain(2,0.0);
  if (fp_name.empty()) {
    cout<<"Error, please specify position file: -p filename"<<endl; 
    return;
  }
  fp.open(fp_name.c_str());
 
  vector < vector<int> > face;
  vector < vector<int> > bond;
  vector < vector<double> > force;
  vector < vector<double> > vel;
  int flag_f=0;
  int flag_v=0;
  int flag_a=0;
  if (!ff_name.empty()) flag_f=1;
  if (!fv_name.empty()) flag_v=1;
  if (!fa_name.empty()) flag_a=1;
  /*if (flag_f){
    ff.open(ff_name.c_str());
  }*/
  if (flag_a){
    ff.open(ff_name.c_str());
  }
  double lxyz[6]={0,0,0,0,0,0};
  if (ftp_name.empty()) {
    cout<<"Error, please specify topology file: -t filename"<<endl; 
    return;
  }
  readTopo(ftp_name,face,bond,lxyz);
  
  int nb = bond.size();
  int na = face.size();
  
  int n,i,j;
  int tmp,timestep;
  string buf;
  vector <int> time;
  vector < vector<double> > pos;
  if (fw_name.empty()) {
    cout<<"Error, please specify output folder and file name, -w dir/filename"<<endl; 
    return;
  }
  stringstream output_filename;
  ofstream of;
  // read position each time step  
  while(fp>>n){

    fp>>buf>>buf>>buf;
    timestep = atoi(buf.c_str());
    time.push_back(timestep);
    
    for (i=0;i<n;i++){
      fp>>tmp>>row[0]>>row[1]>>row[2];
      pos.push_back(row);
    }
    //periodicCorrection(time, pos, ncell, 10242, pxyz, lxyz);
    periodicCorrection(time, pos, ncell,  pxyz, lxyz);

    if (flag_a){
      for (i=0;i<7;i++) ff>>buf;
      ff>>buf;// number of atoms
      for (i=0;i<16;i++) ff>>buf;
      for (i=0;i<n;i++){ 
        //ff>>row_f[0]>>row_f[1]>>row_f[2]>>row_f[3]; // for force
        //force.push_back(row_f);
        ff>>row_strain[0]>>row_strain[1]; // for hydrostatic strain
        force.push_back(row_strain);
      }
    }
    //write coordiantes
    output_filename << fw_name << timestep << ".vtk";
   // output_filename = fw_name + timestep + ".vtk";
    //cout<<output_filename.str()<<endl; 
    of.open(output_filename.str().c_str());
    //of.open(output_filename);
    of << "# vtk DataFile Version 3.0\n";
    of << "Cells\n";
    of << "ASCII\n";
    of << "DATASET POLYDATA\n";

    //n=30726;//for current cells only
    //int n1=5136;//for current cells only
    //int n1=204840;//for current cells only
    of << "POINTS "<<n<< " float\n";
    //of << "POINTS "<<n-n1<< " float\n";
        for (j=0;j<n;j++){
       // for (j=n1;j<n;j++)
          of << pos[j][0]<<" "<< pos[j][1]<<" "<<pos[j][2]<<endl;
        }
        
    if (style.compare("chain")==0){
      of << "LINES "<<nb<<" "<<3*nb<<endl;
      for (j=0;j<nb;j++)
        of <<"2 "<<bond[j][0]-1<<" "<<bond[j][1]-1<<endl;
    }else if(style.compare("cell")==0){
      of << "POLYGONS "<<na<<" "<<4*na<<endl;
      for (j=0;j<na;j++)
        of <<"3 "<<face[j][0]-1<<" "<<face[j][1]-1<<" "<<face[j][2]-1<<endl;
    }
    /*
    // force vector 
    if (flag_f){
      of << "POINT_DATA "<<n<< " \n";
      //of << "POINT_DATA "<<n-n1<< " \n";
      of << "VECTORS force float\n";
      for (j=0;j<n;j++)
      //for (j=n1;j<n;j++)
        of << force[j][1]<<" "<< force[j][2]<<" "<<force[j][3]<<endl;
    }*/
    // scalar for hydrostatic strain 
    if (flag_a){
      of << "POINT_DATA "<<n<< " \n";
      //of << "POINT_DATA "<<n-n1<< " \n";
      of << "SCALARS strain float 1\n";
      of << "LOOKUP_TABLE default\n";
      for (j=0;j<n;j++)
      //for (j=n1;j<n;j++)
        of << force[j][1]<<endl;
    }

    
    of.close();
    of.clear();
    if (flag_a) {
      force.clear();
    }
    //clear the position 
    output_filename.str(string());
    pos.clear();    
    time.clear();
    cout<<"time step "<<timestep<<endl;
  }
  if (flag_a) ff.close();
  fp.close();

}


//-------- Read and Write single time step, Less memory comsumption---------------------//
void readWrite(string ftp_name, string fp_name, string ff_name, string fv_name, string fw_name, string style, int ncell, int pxyz[3] ){
 /* ftp_name: topology file, [-t] e.g., in.cells. 
 *  fp_name: atom position file,[-p] e.g, dump.rbc.xyz
 *  ff_name: force file, each row is a force vector for each atom [-f]
 *  fv_name: velocity file, each row is a velocity vector for eeach atom [-v]
 *  fw_name: the prefix for output file, you can include folder name.[-w] e.g., vtk/cell means all the files will be generated with a name cell+timestep.vtk in folder vtk.
 *  style: two potions: cell or chain. Cell is for cell mesh, chain is for polymer. [-s] 
 *  ncell: number of cells [-c]
 *  pxyz: periodic flags for x, y, z axis so that we can stitch straddling cells in periodic boundary conditions [-x, -y, -z]. 
 *  usage: ./lmp2vtk4s -p dump.rbc.xyz -t in.cells -w vtk/cell -s cell -c 3 -z 1
 * */ 
  ifstream fp;
  ifstream ff;
  vector <double> row(3,0.0);
  vector <double> row_f(4,0.0);
  vector <double> row_strain(2,0.0);
  if (fp_name.empty()) {
    cout<<"Error, please specify position file: -p filename"<<endl; 
    return;
  }
  fp.open(fp_name.c_str());
 
  vector < vector<int> > face;
  vector < vector<int> > bond;
  vector < vector<double> > force;
  vector < vector<double> > vel;
  int flag_f=0;
  int flag_v=0;
  if (!ff_name.empty()) flag_f=1;
  if (!fv_name.empty()) flag_v=1;
  if (flag_f){
    ff.open(ff_name.c_str());
  }

  double lxyz[6]={0,0,0,0,0,0};
  if (ftp_name.empty()) {
    cout<<"Error, please specify topology file: -t filename"<<endl; 
    return;
  }
  readTopo(ftp_name,face,bond,lxyz); //read the topology of the file
  
  int nb = bond.size();
  int na = face.size();
  
  int n,i,j;
  int tmp,timestep;
  string buf;
  vector <int> time;
  vector < vector<double> > pos;
  if (fw_name.empty()) {
    cout<<"Error, please specify output folder and file name, -w dir/filename"<<endl; 
    return;
  }
  stringstream output_filename;
  ofstream of;
  // read position each time step  
  while(fp>>n){

    fp>>buf>>buf>>buf;
    timestep = atoi(buf.c_str());
    time.push_back(timestep);
    
    for (i=0;i<n;i++){
      fp>>tmp>>row[0]>>row[1]>>row[2];
      pos.push_back(row);
    }
    periodicCorrection(time, pos, ncell, 10242, pxyz, lxyz); // atoms_per_cell: 10242 is manually specified
    if (flag_f){
      for (i=0;i<7;i++) ff>>buf;
      ff>>buf;// number of atoms
      for (i=0;i<18;i++) ff>>buf;
      for (i=0;i<n;i++){ 
        //ff>>row_f[0]>>row_f[1]>>row_f[2]>>row_f[3]; // for force
        //force.push_back(row_f);
        ff>>row_strain[0]>>row_strain[1]; // for hydrostatic strain
        force.push_back(row_strain);
      }
    }
    //write coordiantes
    output_filename << fw_name << timestep << ".vtk";
   // output_filename = fw_name + timestep + ".vtk";
    //cout<<output_filename.str()<<endl; 
    of.open(output_filename.str().c_str());
    //of.open(output_filename);
    of << "# vtk DataFile Version 3.0\n";
    of << "Cells\n";
    of << "ASCII\n";
    of << "DATASET POLYDATA\n";
//    n=30726;//Jifu: this line may be useful for other combined cases where both cells and platelets are dumped together.  
    //int n1=5136;//for current cells only
    //int n1=204840;//for current cells only
    of << "POINTS "<<n<< " float\n";
    //of << "POINTS "<<n-n1<< " float\n";
        for (j=0;j<n;j++){
       // for (j=n1;j<n;j++)
          of << pos[j][0]<<" "<< pos[j][1]<<" "<<pos[j][2]<<endl;
        }
        
    if (style.compare("chain")==0){
      of << "LINES "<<nb<<" "<<3*nb<<endl;
      for (j=0;j<nb;j++)
        of <<"2 "<<bond[j][0]-1<<" "<<bond[j][1]-1<<endl;
    }else if(style.compare("cell")==0){
      of << "POLYGONS "<<na<<" "<<4*na<<endl;
      for (j=0;j<na;j++)
        of <<"3 "<<face[j][0]-1<<" "<<face[j][1]-1<<" "<<face[j][2]-1<<endl;
    }
    /*
    // force vector 
    if (flag_f){
      of << "POINT_DATA "<<n<< " \n";
      //of << "POINT_DATA "<<n-n1<< " \n";
      of << "VECTORS force float\n";
      for (j=0;j<n;j++)
      //for (j=n1;j<n;j++)
        of << force[j][1]<<" "<< force[j][2]<<" "<<force[j][3]<<endl;
    }*/
    // scalar for hydrostatic strain 
    if (flag_f){
      of << "POINT_DATA "<<n<< " \n";
      //of << "POINT_DATA "<<n-n1<< " \n";
      of << "SCALARS strain float 1\n";
      of << "LOOKUP_TABLE default\n";
      for (j=0;j<n;j++)
      //for (j=n1;j<n;j++)
        of << force[j][1]<<endl;
    }

    
    of.close();
    of.clear();
    if (flag_f) {
      force.clear();
    }
    //clear the position 
    output_filename.str(string());
    pos.clear();    
    time.clear();
    cout<<"time step "<<timestep<<endl;
  }
  if (flag_f) ff.close();
  fp.close();

}



void readXYZ(string fp_name, vector<int> &time, vector<int> &n_atoms, vector< vector<double> > &pos){
  ifstream fp;
  vector <double> row(3,0.0);
  fp.open(fp_name.c_str());
 
  int n;
  int tmp,timestep;
  string buf;
  while(fp>>n){

    fp>>buf>>buf>>buf;
    timestep = atoi(buf.c_str());
    time.push_back(timestep);
    n_atoms.push_back(n);

    for (int i=0;i<n;i++){
      fp>>tmp>>row[0]>>row[1]>>row[2];
      pos.push_back(row);
    }
  }
  fp.close();
}

void readXYZ(string fp_name, vector<int> &time, vector< vector<double> > &pos){
  ifstream fp;
  vector <double> row(3,0.0);
  fp.open(fp_name.c_str());
 
  int n;
  int tmp,timestep;
  string buf;
  while(fp>>n){

    fp>>buf>>buf>>buf;
    timestep = atoi(buf.c_str());
    time.push_back(timestep);
    
    for (int i=0;i<n;i++){
      fp>>tmp>>row[0]>>row[1]>>row[2];
      pos.push_back(row);
    }
  }
  fp.close();
}

void readForce(string fp_name,int nn, vector< vector<double> > &force){
  ifstream fp;
  vector <double> row(4,0.0);
  fp.open(fp_name.c_str());
  vector <int> n_atoms_local;

  int i,n,id;
  int tmp,ts(0);
  long long sum(0);
  string buf;
  if (fp.is_open()){
    while(!fp.eof()){
      for (i=0;i<7;i++)
        fp>>buf;
      fp>>n;
      //cout<<"n "<<n<<" ts "<<ts<<endl;
      n_atoms_local.push_back(n);
      //sum += n;
      for (i=0;i<18;i++)
        fp>>buf;
      for (i=0;i<n;i++){
        //fp>>id;
        //tmp = sum +id-1;
        //fp>>force[tmp][0]>>force[tmp][1]>>force[tmp][2];
        fp>>row[0]>>row[1]>>row[2]>>row[3];
        force.push_back(row);
      }
      ts++;
    }
  }
  
  fp.close();
/*---debug------------*/
  /*
  ofstream of;
  of.open("forceOut",std::ofstream::app);
  sum =0;
  for (i=0;i<ts-1;i++){
    for (int j=0;j<n_atoms_local[i];j++){
      of << force[sum+j][0]<<" "<< force[sum+j][1]<<" "<<force[sum+j][2]<<endl;
    }
    sum += n_atoms_local[i];
  }
  of.close();
*/
}

void readVelocity(string fp_name,int nn, vector< vector<double> > &vel){
  ifstream fp;
  vector <double> row(4,0.0);
  fp.open(fp_name.c_str());
  vector <int> n_atoms_local;

  int i,n,id;
  int tmp,ts(0);
  long long sum(0);
  string buf;
  if (fp.is_open()){
    while(!fp.eof()){
      for (i=0;i<7;i++)
        fp>>buf;
      fp>>n;
      //cout<<"n "<<n<<" ts "<<ts<<endl;
      n_atoms_local.push_back(n);
      //sum += n;
      for (i=0;i<18;i++)
        fp>>buf;
      for (i=0;i<n;i++){
        //fp>>id;
        //tmp = sum +id-1;
        //fp>>force[tmp][0]>>force[tmp][1]>>force[tmp][2];
        fp>>row[0]>>row[1]>>row[2]>>row[3];
        vel.push_back(row);
      }
      ts++;
    }
  }
  
  fp.close();
/*---debug------------*/
  /*
  ofstream of;
  of.open("forceOut",std::ofstream::app);
  sum =0;
  for (i=0;i<ts-1;i++){
    for (int j=0;j<n_atoms_local[i];j++){
      of << force[sum+j][0]<<" "<< force[sum+j][1]<<" "<<force[sum+j][2]<<endl;
    }
    sum += n_atoms_local[i];
  }
  of.close();
*/
}
void readAdBond(string fp_name,int nb, vector< vector<int> > &bond){
  ifstream fp;
  fp.open(fp_name.c_str());
 
  int i,n,id,i1,i2;
  int tmp,ts(0);
  string buf;
  if (fp.is_open()){
    while(!fp.eof()){
      vector <int> row;//has to be here, wrong if out of this loop
      for (i=0;i<7;i++) fp>>buf;
      fp>>n;//number of entries
      //cout<<"Entries "<<n<<endl;
      for (i=0;i<19;i++) fp>>buf;
      for (i=0;i<n;i++){
        fp>>buf>>id>>i1>>i2>>buf;
        if (id > nb) {
          row.push_back(i1);
          row.push_back(i2);
          row.push_back(id);
        }
      }
      bond.push_back(row);
      ts++;
    }
  }
  fp.close();
 /* 
  ofstream of;
  of.open("adBond",std::ofstream::app);
  for (i=0;i<ts-1;i++){
    n = bond[i].size();
    cout<<"n "<<n<<endl;
    for (int j=0;j<n/2;j++)
      of << bond[i][2*j]<<" "<< bond[i][2*j+1]<<" "<<endl;
  }
  of.close();*/
}

void readTopo(string ftp_name, vector< vector<int> > &face, vector< vector<int> > &bond, double lxyz[6]){
  fstream ftp;
  ftp.open(ftp_name.c_str());
  vector<int> row(3,0); // 3 node for a face
  vector<int> bd0(2,0); // 2 node for a bond
  int na,count,tmp;
  int natom, nbond;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  count=0;
  string buf,s_str;
  size_t loc;
  ftp>>buf;//first line

  while (getline(ftp,buf)){
    loc = buf.find("atoms");
    if (loc != string::npos)  {
      natom = atoi(buf.substr(0,loc).c_str());
      cout<<natom<<" atoms"<<endl;
    }
    loc = buf.find("bonds");
    if (loc != string::npos)  {
      nbond = atoi(buf.substr(0,loc).c_str());
      cout<<nbond<<" bonds"<<endl;
    }
    loc = buf.find("angles");
    if (loc != string::npos)  {
      na = atoi(buf.substr(0,loc).c_str());
      cout<<na<<" angles"<<endl;
    }
    loc = buf.find("xlo");
    if (loc != string::npos){
      s_str=buf.substr(0,loc);
      loc = s_str.find(" ");
      xlo = atof(s_str.substr(0,loc).c_str());
      xhi = atof(s_str.substr(loc+1,s_str.length()).c_str());
      lxyz[0]=xlo;
      lxyz[1]=xhi;
      cout<<xlo<<" "<<xhi<<" xlo xhi"<<endl;
    }
    loc = buf.find("ylo");
    if (loc != string::npos){
      s_str=buf.substr(0,loc);
      loc = s_str.find(" ");
      ylo = atof(s_str.substr(0,loc).c_str());
      yhi = atof(s_str.substr(loc+1,s_str.length()).c_str());
      lxyz[2]=ylo;
      lxyz[3]=yhi;
      cout<<ylo<<" "<<yhi<<" ylo yhi"<<endl;
    }
    loc = buf.find("zlo");
    if (loc != string::npos){
      s_str=buf.substr(0,loc);
      loc = s_str.find(" ");
      zlo = atof(s_str.substr(0,loc).c_str());
      zhi = atof(s_str.substr(loc+1,s_str.length()).c_str());
      lxyz[4]=zlo;
      lxyz[5]=zhi;
      cout<<zlo<<" "<<zhi<<" zlo zhi"<<endl;
    }
    loc = buf.find("Angles");
    if (loc != string::npos)  {
      for (int i=0;i<na;i++){
        ftp>>tmp>>tmp>>row[0]>>row[1]>>row[2];
        face.push_back(row);
      }
    }
    loc = buf.find("Bonds");
    if (loc != string::npos)  {
      for (int i=0;i<nbond;i++){
        ftp>>tmp>>tmp>>bd0[0]>>bd0[1];
        bond.push_back(bd0);
      }
    }
  }
  ftp.close();
}

void writeVTKsolid(string fw_name, string style, vector<int> &time, vector< vector<double> > &pos, vector< vector<int> > &face, vector< vector<int> > &bond ){
      int i,j,k;
      int nn = pos.size()/time.size();
      int nb = bond.size();
      int na = face.size();
      //cout<<"face size "<<face.size()<<"ts "<<time.size()<<endl;
      for (i=0;i<time.size();i++){
        stringstream output_filename;
        output_filename << fw_name << time[i] << ".vtk";
    
        ofstream of;
        of.open(output_filename.str().c_str());
        of << "# vtk DataFile Version 3.0\n";
        of << "Cells\n";
        of << "ASCII\n";
        of << "DATASET POLYDATA\n";
        
        of << "POINTS "<<nn<< " float\n";
        for (j=0;j<nn;j++)
          of << pos[i*nn+j][0]<<" "<< pos[i*nn+j][1]<<" "<<pos[i*nn+j][2]<<endl;
        
        if (style.compare("chain")==0){
          of << "LINES "<<nb<<" "<<3*nb<<endl;
          for (j=0;j<nb;j++)
            of <<"2 "<<bond[j][0]-1<<" "<<bond[j][1]-1<<endl;
        }else if(style.compare("cell")==0){
          of << "POLYGONS "<<na<<" "<<4*na<<endl;
          for (j=0;j<na;j++)
            of <<"3 "<<face[j][0]-1<<" "<<face[j][1]-1<<" "<<face[j][2]-1<<endl;
        }
        of.close();
      }
}

/*--with force vector--*/
void writeVTKsolid(string fw_name, string style,vector<int> &time,vector<int> &n_atoms, vector< vector<double> > &pos, vector< vector<double> > &force,vector< vector<int> > &face, vector< vector<int> > &bond){
  int i,j,k;
  //int nn = pos.size()/time.size();
  int nn;
  int nb = bond.size();
  int na = face.size();
  long long sum(0);
  //cout<<"face size "<<face.size()<<"ts "<<time.size()<<endl;
  for (i=0;i<time.size();i++){
    stringstream output_filename;
    output_filename << fw_name << time[i] << ".vtk";
  
    ofstream of;
    of.open(output_filename.str().c_str());
    of << "# vtk DataFile Version 3.0\n";
    of << "Cells\n";
    of << "ASCII\n";
    of << "DATASET POLYDATA\n";
    
    nn = n_atoms[i];

    of << "POINTS "<<nn<< " float\n";
    for (j=0;j<nn;j++)
      of << pos[sum+j][0]<<" "<< pos[sum+j][1]<<" "<<pos[sum+j][2]<<endl;
    
    if (style.compare("chain")==0){
      of << "LINES "<<nb<<" "<<3*nb<<endl;
      for (j=0;j<nb;j++)
        of <<"2 "<<bond[j][0]-1<<" "<<bond[j][1]-1<<endl;
    }else if(style.compare("cell")==0){
      of << "POLYGONS "<<na<<" "<<4*na<<endl;
      for (j=0;j<na;j++)
        of <<"3 "<<face[j][0]-1<<" "<<face[j][1]-1<<" "<<face[j][2]-1<<endl;
    }

    of << "POINT_DATA "<<nn<< " \n";
    of << "VECTORS force float\n";
    for (j=0;j<nn;j++)
      of << force[sum+j][1]<<" "<< force[sum+j][2]<<" "<<force[sum+j][3]<<endl;
    
    of.close();
    sum += nn;
  }
}

void writeVTKsolid(string fw_name, string style,vector<int> &time,vector<int> &n_atoms, vector< vector<double> > &pos, vector< vector<double> > &force,vector< vector<double> > &vel,vector< vector<int> > &face, vector< vector<int> > &bond){
  int i,j,k;
  //int nn = pos.size()/time.size();
  int nn;
  int nb = bond.size();
  int na = face.size();
  long long sum(0);
  //cout<<"face size "<<face.size()<<"ts "<<time.size()<<endl;
  for (i=0;i<time.size();i++){
    stringstream output_filename;
    output_filename << fw_name << time[i] << ".vtk";
  
    ofstream of;
    of.open(output_filename.str().c_str());
    of << "# vtk DataFile Version 3.0\n";
    of << "Cells\n";
    of << "ASCII\n";
    of << "DATASET POLYDATA\n";
    
    nn = n_atoms[i];

    of << "POINTS "<<nn<< " float\n";
    for (j=0;j<nn;j++)
      of << pos[sum+j][0]<<" "<< pos[sum+j][1]<<" "<<pos[sum+j][2]<<endl;
    
    if (style.compare("chain")==0){
      of << "LINES "<<nb<<" "<<3*nb<<endl;
      for (j=0;j<nb;j++)
        of <<"2 "<<bond[j][0]-1<<" "<<bond[j][1]-1<<endl;
    }else if(style.compare("cell")==0){
      of << "POLYGONS "<<na<<" "<<4*na<<endl;
      for (j=0;j<na;j++)
        of <<"3 "<<face[j][0]-1<<" "<<face[j][1]-1<<" "<<face[j][2]-1<<endl;
    }

    of << "POINT_DATA "<<nn<< " \n";
    of << "VECTORS force float\n";
    for (j=0;j<nn;j++)
      of << force[sum+j][1]<<" "<< force[sum+j][2]<<" "<<force[sum+j][3]<<endl;
    of << "VECTORS velocity float\n";
    for (j=0;j<nn;j++)
      of << vel[sum+j][1]<<" "<< vel[sum+j][2]<<" "<<vel[sum+j][3]<<endl;
    

    of.close();
    sum += nn;
  }
}
void writeVTKsolid(string fw_name, string style,vector<int> &time, vector<int> &n_atoms, vector< vector<double> > &pos, vector< vector<double> > &force,vector< vector<int> > &face, vector< vector<int> > &bond, vector< vector<int> > &adBond){
  int i,j,k;
  //int nn = pos.size()/time.size();
  int nn;
  int nb = bond.size();
  int na = face.size();
  long long sum(0);
  //cout<<"face size "<<face.size()<<"ts "<<time.size()<<endl;
  /*
  for (i=0;i<time.size();i++){
    stringstream output_filename;
    output_filename << fw_name << time[i] << ".vtk";
  
    ofstream of;
    of.open(output_filename.str().c_str());
    of << "# vtk DataFile Version 3.0\n";
    of << "Cells\n";
    of << "ASCII\n";
    of << "DATASET POLYDATA\n";
    
    of << "POINTS "<<nn<< " float\n";
    for (j=0;j<nn;j++)
      of << pos[i*nn+j][0]<<" "<< pos[i*nn+j][1]<<" "<<pos[i*nn+j][2]<<endl;
    
    if (style.compare("chain")==0){
      of << "LINES "<<nb<<" "<<3*nb<<endl;
      for (j=0;j<nb;j++)
        of <<"2 "<<bond[j][0]-1<<" "<<bond[j][1]-1<<endl;
    }else if(style.compare("cell")==0){
      nb = adBond[i].size()/2;//notice factor 1/2
      of << "LINES "<<nb<<" "<<3*nb<<endl;
      for (j=0;j<nb;j++)
        of <<"2 "<<adBond[i][2*j]-1<<" "<<adBond[i][2*j+1]-1<<endl;
      of << "POLYGONS "<<na<<" "<<4*na<<endl;
      for (j=0;j<na;j++)
        of <<"3 "<<face[j][0]-1<<" "<<face[j][1]-1<<" "<<face[j][2]-1<<endl;
    }

    of << "POINT_DATA "<<nn<< " \n";
    of << "VECTORS force float\n";
    for (j=0;j<nn;j++)
      of << force[i*nn+j][0]<<" "<< force[i*nn+j][1]<<" "<<force[i*nn+j][2]<<endl;
    
    of.close();
  }*/
  
  for (i=0;i<time.size();i++){
    stringstream output_filename;
    output_filename << fw_name <<"bond"<< time[i] << ".vtk";
  
    ofstream of;
    of.open(output_filename.str().c_str());
    of << "# vtk DataFile Version 3.0\n";
    of << "Cells\n";
    of << "ASCII\n";
    of << "DATASET UNSTRUCTURED_GRID\n";
    
    nn = n_atoms[i];

    of << "POINTS "<<nn<< " float\n";
    for (j=0;j<nn;j++)
      of << pos[sum+j][0]<<" "<< pos[sum+j][1]<<" "<<pos[sum+j][2]<<endl;

    nb = adBond[i].size()/3;//factor 3
    of << "CELLS "<<nb<<" "<<3*nb<<endl; 
      for (j=0;j<nb;j++)
        of <<"2 "<<adBond[i][3*j]-1<<" "<<adBond[i][3*j+1]-1<<endl;
   
    of <<"CELL_TYPES "<<nb<<endl;
      for (j=0;j<nb;j++)
        of <<"3 "<<endl;

    of << "POINT_DATA "<<nn<< " \n";
    of << "VECTORS force float\n";
    for (j=0;j<nn;j++)
      of << force[sum+j][0]<<" "<< force[sum+j][1]<<" "<<force[sum+j][2]<<endl;
    
    of << "CELL_DATA "<<nb<< " \n";
    of << "SCALARS bondType int 1\n";
    of << "LOOKUP_TABLE default\n";
    for (j=0;j<nb;j++)
      of << adBond[i][3*j+2]<<endl;

    sum += nn;
    
    of.close();
  }

}


void periodicCorrection(vector<int> &time, vector< vector<double> > &pos, int nc, int pxyz[3], double lxyz[6]){
  int i,j,k;
  int nn = pos.size()/time.size();
  int nOneCell = nn/nc;
  //cout<<"node per cell "<<nOneCell<<endl;
  double xmin,ymin,zmin,xmax,ymax,zmax;
  double span;
  for(k=0;k<time.size();k++){//time
    vector <int> flag(nc,0);
    int st = k*nn;
    // check if crossing
    for (i=0;i<nc;i++){
      xmin = xmax = pos[st+i*nOneCell][0];
      ymin = ymax = pos[st+i*nOneCell][1];
      zmin = zmax = pos[st+i*nOneCell][2];
      for (j=1;j<nOneCell;j++){
        if (xmin > pos[st+i*nOneCell+j][0]) xmin = pos[st+i*nOneCell+j][0];    
        if (xmax < pos[st+i*nOneCell+j][0]) xmax = pos[st+i*nOneCell+j][0];    
        if (ymin > pos[st+i*nOneCell+j][1]) ymin = pos[st+i*nOneCell+j][1];    
        if (ymax < pos[st+i*nOneCell+j][1]) ymax = pos[st+i*nOneCell+j][1];    
        if (zmin > pos[st+i*nOneCell+j][2]) zmin = pos[st+i*nOneCell+j][2];    
        if (zmax < pos[st+i*nOneCell+j][2]) zmax = pos[st+i*nOneCell+j][2];    
      }
      if ((xmax-xmin) > 0.5*(lxyz[1]-lxyz[0]) && pxyz[0] ) flag[i]=1;
      if ((ymax-ymin) > 0.5*(lxyz[3]-lxyz[2]) && pxyz[1] ) flag[i]=1;
      if ((zmax-zmin) > 0.5*(lxyz[5]-lxyz[4]) && pxyz[2]) flag[i]=1;
    }
    if (pxyz[0]){//x
      span=lxyz[1]-lxyz[0];
      for(i=0;i<nc;i++){
        if (flag[i]){
            for (j=0;j<nOneCell;j++)//0.5*lxyz assume box starts from 0
              if (pos[st+i*nOneCell+j][0]<0.5*span + lxyz[0]) pos[st+i*nOneCell+j][0] += span;
        }
      } 
    }
    if (pxyz[1]){//y
      span=lxyz[3]-lxyz[2];
      for(i=0;i<nc;i++){
        if (flag[i]){
            for (j=0;j<nOneCell;j++)
              if (pos[st+i*nOneCell+j][1]<0.5*span+lxyz[2]) pos[st+i*nOneCell+j][1] += span;
        }
      } 
    }
    if (pxyz[2]){//z
      span=lxyz[5]-lxyz[4];
      for(i=0;i<nc;i++){
        if (flag[i]){
            for (j=0;j<nOneCell;j++)
              if (pos[st+i*nOneCell+j][2]<0.5*span+lxyz[4]) pos[st+i*nOneCell+j][2] += span;
        }
      } 
    }
  }
} 

void periodicCorrection(vector<int> &time, vector< vector<double> > &pos, int nc, int nOneCell, int pxyz[3], double lxyz[6]){
  int i,j,k;
  int nn = pos.size()/time.size();
  //int nOneCell = nn/nc;
  //cout<<"node per cell "<<nOneCell<<endl;
  double xmin,ymin,zmin,xmax,ymax,zmax;
  double span;
  for(k=0;k<time.size();k++){//time
    vector <int> flag(nc,0);
    int st = k*nn;
    // check if crossing
    for (i=0;i<nc;i++){
      xmin = xmax = pos[st+i*nOneCell][0];
      ymin = ymax = pos[st+i*nOneCell][1];
      zmin = zmax = pos[st+i*nOneCell][2];
      for (j=1;j<nOneCell;j++){
        if (xmin > pos[st+i*nOneCell+j][0]) xmin = pos[st+i*nOneCell+j][0];    
        if (xmax < pos[st+i*nOneCell+j][0]) xmax = pos[st+i*nOneCell+j][0];    
        if (ymin > pos[st+i*nOneCell+j][1]) ymin = pos[st+i*nOneCell+j][1];    
        if (ymax < pos[st+i*nOneCell+j][1]) ymax = pos[st+i*nOneCell+j][1];    
        if (zmin > pos[st+i*nOneCell+j][2]) zmin = pos[st+i*nOneCell+j][2];    
        if (zmax < pos[st+i*nOneCell+j][2]) zmax = pos[st+i*nOneCell+j][2];    
      }
      if ((xmax-xmin) > 0.5*(lxyz[1]-lxyz[0]) && pxyz[0] ) flag[i]=1;
      if ((ymax-ymin) > 0.5*(lxyz[3]-lxyz[2]) && pxyz[1] ) flag[i]=1;
      if ((zmax-zmin) > 0.5*(lxyz[5]-lxyz[4]) && pxyz[2]) flag[i]=1;
    }
    if (pxyz[0]){//x
      span=lxyz[1]-lxyz[0];
      for(i=0;i<nc;i++){
        if (flag[i]){
            for (j=0;j<nOneCell;j++)//0.5*lxyz assume box starts from 0
              if (pos[st+i*nOneCell+j][0]<0.5*span + lxyz[0]) pos[st+i*nOneCell+j][0] += span;
        }
      } 
    }
    if (pxyz[1]){//y
      span=lxyz[3]-lxyz[2];
      for(i=0;i<nc;i++){
        if (flag[i]){
            for (j=0;j<nOneCell;j++)
              if (pos[st+i*nOneCell+j][1]<0.5*span+lxyz[2]) pos[st+i*nOneCell+j][1] += span;
        }
      } 
    }
    if (pxyz[2]){//z
      span=lxyz[5]-lxyz[4];
      for(i=0;i<nc;i++){
        if (flag[i]){
            for (j=0;j<nOneCell;j++)
              if (pos[st+i*nOneCell+j][2]<0.5*span+lxyz[4]) pos[st+i*nOneCell+j][2] += span;
        }
      } 
    }
  }
}


/*write pure xyz file for matlab processing*/
void writeXYZ(string fw_name, vector<int> &time, vector< vector<double> > &pos){
  int i,j,k;
  int nn = pos.size()/time.size();
  //int nb = bond.size();
  //int na = face.size();
  //cout<<"face size "<<face.size()<<"ts "<<time.size()<<endl;
  stringstream output_filename;
  output_filename << fw_name;
  ofstream of;
  of.open(output_filename.str().c_str(),std::ofstream::app);
  for (i=0;i<time.size();i++){
    for (j=0;j<nn;j++)
      of << pos[i*nn+j][0]<<" "<< pos[i*nn+j][1]<<" "<<pos[i*nn+j][2]<<endl;
  }
  of.close();
}
