#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip> 
#include <cmath>
#include <stdlib.h>
#include <vector>

using namespace std;

bool sample_production()
{
  double rate = 0.01;
  double rand = drand48();
  if(drand48()<rate) return true;
  else return false;
}

double sample_pt()
{
  ifstream ptfile;
  ptfile.open("/home/jsm55/Research/Charming/BC2Charm/FONNL/FONNL.dat");
  
  double sigma=0.0;

  string line;
  vector<double> sigma_cdf, pt;
  while(getline(ptfile, line)) 
    {
      line.erase(line.begin(), find_if(line.begin(), line.end(), not1(ptr_fun<int, int>(isspace)))); 
      if(line[0] == '#') continue;
      double pt_val, dsigma_val;
      istringstream(line) >> pt_val >> dsigma_val;
      sigma += dsigma_val;
      sigma_cdf.push_back(sigma);
      pt.push_back(pt_val);
    }
  int len = sigma_cdf.size();
  for(int i=0; i<len; i++) sigma_cdf[i] /= sigma;

  double rand = drand48();
  int ipt = 0;
  while(rand > sigma_cdf[ipt]) ipt++;

  return pt[ipt];
}

int main()
{
  int nevents = 2;
  double mass = 1.35;
  double sqrt_s = 2760;

  ifstream infile;
  ofstream outfile;

  char* cpath = "/home/jsm55/Research/Charming";

  char filename[] = "%s/superMC/data/bctable_event_%d.dat";
  char printname[] = "%s/BC2Charm/data/charmdata_event_%d.dat";
  char buffer1[200], buffer2[200];

  for(int n=1; n<=nevents; n++){
    sprintf(buffer1,filename,cpath,n);
    sprintf(buffer2,printname,cpath,n);
    infile.open(buffer1);
    outfile.open(buffer2,std::ios_base::trunc);
    do{	
      double xx, yy;
      infile >> xx >> yy;
      if(sample_production()){ 
	double pt = sample_pt();
	double phi = drand48()*2.0*M_PI;
	double px = pt*cos(phi);
	double py = pt*sin(phi);
	double pz = pow(pow(sqrt_s/2.0,2.0)-pow(mass,2.0),0.5);
	double E = pow(pow(px,2.0)+pow(py,2.0)+pow(pz,2.0)+pow(mass,2.0),0.5);
	outfile << setprecision(5) << setw(5) << xx << "\t" 
		<< setprecision(5) << setw(5) << yy << "\t" 
		<< setprecision(12) << setw(12) << E << "\t" 
		<< setprecision(12) << setw(12) << px << "\t" 
		<< setprecision(12) << setw(12) << py << "\t" 
		<< setprecision(12) << setw(12) << pz << "\t" 
		<< endl;
      }
    }
    while(!infile.eof());
    infile.close();
    outfile.close();
  }

return 0;
}
