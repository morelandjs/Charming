#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip> 
#include <cmath>
#include <stdlib.h>
#include <vector>
#include "ParameterReader.h"

using namespace std;

bool sample_production(double rate)
{
  if(drand48()<rate) return true;
  else return false;
}

double sample_pt()
{
  ifstream ptfile;
  ptfile.open("/home/jsm55/Research/Charming/BC2Charm/src/FONNL/FONNL.dat");
  
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

void histogram(ParameterReader* _paraRdr)
{
  ParameterReader* paraRdr = _paraRdr;
  double mass = paraRdr->getVal("mass");
  double eta_min = paraRdr->getVal("eta_min");
  double eta_max = paraRdr->getVal("eta_max");

  ifstream infile;
  ofstream outfile;
  outfile.open("/home/jsm55/Research/Charming/BC2Charm/data/histogram.dat",std::ios_base::trunc);

  for(int i=0; i<1000000; i++){
    double pt = sample_pt();
    double phi = drand48()*2.0*M_PI;
    double eta = eta_min + drand48()*(eta_max-eta_min);
    double px = pt*cos(phi);
    double py = pt*sin(phi);
    double pz = pt*sinh(eta);
    double p = pow(px*px+py*py+pz*pz,0.5);
    double E = pow(p*p+mass*mass,0.5);
    outfile << setprecision(12) << setw(20) << E 
            << setprecision(12) << setw(20) << p 
	    << setprecision(12) << setw(20) << px 
	    << setprecision(12) << setw(20) << py  
	    << setprecision(12) << setw(20) << pz  
	    << endl;
  }
  outfile.close();
}

void makeCharm(ParameterReader* _paraRdr)
{
  ParameterReader* paraRdr = _paraRdr;
  int nevents = paraRdr->getVal("nevents");
  double mass = paraRdr->getVal("mass");
  double eta_min = paraRdr->getVal("eta_min");
  double eta_max = paraRdr->getVal("eta_max");
  double rate = paraRdr->getVal("production_rate");

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
      if(sample_production(rate)){ 
	double pt = sample_pt();
	double phi = drand48()*2.0*M_PI;
	double eta = eta_min + drand48()*(eta_max-eta_min);
	double px = pt*cos(phi);
	double py = pt*sin(phi);
	double pz = pt*sinh(eta);
	double p = pow(px*px+py*py+pz*pz,0.5);
	double E = pow(p*p+mass*mass,0.5);
	outfile << setprecision(5) << setw(13) << 0.0  
	        << setprecision(5) << setw(13) << xx  
		<< setprecision(5) << setw(13) << yy 
	        << setprecision(5) << setw(13) << 0.0
		<< setprecision(12) << setw(20) << E 
		<< setprecision(12) << setw(20) << px 
		<< setprecision(12) << setw(20) << py  
		<< setprecision(12) << setw(20) << pz  
		<< endl;
      }
    }
    while(!infile.eof());
    infile.close();
    outfile.close();
  }
}

int main(int argc, char *argv[])
{
  // read parameters
  ParameterReader paraRdr;
  paraRdr.readFromFile("/home/jsm55/Research/Charming/BC2Charm/src/parameters.dat");
  paraRdr.readFromArguments(argc, argv);
  cout << endl; paraRdr.echo();

  // routines
  //histogram(&paraRdr);
  makeCharm(&paraRdr);

  return 0;
}
