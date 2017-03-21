#include <iostream>
#include <cmath>
//#include "TH1.h"
//#include <TMath.h>
//#include <TGraph.h>
//#include <TH2.h>
//#include <TH2D.h>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "gauss.c"
#include "lorentzian.c"
#include "lmrqmin.c"


using namespace std;

int DataSize(ifstream& stream)
{
  int count = 0;
  string line;
  while(!stream.eof())
    {
      std::getline(stream,line);
      count++;
    }
  cout << count;
  return count;
}

double Lorentz(double xval, double a, double xavg)
{
  double yval;
  yval = (1/M_PI)*(a)/(pow((xval-xavg),2) + pow(a,2));
  return yval;
  
}

void Main()
{

  ifstream data;
  data.open("hw3_fitting.dat");
  
  int NumData = DataSize(data);
  
  
  //Coefficients for reading in data
  int num = 0;
  int i = 0;
  double x[NumData],y[NumData],error[NumData];
  while(!data.eof())
    {
      data >> x[i] >> y[i] >> error[i];
      i++;
    }
  /*for(j = 0; j < NumData; j++)*/

  




  
  data.close();  
}
