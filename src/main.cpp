#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <gmp.h>
#include "math.h"

#include "ht.h"
#include "elliptic.h"
#include "zolo.h"
#include "zolorf.h"

int main() {

  std::cout << std::setprecision(12);

  bool ZOLO=true;
  if (ZOLO) { 
    elliptic e;
    zolo z;
    z.setZolo(ftype(1e-3),ftype(10),21); // construct the Zolotarev approximation
                                         // output roots of z(x)=1 and extrema
//  z.writeZolo(200);
    zolorf zrf(z); // construct factor and partial fraction formulations
                   // output approximation coefficents
    int Nout=200;
    zrf.outputRFZoloRange(Nout,true); // output Nout data points of PF formulation 
    zrf.outputRFZoloRange(Nout,false); // output Nout data points of factor formulation
  }

  bool HT=false; 
  if (HT) {
    ht h;
    std::cout << "setHT" << std::endl;
    int Nbase=17;
    int scale=1; // adjust the range of approximation
    int N=Nbase/scale;
    ftype xlow=1e-3;
    ftype xhigh=20;
    h.setHT(N);
    h.writeFactors();
    std::cout << "writeHT" << std::endl;
    int Nplot=200;
    bool REFLECT=true;
    h.writeHT(Nplot,xlow,xhigh,scale,REFLECT);
  }

  return 0;
}
