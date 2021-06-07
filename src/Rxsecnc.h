#include <string>
//---------------------------------------
// routine to calculate NC xsecs---------
// Author : R. Aggarwal (Dated : 15 Feb 2020)
// gives reduced xsection................ 
//----------------------------------------

// Following variables should be provided as external cards
// if user wants to change them (in future)
double ZMass = 91.1876;
double WMass = 80.398; 
double AlphaEM = 7.297352570e-03;
double GFermi=1.16637e-05; 
double TopMass=171.2; 
double  BottomMass=4.20;

double  Vub= 41.2e-3;
double  Vcb=3.93e-3; 

double  Sin2ThetaW=0.23127; 
double  Sin2ThetaC=0.05; 
double  vu=0.19164;
double  vd=-0.34582;
double  ve=-0.03746;

double   au=0.5;
double  ad=-0.5;
double  ae=-0.5;

double pi =3.14159;
// variables for Rxsecnc_xq2 
double  sqrt_s = 318.1; // should be external card for calculation of y
int Lepcharge =1; // should be external card for calculation of Y
int n_points = 429;// should be external card for user




double  F2_LO(double x, double q2){

  double pz= q2/((ZMass*ZMass+q2)*(4*(Sin2ThetaW)*(1-Sin2ThetaW)));

  double Au=4.0/9.0 -2*pz*(2.0/3.0)*(vu)*(ve) + pz*pz*(ve*ve+ae*ae)*(vu*vu+au*au);
  double Ad=1.0/9.0 -2*pz*(-1.0/3.0)*(vd)*(ve) + pz*pz*(ve*ve+ae*ae)*(vd*vd+ad*ad);
  double weights[] = { 0.,Ad,Au,Ad,Au,Ad, 0.,  Ad, Au, Ad, Au, Ad, 0.}; // as in HERAPDF (top set to 0)

  
  double * sum = new double[1];
  double xival[] = {x};
  double q2ival[] = {q2};
  
  QCDNUM::zmstfun(2,weights,xival,q2ival,sum, 1,0);
  
  return sum[0];
}


double  xF3_LO(double x, double q2){

  double pz= q2/((ZMass*ZMass+q2)*(4*(Sin2ThetaW)*(1-Sin2ThetaW)));

  double Bu=-2*(2.0/3.0)*au*ae*pz+4*au*ae*vu*ve*pz*pz;
  double Bd=-2*(-1.0/3.0)*ad*ae*pz+4*ad*ae*vd*ve*pz*pz;
  // double weights[] = { -Bu,-Bd,-Bu,-Bd,-Bu,-Bd, 0.,  Bd, Bu, Bd, Bu, Bd, Bu };
  double weights[] = { 0.,-Bd,-Bu,-Bd,-Bu,-Bd, 0.,  Bd, Bu, Bd, Bu, Bd, 0. };// as in HERAPDF (top set to 0)

  double * sum = new double[1];
  double xival[] = {x};
  double q2ival[] = {q2};
  QCDNUM::zmstfun(3,weights,xival,q2ival,sum, 1,0);
  return sum[0];

}



double FL_LO(double x, double q2){
  double pz= q2/((ZMass*ZMass+q2)*(4*(Sin2ThetaW)*(1-Sin2ThetaW)));

  double Au=4.0/9.0 -2*pz*(2.0/3.0)*(vu)*(ve) + pz*pz*(ve*ve+ae*ae)*(vu*vu+au*au);
  double Ad=1.0/9.0 -2*pz*(-1.0/3.0)*(vd)*(ve) + pz*pz*(ve*ve+ae*ae)*(vd*vd+ad*ad);
  
  double weights[] = { 0.,Ad,Au,Ad,Au,Ad, 0.,  Ad, Au, Ad, Au, Ad, 0.};// as in HERAPDF (top set to 0)
  double * sum = new double[1];
  double xival[] = {x};
  double q2ival[] = {q2};
  QCDNUM::zmstfun(1,weights,xival,q2ival,sum, 1,0);
  return sum[0];
  
}




double Rxsecnc_xq2_i(double x, double q2){
  double  Rxsec =-1.;
  double  y =0.04;
  
  y= q2/sqrt_s/sqrt_s/x;
  double Y_plus = 1+(1-y)*(1-y);
  double Y_minus = 1-(1-y)*(1-y);

  double  F2 = F2_LO(x,q2);
  double  xF3 = xF3_LO(x,q2);
  double FL = FL_LO(x,q2);
  
  if(Lepcharge ==1)Rxsec =  F2-(Y_minus/Y_plus)*xF3-(y*y/Y_plus)*FL;
  else if(Lepcharge ==-1)Rxsec =  F2+(Y_minus/Y_plus)*xF3-(y*y/Y_plus)*FL;
  return Rxsec;
  
}


double* Rxsecnc_xq2(double *x, double *q2){
  
  double * Rxsec = new double[ n_points];
  for(int ix=0; ix<n_points; ix++) {
    Rxsec[ix] = Rxsecnc_xq2_i(x[ix],q2[ix]);
  }
  return Rxsec;
  
}
//===============================================================
// double differential NC cross sections (without pol and order cards)
//===============================================================

double NCPropagator(double q2, double x){
   double y= q2/sqrt_s/sqrt_s/x;
  double Yplus=(1+pow((1-y),2));
  double alpha=AlphaEM;
   double conversion_factor= 0.3894e9;  // GeV^2 -> pb^2
  // double conversion_factor= 1.;//0.3894e9;  // GeV^2 -> pb^2
   return conversion_factor*2*pi*Yplus*alpha*alpha/(x*(q2*q2));
   // return conversion_factor*2*pi*alpha*alpha/(x*(q2*q2));
}

// R. Aggarwal (10 March 2021)..........
// routine for calculation of dd xsec as a function of q2 and x only
// as a fit in the spline calculation for integration
// for y, sqrt_s should be supplied as required, externally
// lepton_charge should be supplied as required, externally
// if in future pol should be supplied, routine should be modified



double dd_xsecnc_xq2_i(double x, double q2){

  double dd_xsec = -1.0;
  dd_xsec = NCPropagator(q2,x)*Rxsecnc_xq2_i(x,q2);
  return dd_xsec;
}




double* dd_xsecnc_xq2(double *x, double *q2){
  
  double * xsec = new double[ n_points];
  for(int ix=0; ix<n_points; ix++) {
    xsec[ix] = dd_xsecnc_xq2_i(x[ix],q2[ix]);
  }
  return xsec;
  
}


double* Rxsecnc(double *x, double *q2, double cmsE, int npoints, string order, string Pol,int Lep_charge){
  // Note : Pol switch not used at the moment ......
  // routine will have to be modified when polarisation of
  // lepton is used as non-zero (if in future)

  // order switch is also not used here ......
  //order switch will change the order of calculations in the qcdnum main code
  //this has to be an external card for the end-user
  n_points =npoints;
  Lepcharge = Lep_charge;
  sqrt_s = cmsE;
  double * Rxsec = new double[n_points];
  Rxsec = Rxsecnc_xq2(x,q2);
  return Rxsec;
  
}

double* dd_xsecnc(double *x, double *q2, double cmsE, int npoints, string order, string Pol,int Lep_charge){
  // Note : Pol switch not used at the moment ......
  // routine will have to be modified when polarisation of
  // lepton is used as non-zero (if in future)

  // order switch is also not used here ......
  //order switch will change the order of calculations in the qcdnum main code
  //this has to be an external card for the end-user
  n_points =npoints;
  Lepcharge = Lep_charge;
  sqrt_s = cmsE;
  double * xsec = new double[n_points];
  xsec = dd_xsecnc_xq2(x,q2);
  return xsec;
  
}


double fun_int(int * ipx, int *ipq, bool *first){

  int ix =*ipx;
  int iq = *ipq;

  double q2 = QCDNUM::qfrmiq(iq);
  double x =QCDNUM::xfrmix(ix);

  // double xsec =  1./(x);//*q2*q2);
  double xsec = dd_xsecnc_xq2_i(x,q2);//x*q2;
  // cout <<ix <<","<<iq<< " = "<< x << " , " <<q2 << " Xsec " << xsec << endl;
  return  xsec;

}

// routine for calculation of dd xsec as a function of iq2 and ix only
// as a fit in the spline calculation for integration
