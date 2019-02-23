#pragma once
class FindFermi
{
private:
	
	int func(double *params, Ipp64f * argkz, Ipp64f * argCos, Ipp64f * argSin, Ipp64f *r, int length, Ipp64f *temp, Ipp64f *out);
	//int nPoints;
	Ipp64f *subMaxR;
	Ipp64f *temp1;
	Ipp64f *argCos;
	Ipp64f *argSin;
	Ipp64f *argkz;
	Ipp64f *gtheta;
	Ipp64f *kz;
	Ipp64f *funcval;
	Ipp64f *dfunc1;
	Ipp64f *dfunc2;
	Ipp64f *dfunc;
	Ipp64f *r;
	Ipp64f *rtemp;
	Ipp64f *kx;
	Ipp64f *ky;
	//Ipp64f *startpoint ;
	Ipp64f tol;
	int fineN;
	int gridN;
	int cdevs;
	int nfinepoint;
	Ipp64f *finthe;
	double params[9];
public:
	

	int nPoints;
	FindFermi(double * param,int cdev,int gridNpoint);
	int UpdatePar(double *param);
	int PrintPar();
	int ReturnStart(Ipp64f *startpoint);
	int ReturnNumPoint();
	int interfunc(Ipp64f *theta, Ipp64f *kx, Ipp64f *ky, int Laylength, int clength, int flength, Ipp64f *temp, Ipp64f *temp2, Ipp64f *out);
	virtual ~FindFermi();
};

