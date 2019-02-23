#include "stdafx.h"
#include "DataExtractor.h"
#include "FindFermi.h"
using namespace std;

//int FindFermi::func(Ipp64f * params, Ipp64f * argkz, Ipp64f * argCos, Ipp64f * argSin, Ipp64f * r, int length, Ipp64f * temp, Ipp64f * out)
int FindFermi::func(double *params, Ipp64f * argkz, Ipp64f * argCos, Ipp64f * argSin, Ipp64f *r, int length, Ipp64f *temp, Ipp64f *out) {

	ippsMulC_64f(argCos, 3.747665940, temp, length);
	ippsMul_64f_I(r, temp, length);
	vdCos(length, temp, &temp[1 * length]); // cos cos
	ippsMulC_64f(argSin, 3.747665940, temp, length);
	ippsMul_64f_I(r, temp, length);
	vdCos(length, temp, &temp[2 * length]); // cos sin
	ippsMulC_64f(argCos, 3.747665940 / 2, temp, length);
	ippsMul_64f_I(r, temp, length);
	vdCos(length, temp, &temp[3 * length]); // cos cos/2
	ippsMulC_64f(argSin, 3.747665940 / 2, temp, length);
	ippsMul_64f_I(r, temp, length);
	vdCos(length, temp, &temp[4 * length]); // cos sin/2
	ippsMulC_64f(argCos, 3.747665940 * 2, temp, length);
	ippsMul_64f_I(r, temp, length);
	vdCos(length, temp, &temp[5 * length]); // cos 2 cos
	ippsMulC_64f(argSin, 3.747665940 * 2, temp, length);
	ippsMul_64f_I(r, temp, length);
	vdCos(length, temp, &temp[6 * length]); // cos 2 sin
	ippsMulC_64f(argkz, 0.5, &temp[8 * length], length); //kzc / 2
	vdCos(length, &temp[8 * length], &temp[7 * length]); // cos kzc/2
	vdCos(length, &temp[8 * length], &temp[9 * length]); // cos kzc ***made same as above, cos kzc/2

	ippsAdd_64f(&temp[5 * length], &temp[6 * length], temp, length);// param 5
	ippsMulC_64f(temp, -35164.83516*params[5 - 1], out, length);

	ippsMul_64f(&temp[1 * length], &temp[2 * length], temp, length);// param 4
	ippsMulC_64f_I(-35164.83516 * 2 * params[4 - 1], temp, length);
	ippsAdd_64f_I(temp, out, length);

	ippsAdd_64f(&temp[1 * length], &temp[2 * length], temp, length);// param 3
	ippsMulC_64f_I(-35164.83516 * params[3 - 1], temp, length);
	ippsAdd_64f_I(temp, out, length);
	ippsAddC_64f_I(35164.83516 / 2 * params[2 - 1], out, length);// param 2
	ippsSub_64f(&temp[2 * length], &temp[1 * length], temp, length);// param 6
	ippsSqr_64f_I(temp, length); //square
	ippsMul_64f_I(&temp[3 * length], temp, length); // mult by cos cos/2
	ippsMul_64f_I(&temp[4 * length], temp, length); // mult by cos sin/2
	ippsMul_64f_I(&temp[7 * length], temp, length); // mult by cos  kz/2
	ippsMulC_64f_I(-35164.83516 * params[6 - 1], temp, length);
	ippsMulC_64f_I(-35164.83516 *params[7 - 1], &temp[9 * length], length);
	ippsAdd_64f_I(temp, out, length);
	ippsAdd_64f_I(&temp[9 * length], out, length);

	return 0;
}
int FindFermi::interfunc(Ipp64f *theta, Ipp64f *kx, Ipp64f *ky, int Laylength, int clength, int flength, Ipp64f *temp, Ipp64f *temp2, Ipp64f *out) {
	int count;
	for (int i = 0; i < clength; ++i) {
		ippsSub_64f(&kx[i*Laylength], &kx[i*Laylength + 1], &temp[2 * 0 * (Laylength - 1)], Laylength - 1);//kx[i+1]-kx[i]
		ippsSub_64f(&ky[i*Laylength], &ky[i*Laylength + 1], &temp[(2 * 0 + 1) * (Laylength - 1)], Laylength - 1);//ky[i+1]-ky[i]
		ippsMul_64f_I(&temp[2 * 0 * (Laylength - 1)], &temp[2 * 0 * (Laylength)], Laylength - 1);//(kx[i+1]-kx[i])^2
		ippsMul_64f_I(&temp[(2 * 0 + 1) * (Laylength - 1)], &temp[(2 * 0 + 1) * (Laylength - 1)], Laylength - 1);//(ky[i+1]-ky[i])^2
		ippsAdd_64f_I(&temp[2 * 0 * (Laylength - 1)], &temp[(2 * 0 + 1) * (Laylength - 1)], Laylength - 1);//(kx[i+1]-kx[i])^2+(ky[i+1]-ky[i])^2
		ippsSqrt_64f_I(&temp[(2 * 0 + 1) * (Laylength - 1)], Laylength - 1);//Sqrt[(kx[i+1]-kx[i])^2+(ky[i+1]-ky[i])^2]
		for (int j = 0; j < Laylength; ++j) {
			temp[j] = 0;

			for (int k = 0; k < j; ++k) {
				temp[j] = temp[j] + temp[(2 * 0 + 1) * (Laylength - 1) + k];//cumulative sum
			}
		}
		for (int j = 0; j < flength; ++j) {
			temp2[j] = temp[Laylength - 1] / flength * j;
		}
		count = 1;
		out[i * flength] = 0;
		for (int k = 0; k < Laylength; ++k) {
			if (count > flength - 1) break;
			if (temp[k] >= temp2[count]) {
				out[i * flength + count] = theta[k - 1] + (theta[k] - theta[k - 1])*(temp2[count] - temp[k - 1]) / (temp[k] - temp[k - 1]);
				++count;

			}
		}




	}
}

FindFermi::FindFermi(double * param,int cedv,int gridNpoint)
{
	//DataExtractor extractor(name);
	//Ipp64f params[9] = { 0.074, 475, 525, -60, 16, 1000, 0.5, 17, 8 };
	UpdatePar(param);
	fineN = 3000;//innitial grid inplane
	gridN = gridNpoint;//actual grid inplane
	cdevs = cedv;//Kz grid
	//Ipp64f * starts = extractor.getDataArray();
	//int nPoints = floor((extractor.getNumberOfLines()) / 2);
	nPoints = gridN * cdevs;
	
	nfinepoint = (fineN)*cdevs;
	tol = 1E-10;
	subMaxR = new Ipp64f[1];
	temp1 = new Ipp64f[20 * nfinepoint];
	argCos = new Ipp64f[nfinepoint];
	argSin = new Ipp64f[nfinepoint];
	argkz = new Ipp64f[nfinepoint];
	gtheta = new Ipp64f[nfinepoint];
	kz = new Ipp64f[nfinepoint];
	funcval = new Ipp64f[nfinepoint];
	dfunc1 = new Ipp64f[nfinepoint];
	dfunc2 = new Ipp64f[nfinepoint];
	dfunc = new Ipp64f[nfinepoint];
	r = new Ipp64f[nfinepoint];
	rtemp = new Ipp64f[nfinepoint];
	kx = new Ipp64f[nfinepoint];
	ky = new Ipp64f[nfinepoint];
	finthe = new Ipp64f[nfinepoint];
	//Ipp64f *startpoint = new Ipp64f[3 * nPoints];
	for (int i = 0; i < cdevs; i++)
	{
		for (int j = 0; j < fineN; j++)
		{
			/*cout << starts[2 * i] << " " << starts[2 * i + 1] << endl;*/
			gtheta[i*fineN + j] = 2 * 3.1415926 / (fineN - 1) * j;
			kz[i*fineN + j] = -2 * 3.1415926 / 13.2 + 4 * 3.1415926 / 13.2 / cdevs / 2 + 4 * 3.1415926 / 13.2 / cdevs * i;
			r[i*fineN + j] = 0.3;

		}
	}
	vdSin(nfinepoint, gtheta, argSin);
	vdCos(nfinepoint, gtheta, argCos);
	ippsMulC_64f(kz, 13.2, argkz, nfinepoint);
	//cout << "parameter: ";
	//for (int i = 0; i < 10; ++i) {
	//	cout << params[i] << ", ";
	//}
	//cout << endl;
	
}

int FindFermi::UpdatePar(double * param)
{
	
	for (int i = 0; i < 9; ++i) {
		params[i] = param[i];
	}
	//cout << nPoints << endl;
	return 0;
}

int FindFermi::PrintPar()
{
	cout << "parameter: ";
	for (int i = 0; i < 10; ++i) {
		cout << params[i] << ", ";
	}
	cout << endl;
	return 0;
}



int FindFermi::ReturnStart(Ipp64f * startpoint)
{
	//cout << nPoints << endl;
	for (int i = 0; i < 50; i++)
	{  /*vdSin(nPoints, theta, argSin);
		vdCos(nPoints, theta, argCos);
		ippsMulC_64f_I(3.74767, argSin, nPoints);
		ippsMulC_64f_I(3.74767, argCos, nPoints);
		ippsMulC_64f(kz, 3.3, argkz, nPoints);
		func(argkz, argCos, argSin, r, nPoints, temp1, temp2, temp3, funcval);

		ippsAddC_64f(r, 1E-5, rtemp, nPoints);
		func(argkz, argCos, argSin, rtemp, nPoints, temp1, temp2, temp3, dfunc1);
		ippsAddC_64f(r, -1E-5, rtemp, nPoints);
		func(argkz, argCos, argSin, rtemp, nPoints, temp1, temp2, temp3, dfunc2);
		ippsMulC_64f_I(-1, dfunc2, nPoints);
		ippsAdd_64f(dfunc1, dfunc2, dfunc, nPoints);
		ippsDivC_64f_I(2E-5, dfunc, nPoints);

		ippsDiv_64f_I(dfunc, funcval, nPoints);
		ippsSub_64f_I(funcval, r, nPoints);*/

		//caxis lattice parameter =3.3  k_z * c
		func(params, argkz, argCos, argSin, r, nfinepoint, temp1, funcval);
		//func(nPoints, temp1,funcval);
		//cout << argCos[1] << endl;
		ippsAbs_64f(funcval, temp1, nfinepoint);
		ippsMax_64f(temp1, nfinepoint, subMaxR);
		//cout << i << endl;
		//cout << funcval[0] << endl;
		//cout << subMaxR[0] << endl;
		if (subMaxR[0] <= tol)
		{

			break;
		}
		ippsAddC_64f(r, 1E-6, rtemp, nfinepoint);
		func(params, argkz, argCos, argSin, rtemp, nfinepoint, temp1, dfunc1);
		//func(nPoints, temp1,dfunc1);
		ippsAddC_64f(r, -1E-6, rtemp, nfinepoint);
		func(params, argkz, argCos, argSin, rtemp, nfinepoint, temp1, dfunc2);
		//func(nPoints, temp1, dfunc2);
		ippsMulC_64f_I(-1, dfunc2, nfinepoint);
		ippsAdd_64f(dfunc1, dfunc2, dfunc, nfinepoint);
		ippsDivC_64f_I(2E-6, dfunc, nfinepoint);

		ippsDiv_64f_I(dfunc, funcval, nfinepoint);
		ippsSub_64f_I(funcval, r, nfinepoint);

	}

	//for (int i = 0; i < nPoints; i++)
	//{
	//	cout << r[i] << endl;
	//}
	//ippsMul_64f(r, argCos, kx, nPoints);
	//ippsMul_64f(r, argSin, ky, nPoints);

	ippsMul_64f(r, argCos, kx, nfinepoint);
	ippsMul_64f(r, argSin, ky, nfinepoint);
	interfunc(gtheta, kx, ky, fineN, cdevs, gridN, temp1, rtemp, finthe);

	//get the grid for the uniform grid

	for (int i = 0; i < cdevs; i++)
	{
		for (int j = 0; j < gridN; j++)
		{

			gtheta[i*gridN + j] = finthe[i*gridN + j];
			kz[i*gridN + j] = -2 * 3.1415926 / 13.2 + 4 * 3.1415926 / 13.2 / cdevs / 2 + 4 * 3.1415926 / 13.2 / cdevs * i;
			r[i*gridN + j] = 0.3;
		}
	}
	
	vdSin(nPoints, gtheta, argSin);
	vdCos(nPoints, gtheta, argCos);
	ippsMulC_64f(kz, 13.2, argkz, nPoints);
	for (int i = 0; i < 40; i++)
	{


		func(params, argkz, argCos, argSin, r, nPoints, temp1, funcval);
		//func(nPoints, temp1,funcval);
		//cout << argCos[1] << endl;
		ippsAbs_64f(funcval, temp1, nPoints);
		ippsMax_64f(temp1, nPoints, subMaxR);
		//cout << i << endl;
		//cout << funcval[0] << endl;
		//cout << subMaxR[0] << endl;
		if (subMaxR[0] <= tol)
		{

			break;
		}
		ippsAddC_64f(r, 1E-6, rtemp, nPoints);
		func(params, argkz, argCos, argSin, rtemp, nPoints, temp1, dfunc1);
		//func(nPoints, temp1,dfunc1);
		ippsAddC_64f(r, -1E-6, rtemp, nPoints);
		func(params, argkz, argCos, argSin, rtemp, nPoints, temp1, dfunc2);
		//func(nPoints, temp1, dfunc2);
		ippsMulC_64f_I(-1, dfunc2, nPoints);
		ippsAdd_64f(dfunc1, dfunc2, dfunc, nPoints);
		ippsDivC_64f_I(2E-6, dfunc, nPoints);

		ippsDiv_64f_I(dfunc, funcval, nPoints);
		ippsSub_64f_I(funcval, r, nPoints);










	}



	ippsMul_64f(r, argCos, kx, nPoints);
	ippsMul_64f(r, argSin, ky, nPoints);

	for (int i = 0; i < nPoints; i++) {
		startpoint[i * 3] = kx[i];
		startpoint[i * 3 + 1] = ky[i];
		startpoint[i * 3 + 2] = kz[i];
	}
	//ofstream fout;
	//fout.open("FindFermi.dat");
	//fout.precision(15);

	//for (int i = 0; i < nPoints; i++) {

	//	fout << startpoint[i * 3] << "\t" << startpoint[i * 3 + 1] << "\t" << startpoint[i * 3 + 2] << endl;
		//	fout << theta[i] << "\t" << r[i] << "\t" << kz[i] << endl;
	//		//cout << kx[i] << "\t" << ky[i] << "\t" << kz[i] << endl;
	//}

	//fout.close();
	//while (true);
	return 0;
}

int FindFermi::ReturnNumPoint()
{
	cout << "nPoints = " << &nPoints << "  " << nPoints << endl;
	return nPoints;
}


FindFermi::~FindFermi()
{
	delete subMaxR;
	delete temp1;
	delete argCos;
	delete argSin;
	delete gtheta;
	delete argkz;
	delete kz;
	delete funcval;
	delete dfunc1;
	delete dfunc2;
	delete dfunc;
	delete r;
	delete rtemp;
	delete kx;
	delete ky;
	delete finthe;
}
