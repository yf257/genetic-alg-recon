#include "stdafx.h"
#include "calculateAMRO.h"
using namespace std;




calculateAMRO::calculateAMRO(double *data, double * param, Ipp64f * theta, int cdev, int gridN,int dataLeng,Ipp64f * phi)
{
	
	//phi = 0;
	
	
	_dataLeng = dataLeng/4;
	thetas = new Ipp64f[_dataLeng];
	phis = new Ipp64f[4];
	_data = new Ipp64f[_dataLeng*4];
	updatetheta(theta, _dataLeng);
	updatedata(data, _dataLeng * 4);
	updatephi(phi, 4);
	//cout << dataLeng << endl;
	params = new Ipp64f[9];
	updatepara(param,9);

	//ippsDivC_64f_I(data[0], data, dataLeng);
	//cout << dataLeng << endl;
	condout = new Ipp64f[_dataLeng];
	//tau = .5;
	final = 8*param[0];//time final?
	steps = 500;//number of time steps?
	h = final / steps;
	field45 = 7.91209; // 45 tesla in appropriate units
	Fermi =new FindFermi( param,cdev,gridN);

	//DataExtractor extractor("starts.dat");
	//Ipp64f * starts = extractor.getDataArray();
	//int nPoints = floor((extractor.getNumberOfLines()) / 3);//x,y,z coordinates?

	nPoints = Fermi->nPoints;
	
	starts = new Ipp64f[nPoints * 3];
	//Fermi->ReturnStart(starts);
	output = new Ipp64f[nPoints*steps * 3]; //stores evolution of orbit around Fermi surface
	times = new Ipp64f[steps*nPoints]; //time steps

	//std::clock_t startT;
	//duration;

	field = new Ipp64f[3];
	vzStorage = new Ipp64f[steps*nPoints];
	vz0Storage = new Ipp64f[nPoints];
	DOS = new Ipp64f[nPoints];
	ones = new Ipp64f[nPoints];//for inverting
	ippsSet_64f(1, ones, nPoints);
	taus = new Ipp64f[nPoints];//phi dependent taus

	vx = new Ipp64f[nPoints];
	vy = new Ipp64f[nPoints];
	vz = new Ipp64f[nPoints];

	argx = new Ipp64f[nPoints];
	argy = new Ipp64f[nPoints];
	argz = new Ipp64f[nPoints];

	tempx = new Ipp64f[20 * nPoints];
	tempy = new Ipp64f[20 * nPoints];
	tempz = new Ipp64f[20 * nPoints];

	k1x = new Ipp64f[nPoints];
	k1y = new Ipp64f[nPoints];
	k1z = new Ipp64f[nPoints];
	k2x = new Ipp64f[nPoints];
	k2y = new Ipp64f[nPoints];
	k2z = new Ipp64f[nPoints];
	k3x = new Ipp64f[nPoints];
	k3y = new Ipp64f[nPoints];
	k3z = new Ipp64f[nPoints];
	k4x = new Ipp64f[nPoints];
	k4y = new Ipp64f[nPoints];
	k4z = new Ipp64f[nPoints];//what are these? why do we need k1-4?
	tempdif = new Ipp64f[_dataLeng*3];
	//Fermi->ReturnStart(starts);
	exptau = new Ipp64f[steps*nPoints];
	minDos = 0;
    maxDos = 0;
}
calculateAMRO::~calculateAMRO()
{
	delete params;
	delete thetas;
	
	delete condout;
	delete Fermi;
	delete output;
	delete starts;
	delete times;
	delete field;
	delete vzStorage;
	delete vz0Storage;
	delete DOS;
	delete ones;//for inverting
	
	delete taus;//phi dependent taus

	delete vx;
	delete vy;
	delete vz;

	delete argx;
	delete argy;
	delete argz;

	delete tempx;
	delete tempy;
	delete tempz;

	delete k1x;
	delete k1y;
	delete k1z;
	delete k2x;
	delete k2y;
	delete k2z;
	delete k3x;
	delete k3y;
	delete k3z;
	delete k4x;
	delete k4y;
	delete k4z;//what are these? why do we need k1-4?
	delete tempdif;
	delete _data;
	delete phis;
	delete exptau;

}
int calculateAMRO::updatedata(Ipp64f *data, int Length) {
	for (int i = 0; i < Length; ++i) {
		_data[i] = data[i];
	}
	return 0;
}
int calculateAMRO::updatetheta(Ipp64f * theta, int Length)
{
	for (int i = 0; i < Length; ++i) {
		thetas[i] = theta[i];
	}
	return 0;
}

int calculateAMRO::updatepara(double * param, int Length)
{

	for (int i = 0; i < Length; ++i) {
		params[i] = param[i];
		
	}

	return 0;
}
int calculateAMRO::updatephi(Ipp64f * phi, int Length)
{

	for (int i = 0; i < Length; ++i) {
		phis[i] = phi[i];

	}

	return 0;
}

Ipp64f calculateAMRO::returnvalue(double * param)
{
	
	total = 0;
	sum_resdual = 0;
	std::clock_t startT = std::clock();
	updatepara(param, 9);

	Fermi->UpdatePar(param);
	//cout << Fermi->nPoints << endl;
	Fermi->ReturnStart(starts);
	//Fermi->ReturnStart(starts);
	for (int phi_i = 0; phi_i < 4; phi_i++)
	{

		for (int j = 0; j < nPoints; j++) { //initialize Fermi surface to starting grid, only needs to be done once
			output[0 * nPoints + j] = starts[j * 3 + 0];
			output[1 * nPoints + j] = starts[j * 3 + 1];
			output[2 * nPoints + j] = starts[j * 3 + 2];
		}

		ippsCopy_64f(&output[nPoints * (0)], argx, nPoints);//initial velocities for DOS calc;
		ippsCopy_64f(&output[nPoints * (1)], argy, nPoints);
		ippsCopy_64f(&output[nPoints * (2)], argz, nPoints);
		veloX(params, argx, argy, argz, nPoints, tempx, vx); //velocities for DOS are stored in vx, vy, and vz buffers.
		veloY(params, argx, argy, argz, nPoints, tempy, vy);
		veloZ(params, argx, argy, argz, nPoints, tempz, vz);

		ippsSqr_64f_I(vx, nPoints);//in-place square of velocities
		ippsSqr_64f_I(vy, nPoints);
		ippsSqr_64f_I(vz, nPoints);

		ippsAdd_64f(vx, vy, tempx, nPoints);//add all square velocities
		ippsAdd_64f_I(vz, tempx, nPoints);
		ippsSqrt_64f_I(tempx, nPoints);//square root
		ippsDiv_64f(tempx, ones, DOS, nPoints);
		ippsMax_64f(DOS, nPoints, &maxDos);
		ippsMin_64f(DOS, nPoints, &minDos);

		/*argx[0] = 1;
		argy[0] = 2;
		argz[0] = 1;
		veloZ(params, argx, argy, argz, 1, tempz, vz);
		cout << vz[0] << endl;
		*/
		for (int th = 0; th < _dataLeng; th++) {

			//for (int p = 0; p < steps; p++) { //re-initialize times SLOW STEP CREATE TEMP VARIABLE
										  //times[p] = -p;//why -p?
			//	ippsSet_64f(-p, &times[nPoints * p], nPoints);
			//}
			ippsSet_64f(-1, times, steps*nPoints);
			ippsMulC_64f_I(h, times, steps*nPoints);//time stamps

			field[0] = field45 * sin(thetas[th])*cos(phis[phi_i]);  //set field(phi == 0?)
			field[1] = field45 * sin(thetas[th])*sin(phis[phi_i]);
			field[2] = field45 * cos(thetas[th]);
			//for (int z = 0; z < nPoints; z++) {
			//	if (th == 0 ) {
			//		cout << setprecision(20) << times[nPoints*2 + z] << " " << endl;
			//	}
			//}


			for (int i = 1; i < steps; i++) {

				ippsCopy_64f(&output[nPoints * (3 * (i - 1) + 0)], argx, nPoints);//copy arguments for k1;
				ippsCopy_64f(&output[nPoints * (3 * (i - 1) + 1)], argy, nPoints);
				ippsCopy_64f(&output[nPoints * (3 * (i - 1) + 2)], argz, nPoints);
				veloX(params, argx, argy, argz, nPoints, tempx, vx); //calculate velocities;
				veloY(params, argx, argy, argz, nPoints, tempy, vy);
				veloZ(params, argx, argy, argz, nPoints, tempz, vz);
				ippsCopy_64f(vz, &vzStorage[nPoints * (i - 1)], nPoints);//store vz for conductivity later

				//taufun(params, argx, argy, nPoints, tempx, taus, ones);// calculate k dependent tau
				taufundos(params, minDos, maxDos, argx, argy, argz, nPoints, tempx, taus, ones);
				ippsDiv_64f_I(taus, &times[nPoints * (i - 1)], nPoints);
				//ippsDivC_64f_I(tau, &times[nPoints * (i - 1)], nPoints);
				//ippsExp_64f_I(&times[nPoints * (i - 1)], nPoints);
				//ippsMulC_64f_I((1E-12) * h, &times[nPoints * (i - 1)], nPoints);

				//taufun(params, argx, argy, nPoints, tempx, taus);
				//for (int z = 0; z < nPoints; z++) {
				//	if (th == 0 && i == 30) {
				//		cout << setprecision(20) << taus[z] << " " << endl;
				//	}
				//}
				/*for (int z = 0; z < nPoints; z++) {
				if (th == 1 && i == 1) {
				cout << setprecision(20)<< argx[z] << " " << endl;
				}
				}*/


				fx(field, vy, vz, nPoints, tempx, k1x); //calculate evolution in k and store in k1
				fy(field, vx, vz, nPoints, tempy, k1y);
				fz(field, vx, vy, nPoints, tempz, k1z);

				ippsMulC_64f(k1x, h / 2, tempx, nPoints); //prep evolved k step for k2
				ippsMulC_64f(k1y, h / 2, tempy, nPoints);
				ippsMulC_64f(k1z, h / 2, tempz, nPoints);
				ippsAdd_64f(tempx, &output[nPoints * (3 * (i - 1) + 0)], argx, nPoints); //add step to previous k point, load into arguments for k2;
				ippsAdd_64f(tempy, &output[nPoints * (3 * (i - 1) + 1)], argy, nPoints);
				ippsAdd_64f(tempz, &output[nPoints * (3 * (i - 1) + 2)], argz, nPoints);
				veloX(params, argx, argy, argz, nPoints, tempx, vx); //calculate velocities;
				veloY(params, argx, argy, argz, nPoints, tempy, vy);
				veloZ(params, argx, argy, argz, nPoints, tempz, vz);
				fx(field, vy, vz, nPoints, tempx, k2x); //calculate evolution in k and store in k2
				fy(field, vx, vz, nPoints, tempy, k2y);
				fz(field, vx, vy, nPoints, tempz, k2z);

				ippsMulC_64f(k2x, h / 2, tempx, nPoints); //prep evolved k step for k3
				ippsMulC_64f(k2y, h / 2, tempy, nPoints);
				ippsMulC_64f(k2z, h / 2, tempz, nPoints);
				ippsAdd_64f(tempx, &output[nPoints * (3 * (i - 1) + 0)], argx, nPoints); //add step to previous k point, load into arguments for k3;
				ippsAdd_64f(tempy, &output[nPoints * (3 * (i - 1) + 1)], argy, nPoints);
				ippsAdd_64f(tempz, &output[nPoints * (3 * (i - 1) + 2)], argz, nPoints);
				veloX(params, argx, argy, argz, nPoints, tempx, vx); //calculate velocities;
				veloY(params, argx, argy, argz, nPoints, tempy, vy);
				veloZ(params, argx, argy, argz, nPoints, tempz, vz);
				fx(field, vy, vz, nPoints, tempx, k3x); //calculate evolution in k and store in k3
				fy(field, vx, vz, nPoints, tempy, k3y);
				fz(field, vx, vy, nPoints, tempz, k3z);

				ippsMulC_64f(k3x, h, tempx, nPoints); //prep evolved k step for k4
				ippsMulC_64f(k3y, h, tempy, nPoints);
				ippsMulC_64f(k3z, h, tempz, nPoints);
				ippsAdd_64f(tempx, &output[nPoints * (3 * (i - 1) + 0)], argx, nPoints); //add step to previous k point, load into arguments for k4;
				ippsAdd_64f(tempy, &output[nPoints * (3 * (i - 1) + 1)], argy, nPoints);
				ippsAdd_64f(tempz, &output[nPoints * (3 * (i - 1) + 2)], argz, nPoints);
				veloX(params, argx, argy, argz, nPoints, tempx, vx); //calculate velocities;
				veloY(params, argx, argy, argz, nPoints, tempy, vy);
				veloZ(params, argx, argy, argz, nPoints, tempz, vz);
				fx(field, vy, vz, nPoints, tempx, k4x); //calculate evolution in k and store in k4
				fy(field, vx, vz, nPoints, tempy, k4y);
				fz(field, vx, vy, nPoints, tempz, k4z);

				ippsMulC_64f_I(2, k2x, nPoints); //scale k2
				ippsMulC_64f_I(2, k2y, nPoints);
				ippsMulC_64f_I(2, k2z, nPoints);
				ippsMulC_64f_I(2, k3x, nPoints); //scale k3
				ippsMulC_64f_I(2, k3y, nPoints);
				ippsMulC_64f_I(2, k3z, nPoints);

				ippsAdd_64f(k1x, k2x, tempx, nPoints); //add k1 + k2 to temp
				ippsAdd_64f(k1y, k2y, tempy, nPoints);
				ippsAdd_64f(k1z, k2z, tempz, nPoints);

				ippsAdd_64f_I(k3x, tempx, nPoints); //add in k3
				ippsAdd_64f_I(k3y, tempy, nPoints);
				ippsAdd_64f_I(k3z, tempz, nPoints);

				ippsAdd_64f_I(k4x, tempx, nPoints); //add in k4
				ippsAdd_64f_I(k4y, tempy, nPoints);
				ippsAdd_64f_I(k4z, tempz, nPoints);

				ippsMulC_64f_I(h / 6, tempx, nPoints); //scale the entire sum
				ippsMulC_64f_I(h / 6, tempy, nPoints); //scale the entire sum
				ippsMulC_64f_I(h / 6, tempz, nPoints); //scale the entire sum

				ippsAdd_64f(&output[nPoints * (3 * (i - 1) + 0)], tempx, &output[nPoints * (3 * i + 0)], nPoints); //add sum to previous output and store
				ippsAdd_64f(&output[nPoints * (3 * (i - 1) + 1)], tempy, &output[nPoints * (3 * i + 1)], nPoints);
				ippsAdd_64f(&output[nPoints * (3 * (i - 1) + 2)], tempz, &output[nPoints * (3 * i + 2)], nPoints);
			}


			ippsCopy_64f(&output[nPoints * (3 * (steps - 1) + 0)], argx, nPoints);//get velocity for last point
			ippsCopy_64f(&output[nPoints * (3 * (steps - 1) + 1)], argy, nPoints);
			ippsCopy_64f(&output[nPoints * (3 * (steps - 1) + 2)], argz, nPoints);
			veloZ(params, argx, argy, argz, nPoints, tempz, vz);
			ippsCopy_64f(vz, &vzStorage[nPoints * (steps - 1)], nPoints);

			//taufun(params, argx, argy, nPoints, tempx, taus, ones);// calculate k dependent tau for last point
			taufundos(params, minDos, maxDos, argx, argy, argz, nPoints, tempx, taus, ones);
			ippsDiv_64f_I(taus, &times[nPoints * (steps - 1)], nPoints);
			//ippsDivC_64f_I(tau, &times[nPoints * (steps - 1)], nPoints);
			//ippsExp_64f_I(&times[nPoints * (steps - 1)], nPoints);
			//ippsMulC_64f_I((1E-12)*h, &times[nPoints * (steps - 1)], nPoints);
/*
			ippsCopy_64f(&output[nPoints * (0)], argx, nPoints);//initial velocities for DOS calc;
			ippsCopy_64f(&output[nPoints * (1)], argy, nPoints);
			ippsCopy_64f(&output[nPoints * (2)], argz, nPoints);
			veloX(params, argx, argy, argz, nPoints, tempx, vx); //velocities for DOS are stored in vx, vy, and vz buffers.
			veloY(params, argx, argy, argz, nPoints, tempy, vy);
			veloZ(params, argx, argy, argz, nPoints, tempz, vz);

			ippsSqr_64f_I(vx, nPoints);//in-place square of velocities
			ippsSqr_64f_I(vy, nPoints);
			ippsSqr_64f_I(vz, nPoints);

			ippsAdd_64f(vx, vy, tempx, nPoints);//add all square velocities
			ippsAdd_64f_I(vz, tempx, nPoints);
			ippsSqrt_64f_I(tempx, nPoints);//square root
			ippsDiv_64f(tempx, ones, DOS, nPoints);
*/
			//need to change!!
			  //ippsDiv_64f_I(taus, times, steps);//exponential stuff, negative tau is taken care of in time
			  //ippsExp_64f_I(times, steps);
			  //ippsMulC_64f_I((1E-12)*h, times, steps);

			ippsCopy_64f(&vzStorage[0], vz0Storage, nPoints);//save initial velocity before exp

			for (int i = 0; i < steps; i++) {
				if (i > 0) ippsAdd_64f_I(&times[(i - 1) * nPoints], &times[i * nPoints], nPoints);//integration of (1/tau)
				ippsExp_64f(&times[nPoints * (i)], &exptau[nPoints * (i)], nPoints);
				ippsMulC_64f_I((1E-12) * h, &exptau[nPoints * (i)], nPoints);
				ippsMul_64f_I(&exptau[nPoints * (i)], &vzStorage[i * nPoints], nPoints); //multiply velocities by exp time factor

			}

			for (int i = 0; i < (steps - 1); i++) {
				ippsAdd_64f_I(&vzStorage[i*nPoints], &vzStorage[(i + 1)*nPoints], nPoints); //add all and accumulate in last vector
			}

			ippsMul_64f_I(DOS, &vzStorage[(steps - 1)*nPoints], nPoints);
			ippsMul_64f_I(vz0Storage, &vzStorage[(steps - 1)*nPoints], nPoints);//multiply by initial velocities

			ippsSum_64f(&vzStorage[(steps - 1)*nPoints], nPoints, &total);//sum all elements of velocity vector

			condout[th] = total;
		}
		ippsDivC_64f(condout, condout[0], &tempx[2 * _dataLeng], _dataLeng);//normalize conductivity
		ippsDiv_64f(&tempx[2 * _dataLeng], ones, tempx, _dataLeng);// 1/conductivity to get resistivity 
		ippsSub_64f(tempx, _data, &tempx[_dataLeng], _dataLeng);// (cal[theta]-dat[theta])
		ippsMul_64f_I(&tempx[_dataLeng], &tempx[_dataLeng], _dataLeng);// (cal[theta]-dat[theta])^2
		ippsSum_64f(&tempx[_dataLeng], _dataLeng, &Resdual);//sum( (cal[theta]-dat[theta])^2)
		sum_resdual = sum_resdual + Resdual;
	}



	//cout << sum_resdual << endl;
	return sum_resdual;
}

int calculateAMRO::writefile( double * param)
{
	string filename = "";
	total = 0;
	sum_resdual = 0;
	std::clock_t startT = std::clock();
	updatepara(param, 9);

	Fermi->UpdatePar(param);
	//cout << Fermi->nPoints << endl;
	Fermi->ReturnStart(starts);
	//Fermi->ReturnStart(starts);
	for (int phi_i = 0; phi_i < 4; phi_i++)
	{

		for (int j = 0; j < nPoints; j++) { //initialize Fermi surface to starting grid, only needs to be done once
			output[0 * nPoints + j] = starts[j * 3 + 0];
			output[1 * nPoints + j] = starts[j * 3 + 1];
			output[2 * nPoints + j] = starts[j * 3 + 2];
		}
		ippsCopy_64f(&output[nPoints * (0)], argx, nPoints);//initial velocities for DOS calc;
		ippsCopy_64f(&output[nPoints * (1)], argy, nPoints);
		ippsCopy_64f(&output[nPoints * (2)], argz, nPoints);
		veloX(params, argx, argy, argz, nPoints, tempx, vx); //velocities for DOS are stored in vx, vy, and vz buffers.
		veloY(params, argx, argy, argz, nPoints, tempy, vy);
		veloZ(params, argx, argy, argz, nPoints, tempz, vz);

		ippsSqr_64f_I(vx, nPoints);//in-place square of velocities
		ippsSqr_64f_I(vy, nPoints);
		ippsSqr_64f_I(vz, nPoints);

		ippsAdd_64f(vx, vy, tempx, nPoints);//add all square velocities
		ippsAdd_64f_I(vz, tempx, nPoints);
		ippsSqrt_64f_I(tempx, nPoints);//square root
		ippsDiv_64f(tempx, ones, DOS, nPoints);
		ippsMax_64f(DOS, nPoints, &maxDos);
		ippsMin_64f(DOS, nPoints, &minDos);

		/*argx[0] = 1;
		argy[0] = 2;
		argz[0] = 1;
		veloZ(params, argx, argy, argz, 1, tempz, vz);
		cout << vz[0] << endl;
		*/
		for (int th = 0; th < _dataLeng; th++) {

			//for (int p = 0; p < steps; p++) { //re-initialize times SLOW STEP CREATE TEMP VARIABLE
										  //times[p] = -p;//why -p?
			//	ippsSet_64f(-p, &times[nPoints * p], nPoints);
			//}
			ippsSet_64f(-1, times, steps*nPoints);
			ippsMulC_64f_I(h, times, steps*nPoints);//time stamps

			field[0] = field45 * sin(thetas[th])*cos(phis[phi_i]);  //set field(phi == 0?)
			field[1] = field45 * sin(thetas[th])*sin(phis[phi_i]);
			field[2] = field45 * cos(thetas[th]);
			//for (int z = 0; z < nPoints; z++) {
			//	if (th == 0 ) {
			//		cout << setprecision(20) << times[nPoints*2 + z] << " " << endl;
			//	}
			//}


			for (int i = 1; i < steps; i++) {

				ippsCopy_64f(&output[nPoints * (3 * (i - 1) + 0)], argx, nPoints);//copy arguments for k1;
				ippsCopy_64f(&output[nPoints * (3 * (i - 1) + 1)], argy, nPoints);
				ippsCopy_64f(&output[nPoints * (3 * (i - 1) + 2)], argz, nPoints);
				veloX(params, argx, argy, argz, nPoints, tempx, vx); //calculate velocities;
				veloY(params, argx, argy, argz, nPoints, tempy, vy);
				veloZ(params, argx, argy, argz, nPoints, tempz, vz);
				ippsCopy_64f(vz, &vzStorage[nPoints * (i - 1)], nPoints);//store vz for conductivity later

				//taufun(params, argx, argy, nPoints, tempx, taus, ones);// calculate k dependent tau
				taufundos(params, minDos, maxDos, argx, argy, argz, nPoints, tempx, taus, ones);
				ippsDiv_64f_I(taus, &times[nPoints * (i - 1)], nPoints);
				//ippsDivC_64f_I(tau, &times[nPoints * (i - 1)], nPoints);
				//ippsExp_64f_I(&times[nPoints * (i - 1)], nPoints);
				//ippsMulC_64f_I((1E-12) * h, &times[nPoints * (i - 1)], nPoints);

				//taufun(params, argx, argy, nPoints, tempx, taus);
				//for (int z = 0; z < nPoints; z++) {
				//	if (th == 0 && i == 30) {
				//		cout << setprecision(20) << taus[z] << " " << endl;
				//	}
				//}
				/*for (int z = 0; z < nPoints; z++) {
				if (th == 1 && i == 1) {
				cout << setprecision(20)<< argx[z] << " " << endl;
				}
				}*/


				fx(field, vy, vz, nPoints, tempx, k1x); //calculate evolution in k and store in k1
				fy(field, vx, vz, nPoints, tempy, k1y);
				fz(field, vx, vy, nPoints, tempz, k1z);

				ippsMulC_64f(k1x, h / 2, tempx, nPoints); //prep evolved k step for k2
				ippsMulC_64f(k1y, h / 2, tempy, nPoints);
				ippsMulC_64f(k1z, h / 2, tempz, nPoints);
				ippsAdd_64f(tempx, &output[nPoints * (3 * (i - 1) + 0)], argx, nPoints); //add step to previous k point, load into arguments for k2;
				ippsAdd_64f(tempy, &output[nPoints * (3 * (i - 1) + 1)], argy, nPoints);
				ippsAdd_64f(tempz, &output[nPoints * (3 * (i - 1) + 2)], argz, nPoints);
				veloX(params, argx, argy, argz, nPoints, tempx, vx); //calculate velocities;
				veloY(params, argx, argy, argz, nPoints, tempy, vy);
				veloZ(params, argx, argy, argz, nPoints, tempz, vz);
				fx(field, vy, vz, nPoints, tempx, k2x); //calculate evolution in k and store in k2
				fy(field, vx, vz, nPoints, tempy, k2y);
				fz(field, vx, vy, nPoints, tempz, k2z);

				ippsMulC_64f(k2x, h / 2, tempx, nPoints); //prep evolved k step for k3
				ippsMulC_64f(k2y, h / 2, tempy, nPoints);
				ippsMulC_64f(k2z, h / 2, tempz, nPoints);
				ippsAdd_64f(tempx, &output[nPoints * (3 * (i - 1) + 0)], argx, nPoints); //add step to previous k point, load into arguments for k3;
				ippsAdd_64f(tempy, &output[nPoints * (3 * (i - 1) + 1)], argy, nPoints);
				ippsAdd_64f(tempz, &output[nPoints * (3 * (i - 1) + 2)], argz, nPoints);
				veloX(params, argx, argy, argz, nPoints, tempx, vx); //calculate velocities;
				veloY(params, argx, argy, argz, nPoints, tempy, vy);
				veloZ(params, argx, argy, argz, nPoints, tempz, vz);
				fx(field, vy, vz, nPoints, tempx, k3x); //calculate evolution in k and store in k3
				fy(field, vx, vz, nPoints, tempy, k3y);
				fz(field, vx, vy, nPoints, tempz, k3z);

				ippsMulC_64f(k3x, h, tempx, nPoints); //prep evolved k step for k4
				ippsMulC_64f(k3y, h, tempy, nPoints);
				ippsMulC_64f(k3z, h, tempz, nPoints);
				ippsAdd_64f(tempx, &output[nPoints * (3 * (i - 1) + 0)], argx, nPoints); //add step to previous k point, load into arguments for k4;
				ippsAdd_64f(tempy, &output[nPoints * (3 * (i - 1) + 1)], argy, nPoints);
				ippsAdd_64f(tempz, &output[nPoints * (3 * (i - 1) + 2)], argz, nPoints);
				veloX(params, argx, argy, argz, nPoints, tempx, vx); //calculate velocities;
				veloY(params, argx, argy, argz, nPoints, tempy, vy);
				veloZ(params, argx, argy, argz, nPoints, tempz, vz);
				fx(field, vy, vz, nPoints, tempx, k4x); //calculate evolution in k and store in k4
				fy(field, vx, vz, nPoints, tempy, k4y);
				fz(field, vx, vy, nPoints, tempz, k4z);

				ippsMulC_64f_I(2, k2x, nPoints); //scale k2
				ippsMulC_64f_I(2, k2y, nPoints);
				ippsMulC_64f_I(2, k2z, nPoints);
				ippsMulC_64f_I(2, k3x, nPoints); //scale k3
				ippsMulC_64f_I(2, k3y, nPoints);
				ippsMulC_64f_I(2, k3z, nPoints);

				ippsAdd_64f(k1x, k2x, tempx, nPoints); //add k1 + k2 to temp
				ippsAdd_64f(k1y, k2y, tempy, nPoints);
				ippsAdd_64f(k1z, k2z, tempz, nPoints);

				ippsAdd_64f_I(k3x, tempx, nPoints); //add in k3
				ippsAdd_64f_I(k3y, tempy, nPoints);
				ippsAdd_64f_I(k3z, tempz, nPoints);

				ippsAdd_64f_I(k4x, tempx, nPoints); //add in k4
				ippsAdd_64f_I(k4y, tempy, nPoints);
				ippsAdd_64f_I(k4z, tempz, nPoints);

				ippsMulC_64f_I(h / 6, tempx, nPoints); //scale the entire sum
				ippsMulC_64f_I(h / 6, tempy, nPoints); //scale the entire sum
				ippsMulC_64f_I(h / 6, tempz, nPoints); //scale the entire sum

				ippsAdd_64f(&output[nPoints * (3 * (i - 1) + 0)], tempx, &output[nPoints * (3 * i + 0)], nPoints); //add sum to previous output and store
				ippsAdd_64f(&output[nPoints * (3 * (i - 1) + 1)], tempy, &output[nPoints * (3 * i + 1)], nPoints);
				ippsAdd_64f(&output[nPoints * (3 * (i - 1) + 2)], tempz, &output[nPoints * (3 * i + 2)], nPoints);
			}


			ippsCopy_64f(&output[nPoints * (3 * (steps - 1) + 0)], argx, nPoints);//get velocity for last point
			ippsCopy_64f(&output[nPoints * (3 * (steps - 1) + 1)], argy, nPoints);
			ippsCopy_64f(&output[nPoints * (3 * (steps - 1) + 2)], argz, nPoints);
			veloZ(params, argx, argy, argz, nPoints, tempz, vz);
			ippsCopy_64f(vz, &vzStorage[nPoints * (steps - 1)], nPoints);

			//taufun(params, argx, argy, nPoints, tempx, taus, ones);// calculate k dependent tau for last point
			taufundos(params, minDos, maxDos, argx, argy, argz, nPoints, tempx, taus, ones);
			ippsDiv_64f_I(taus, &times[nPoints * (steps - 1)], nPoints);
			//ippsDivC_64f_I(tau, &times[nPoints * (steps - 1)], nPoints);
			//ippsExp_64f_I(&times[nPoints * (steps - 1)], nPoints);
			//ippsMulC_64f_I((1E-12)*h, &times[nPoints * (steps - 1)], nPoints);
			/*
			ippsCopy_64f(&output[nPoints * (0)], argx, nPoints);//initial velocities for DOS calc;
			ippsCopy_64f(&output[nPoints * (1)], argy, nPoints);
			ippsCopy_64f(&output[nPoints * (2)], argz, nPoints);
			veloX(params, argx, argy, argz, nPoints, tempx, vx); //velocities for DOS are stored in vx, vy, and vz buffers.
			veloY(params, argx, argy, argz, nPoints, tempy, vy);
			veloZ(params, argx, argy, argz, nPoints, tempz, vz);

			ippsSqr_64f_I(vx, nPoints);//in-place square of velocities
			ippsSqr_64f_I(vy, nPoints);
			ippsSqr_64f_I(vz, nPoints);

			ippsAdd_64f(vx, vy, tempx, nPoints);//add all square velocities
			ippsAdd_64f_I(vz, tempx, nPoints);
			ippsSqrt_64f_I(tempx, nPoints);//square root
			ippsDiv_64f(tempx, ones, DOS, nPoints);
			*/
			//need to change!!
			  //ippsDiv_64f_I(taus, times, steps);//exponential stuff, negative tau is taken care of in time
			  //ippsExp_64f_I(times, steps);
			  //ippsMulC_64f_I((1E-12)*h, times, steps);

			ippsCopy_64f(&vzStorage[0], vz0Storage, nPoints);//save initial velocity before exp

			for (int i = 0; i < steps; i++) {
				if (i > 0) ippsAdd_64f_I(&times[(i - 1) * nPoints], &times[i * nPoints], nPoints);//integration of (1/tau)
				ippsExp_64f(&times[nPoints * (i)], &exptau[nPoints * (i)], nPoints);
				ippsMulC_64f_I((1E-12) * h, &exptau[nPoints * (i)], nPoints);
				ippsMul_64f_I(&exptau[nPoints * (i)], &vzStorage[i * nPoints], nPoints); //multiply velocities by exp time factor

			}

			for (int i = 0; i < (steps - 1); i++) {
				ippsAdd_64f_I(&vzStorage[i*nPoints], &vzStorage[(i + 1)*nPoints], nPoints); //add all and accumulate in last vector
			}

			ippsMul_64f_I(DOS, &vzStorage[(steps - 1)*nPoints], nPoints);
			ippsMul_64f_I(vz0Storage, &vzStorage[(steps - 1)*nPoints], nPoints);//multiply by initial velocities

			ippsSum_64f(&vzStorage[(steps - 1)*nPoints], nPoints, &total);//sum all elements of velocity vector

			condout[th] = total;
		}
		ippsDivC_64f(condout, condout[0], &tempdif[2 * _dataLeng], _dataLeng);//normalize conductivity
		ippsDiv_64f(&tempdif[2 * _dataLeng], ones, tempdif, _dataLeng);// 1/conductivity to get resistivity 
		ippsSub_64f(tempdif, &_data[phi_i*_dataLeng], &tempdif[_dataLeng], _dataLeng);// (cal[theta]-dat[theta])
		ippsMul_64f_I(&tempdif[_dataLeng], &tempdif[_dataLeng], _dataLeng);// (cal[theta]-dat[theta])^2
		ippsSum_64f(&tempdif[_dataLeng], _dataLeng, &Resdual);//sum( (cal[theta]-dat[theta])^2)
		sum_resdual = sum_resdual + Resdual;

		filename = "";
		filename.append("conductivity ");
		filename.append(to_string(phis[phi_i]));
		filename.append(".dat");
		ofstream fout;
		fout.open(filename);
		fout.precision(20);

		for (int i = 0; i < _dataLeng; i++) {

			fout << thetas[i] << '\t' << tempdif[i] << endl;
			//fout <<  condout[i] << endl;
			//cout << thetas[i] << "\t" << condout[i] << endl;
		}

		fout.close();
	}
	//duration = (std::clock() - startT) / (Ipp64f)CLOCKS_PER_SEC;
	return 0;
	
	
}


int calculateAMRO::veloX(double *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out)
{
	ippsMulC_64f(kx, 3.74767, temp, length); //term for sin(kx), param3
	vdSin(length, temp, &temp[1 * length]);
	ippsMulC_64f(&temp[1 * length], params[3 - 1] * 11.4215, out, length);

	ippsMulC_64f(kx, 3.74767, temp, length); //term for sin(kx)cos(ky), param4
	vdSin(length, temp, &temp[1 * length]);
	ippsMulC_64f(ky, 3.74767, temp, length);
	vdCos(length, temp, &temp[2 * length]);
	ippsMul_64f_I(&temp[1 * length], &temp[2 * length], length);
	ippsMulC_64f(&temp[2 * length], 22.8429*params[4 - 1], temp, length);
	ippsAdd_64f_I(temp, out, length);

	ippsMulC_64f(kx, 2 * 3.74767, temp, length); //term for sin(2 kx), param5
	vdSin(length, temp, &temp[1 * length]);
	ippsMulC_64f(&temp[1 * length], params[5 - 1] * 11.4215 * 2, temp, length);
	ippsAdd_64f_I(temp, out, length);

	ippsMulC_64f(kx, 3.74767, temp, length); //term for long complicated kz term, param6
	vdSin(length, temp, &temp[1 * length]); // sin kx
	vdCos(length, temp, &temp[2 * length]); // cos kx
	ippsMulC_64f(ky, 3.74767, temp, length);
	vdSin(length, temp, &temp[3 * length]); // sin ky
	vdCos(length, temp, &temp[4 * length]); // cos ky
	ippsMulC_64f(kx, 3.74767 / 2, temp, length);
	vdSin(length, temp, &temp[5 * length]); // sin kx/2
	vdCos(length, temp, &temp[6 * length]); // cos kx/2
	ippsMulC_64f(ky, 3.74767 / 2, temp, length);
	vdSin(length, temp, &temp[7 * length]); // sin ky/2
	vdCos(length, temp, &temp[8 * length]); // cos ky/2
	ippsMulC_64f(kz, 6.6, temp, length);//kz*c/2(c=13.2)
	vdCos(length, temp, &temp[9 * length]); // cos kz/2
	ippsMul_64f(&temp[6 * length], &temp[1 * length], &temp[10 * length], length);//cos kx/2 * sin kx
	ippsMulC_64f_I(65893 * 4, &temp[10 * length], length); // mult by constant
	ippsSub_64f(&temp[4 * length], &temp[2 * length], &temp[11 * length], length);//cos kx - cos ky
	ippsMulC_64f(&temp[11 * length], 65893, &temp[12 * length], length);//mult by constant
	ippsMul_64f(&temp[5 * length], &temp[12 * length], &temp[13 * length], length);// mult by sin kx/2
	ippsAdd_64f(&temp[13 * length], &temp[10 * length], temp, length);//add those two together
	ippsMul_64f_I(&temp[9 * length], temp, length);//mult by cos kz/2
	ippsMul_64f_I(&temp[11 * length], temp, length);//mult by cos kx - cos ky
	ippsMul_64f_I(&temp[8 * length], temp, length);//mult by cos ky/2
	ippsMulC_64f_I(params[6 - 1] * 0.0000866667, temp, length);
	//ippsMulC_64f_I(params[6 - 1], temp, length);
	ippsAdd_64f_I(temp, out, length);

	return 0;
}


int calculateAMRO::veloY(double *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out) {
	ippsMulC_64f(ky, 3.74767, temp, length); //term for sin(ky), param3
	vdSin(length, temp, &temp[1 * length]);
	ippsMulC_64f(&temp[1 * length], params[3 - 1] * 11.4215, out, length);

	ippsMulC_64f(kx, 3.74767, temp, length); //term for cos(kx)sin(ky), param4
	vdCos(length, temp, &temp[1 * length]);
	ippsMulC_64f(ky, 3.74767, temp, length);
	vdSin(length, temp, &temp[2 * length]);
	ippsMul_64f_I(&temp[1 * length], &temp[2 * length], length);
	ippsMulC_64f(&temp[2 * length], 22.8429*params[4 - 1], temp, length);
	ippsAdd_64f_I(temp, out, length);

	ippsMulC_64f(ky, 2 * 3.74767, temp, length); //term for sin(2 ky), param5
	vdSin(length, temp, &temp[1 * length]);
	ippsMulC_64f(&temp[1 * length], params[5 - 1] * 11.4215 * 2, temp, length);
	ippsAdd_64f_I(temp, out, length);

	ippsMulC_64f(kx, 3.74767, temp, length); //term for long complicated kz term, param6
	vdSin(length, temp, &temp[1 * length]); // sin kx
	vdCos(length, temp, &temp[2 * length]); // cos kx
	ippsMulC_64f(ky, 3.74767, temp, length);
	vdSin(length, temp, &temp[3 * length]); // sin ky
	vdCos(length, temp, &temp[4 * length]); // cos ky
	ippsMulC_64f(kx, 3.74767 / 2, temp, length);
	vdSin(length, temp, &temp[5 * length]); // sin kx/2
	vdCos(length, temp, &temp[6 * length]); // cos kx/2
	ippsMulC_64f(ky, 3.74767 / 2, temp, length);
	vdSin(length, temp, &temp[7 * length]); // sin ky/2
	vdCos(length, temp, &temp[8 * length]); // cos ky/2
	ippsMulC_64f(kz, 6.6, temp, length);//kz*c/2(c=13.2)
	vdCos(length, temp, &temp[9 * length]); // cos kz/2
	ippsMul_64f(&temp[8 * length], &temp[3 * length], &temp[10 * length], length);//cos ky/2 * sin ky
	ippsMulC_64f_I(65893 * 4, &temp[10 * length], length); // mult by constant
	ippsSub_64f(&temp[4 * length], &temp[2 * length], &temp[11 * length], length);//cos kx - cos ky CHECK SIGN
	ippsMulC_64f(&temp[11 * length], 65893, &temp[12 * length], length);//mult by constant
	ippsMul_64f(&temp[7 * length], &temp[12 * length], &temp[13 * length], length);// mult by sin kx/2
	ippsSub_64f(&temp[10 * length], &temp[13 * length], temp, length);//add those two together (neg sign)	
	ippsMul_64f_I(&temp[9 * length], temp, length);//mult by cos kz/2
	ippsMul_64f_I(&temp[11 * length], temp, length);//mult by cos kx - cos ky
	ippsMul_64f_I(&temp[6 * length], temp, length);//mult by cos kx/2	
	ippsMulC_64f_I(params[6 - 1] * 0.0000866667, temp, length);
	//ippsMulC_64f_I(params[6 - 1], temp, length);
	ippsAdd_64f_I(temp, out, length);
	return 0;
}

int calculateAMRO::veloZ(double *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out) {
	ippsMulC_64f(kx, 3.74767, temp, length); //term for long complicated kz term, param6
	vdCos(length, temp, &temp[2 * length]); // cos kx
	ippsMulC_64f(ky, 3.74767, temp, length);
	vdCos(length, temp, &temp[4 * length]); // cos ky
	ippsMulC_64f(kx, 3.74767 / 2, temp, length);
	vdCos(length, temp, &temp[6 * length]); // cos kx/2
	ippsMulC_64f(ky, 3.74767 / 2, temp, length);
	vdCos(length, temp, &temp[8 * length]); // cos ky/2
	ippsMulC_64f(kz, 6.6, temp, length);
	vdSin(length, temp, &temp[9 * length]); // sin kzc/2
	ippsMulC_64f(kz, 6.6, temp, length); // changed to kz c/2
	vdSin(length, temp, &temp[13 * length]);// sin kzc **fixed to kz c/2

	ippsSub_64f(&temp[4 * length], &temp[2 * length], &temp[11 * length], length);//cos kx - cos ky
	ippsSqr_64f_I(&temp[11 * length], length);//square it
	ippsMul_64f_I(&temp[9 * length], &temp[11 * length], length);// times sin kz/2
	ippsMul_64f_I(&temp[8 * length], &temp[11 * length], length);// times cos ky/2
	ippsMul_64f_I(&temp[6 * length], &temp[11 * length], length);// times cos ky/2
	ippsMulC_64f(&temp[11 * length], params[6 - 1] * 10.0571*2, out, length);
	ippsMulC_64f(&temp[13 * length], params[7 - 1] * 10.0571*2, &temp[12 * length], length);//h7 term ***removed factor of two when switching to cos kz from cos kz/2
	ippsAdd_64f_I(&temp[12 * length], out, length);


	return 0;
}

int calculateAMRO::fx(Ipp64f * field, Ipp64f *vy, Ipp64f *vz, int length, Ipp64f *temp, Ipp64f *out) {
	/*return -1 / (11538.5) * (vy * field[2] - vz*field[1]);*/
	ippsMulC_64f(vy, field[2], temp, length);
	ippsMulC_64f(vz, field[1], &temp[length], length);
	ippsMulC_64f_I(-1, &temp[length], length);
	ippsAdd_64f_I(&temp[length], temp, length);
	ippsMulC_64f(temp, 1 / (11538.5), out, length);//+1 to run back in time

	return 0;

}


int calculateAMRO::fy(Ipp64f * field, Ipp64f *vx, Ipp64f *vz, int length, Ipp64f *temp, Ipp64f *out) {
	//return -1 / (11538.5) * (vz * field[0] - vx*field[2]);
	ippsMulC_64f(vz, field[0], temp, length);
	ippsMulC_64f(vx, field[2], &temp[length], length);
	ippsMulC_64f_I(-1, &temp[length], length);
	ippsAdd_64f_I(&temp[length], temp, length);
	ippsMulC_64f(temp, 1 / (11538.5), out, length);//+1 to run back in time

	return 0;
}

int calculateAMRO::fz(Ipp64f * field, Ipp64f *vx, Ipp64f *vy, int length, Ipp64f *temp, Ipp64f *out) {
	//return -1 / (11538.5) * (vx * field[1] - vy*field[0]);
	ippsMulC_64f(vx, field[1], temp, length);
	ippsMulC_64f(vy, field[0], &temp[length], length);
	ippsMulC_64f_I(-1, &temp[length], length);
	ippsAdd_64f_I(&temp[length], temp, length);
	ippsMulC_64f(temp, 1 / (11538.5), out, length);//+1 to run back in time

	return 0;
}


int calculateAMRO::taufun(double *params, Ipp64f *kx, Ipp64f *ky, int length, Ipp64f *temp, Ipp64f *out, Ipp64f *ones) {
	ippsDiv_64f(kx, ky, temp, length);
	//cout << temp[0] << endl;
	ippsAtan_64f_A50(temp, &temp[length], length);
	vdSin(length, &temp[length], &temp[2 * length]);//sin(arctan(ky/kx))
	ippsSqr_64f_I(&temp[2 * length], length);
	vdCos(length, &temp[length], &temp[3 * length]);//cos(arctan(ky/kx))
	ippsSqr_64f_I(&temp[3 * length], length);
	ippsSub_64f(&temp[2 * length], &temp[3 * length], &temp[4 * length], length);//sin(arctan(ky/kx))^2-cos(arctan(ky/kx))^2
	//cout << temp[2 * length] << endl;
	ippsPowx_64f_A50(&temp[4 * length], params[9 - 1], &temp[5 * length], length);
	//cout << temp[5 * length] << endl;
	ippsMulC_64f_I(params[8 - 1], &temp[5 * length], length);
	//cout << temp[5 * length] << endl;
	ippsAddC_64f_I(1 / params[1 - 1], &temp[5 * length], length);
	//cout << temp[5 * length] << endl;

	//cout << temp[6 * length] << endl;


	ippsDiv_64f(&temp[5 * length], ones, out, length);
	//cout << temp[6 * length] << endl;





	return 0;
}

/*
int  calculateAMRO::taufun(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, int length, Ipp64f *temp, Ipp64f *out, Ipp64f *ones) {
	ippsDiv_64f(kx, ky, temp, length);
	//cout << temp[0] << endl;
	ippsAtan_64f_A50(temp, &temp[length], length);
	vdSin(length, &temp[length], &temp[2 * length]);//sin(arctan(ky/kx))
	ippsSqr_64f(&temp[2 * length], &temp[3 * length], length);
	vdCos(length, &temp[length], &temp[4 * length]);//cos(arctan(ky/kx))
	ippsSqr_64f(&temp[4 * length], &temp[5 * length], length);
	ippsSub_64f(&temp[3 * length], &temp[5 * length], &temp[6 * length], length);//sin(arctan(ky/kx))^2-cos(arctan(ky/kx))^2
	ippsSqr_64f_I(&temp[6 * length], length);//(sin(arctan(ky/kx))^2-cos(arctan(ky/kx))^2)
	//cout << temp[2 * length] << endl;
	ippsMul_64f(&temp[2 * length], &temp[4 * length], &temp[7 * length], length);
	ippsMulC_64f_I(2, &temp[7 * length], length);//2*sin(arctan(ky / kx))cos(arctan(ky/kx))
	ippsSqr_64f_I(&temp[7 * length], length);//4*sin(arctan(ky / kx))^2cos(arctan(ky/kx))^2
	ippsSub_64f_I(&temp[7 * length], &temp[6 * length], length);// sin(arctan(ky / kx)) ^ 2 - cos(arctan(ky / kx)) ^ 2-4*sin(arctan(ky / kx))^2cos(arctan(ky/kx))^2
	ippsMulC_64f_I(-1 * params[7], &temp[6 * length], length);
	ippsAddC_64f(&temp[6 * length], params[0], out, length);
	//cout << temp[6 * length] << endl;





	return 0;
}
*/
int calculateAMRO::taufundos(Ipp64f *params, Ipp64f minDos, Ipp64f maxDos, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out, Ipp64f *ones) {
	veloX(params, kx, ky, kz, length, &temp[3 * length], temp); //velocities for DOS are stored in vx, vy, and vz buffers.
	veloY(params, kx, ky, kz, length, &temp[3 * length], &temp[length]);
	veloZ(params, kx, ky, kz, length, &temp[3 * length], &temp[2 * length]);

	ippsSqr_64f_I(temp, length);//in-place square of velocities
	ippsSqr_64f_I(&temp[length], length);
	ippsSqr_64f_I(&temp[2 * length], length);

	ippsAdd_64f(temp, &temp[length], &temp[3 * length], length);//add all square velocities
	ippsAdd_64f_I(&temp[2 * length], &temp[3 * length], length);
	ippsSqrt_64f_I(&temp[3 * length], length);//square root
	//ippsMulC_64f_I(1 / temp[3 * length], &temp[3 * length], length);//density of state
	//ippsMulC_64f_I(params[8 - 1], &temp[3 * length], length);
	//ippsAddC_64f(&temp[3 * length], params[1 - 1], out, length);
	ippsDiv_64f(&temp[3 * length], ones, &temp[4 * length], length);//density of state
	ippsAddC_64f_I(-maxDos, &temp[4 * length], length);

	ippsMulC_64f_I((1 / params[8 - 1] - 1 / params[1 - 1]) / (maxDos - minDos), &temp[4 * length], length);
	//	ippsMulC_64f_I(1 / params[1 - 1], &temp[4 * length], length);
	ippsAddC_64f_I(1 / params[8 - 1], &temp[4 * length], length);
	ippsDiv_64f(&temp[4 * length], ones, out, length);
	//cout << *minDos << endl;
	//delete minDos;






	return 0;
}
int calculateAMRO::printPar()
{
	cout << "parameters:";
	for (int i = 0; i < 9; ++i) {
		cout << params[i]<<"  " ;

	}
	cout << endl;
	return 0;
}
