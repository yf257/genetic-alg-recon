#include "stdafx.h"
#include "GeneticAlgorithm.h"




GeneticAlgorithm::GeneticAlgorithm(double* dataSet, int dataSetLength,  int nPopulation, double scaleFactor, double crossingProbability, int gridN, int cdev,Ipp64f* thetas,Ipp64f* phi){
	
	totalTime = 0;  // for timing things
	cdevNum = cdev;
	gridNum = gridN;
	QueryPerformanceFrequency(&freq);

	_nPopulation = nPopulation;
	initializeRandomNumberGenerators();
	
	


	initializeParameters(dataSet, dataSetLength, nPopulation, scaleFactor, crossingProbability,thetas,phi);	

	
	
}

GeneticAlgorithm::~GeneticAlgorithm(void){
	delete [] _populationParametersNew;
	delete [] _populationParametersOld;
}

void GeneticAlgorithm::calculateMinimum(){
	for(int i  = 0; i < _nPopulation; i++){
		if(_populationParametersOld[i].chiSq < _minimumParameters.chiSq){
			_minimumParameters.h1 = _populationParametersOld[i].h1;
			_minimumParameters.h2 = _populationParametersOld[i].h2;
			_minimumParameters.h3 = _populationParametersOld[i].h3;
			_minimumParameters.h4 = _populationParametersOld[i].h4;
			_minimumParameters.h5 = _populationParametersOld[i].h5;
			_minimumParameters.h6 = _populationParametersOld[i].h6;
			_minimumParameters.h7 = _populationParametersOld[i].h7;
			_minimumParameters.h8 = _populationParametersOld[i].h8;
			_minimumParameters.h9 = _populationParametersOld[i].h9;			
			_minimumParameters.area = _populationParametersOld[i].area;
			_minimumParameters.chiSq = _populationParametersOld[i].chiSq;
			
		}
	}
}



double GeneticAlgorithm::calculateResidual(Parameters::fitParameters * parameters, int threadID){
	
	
	double residual  = 0;	
	double * paramPointer = &(parameters->h1);	
	Ipp64f area = 0;
	for(int i = 0; i < nParams; i++){
		_paramArray[threadID][i] = *paramPointer;
	
		paramPointer++;
	}
	if (validparameterQ(_paramArray[threadID])) {
		
		//residual = calculateAMRO1.returnvalue(_paramArray[threadID]);
		//std::cout<< "valid P" << std::endl;
		
		
		area = _paramArray[threadID][9];
		if (area > 0.2 && area < 0.35) {
			calculateAMRO calculateAMRO1(_dataSet, _paramArray[threadID], _thetas, cdevNum, gridNum, _dataSetLength, _phis);
			residual = calculateAMRO1.returnvalue(_paramArray[threadID]);
			//std::cout << std::left << std::setfill(' ') << std::setw(10) << _paramArray[threadID][0] << "," << _paramArray[threadID][1] << "," << _paramArray[threadID][2] << "," << _paramArray[threadID][3] << "," << _paramArray[threadID][4] << "," << _paramArray[threadID][5] << "," << _paramArray[threadID][6] << "," << _paramArray[threadID][7] << "," << _paramArray[threadID][8] << std::endl;
			//std::cout << "valid dopping: " << area << std::setfill(' ')<<"residual: "<< residual <<std::endl;
			//std::cout << "valid P & doping" << std::endl;
		}
		else {
			residual = std::numeric_limits<double>::infinity();
		//	std::cout << "invalid dopping: " << area <<std::endl;
		}
		//std::cout << "dopping: " << area <<std::endl;
	}
	else {
		residual = std::numeric_limits<double>::infinity();
	//	std::cout << "invalid P" << std::endl;
		//std::cout << std::left << std::setfill(' ') << std::setw(10) << _paramArray[threadID][0] << "," << _paramArray[threadID][1] << "," << _paramArray[threadID][2] << "," << _paramArray[threadID][3] << "," << _paramArray[threadID][4] << "," << _paramArray[threadID][5] << "," << _paramArray[threadID][6] << "," << _paramArray[threadID][7] << "," << _paramArray[threadID][8] << std::endl;
	}
	//double * frequencies = calculateFrequencies(_paramArray[threadID]);//<-directly change this one to conductivity 
		//need to look at this part ****************************************************************

	

	
	
	
	return residual;
}


double GeneticAlgorithm::calculateArea(Parameters::fitParameters * parameters, int threadID) {


	//double residual = 0;
	double * paramPointer = &(parameters->h1);
	Ipp64f area = 0;
	for (int i = 0; i < nParams; i++) {
		_paramArray[threadID][i] = *paramPointer;

		paramPointer++;
	}
	//if (validparameterQ(_paramArray[threadID])) {
		
		//residual = calculateAMRO1.returnvalue(_paramArray[threadID]);
		//std::cout<< "valid P" << std::endl;


		area = 0.5*(return_area(_paramArray[threadID],0)+ return_area(_paramArray[threadID], 0.237999));
		//if (area > 0.2 && area < 0.5) {

			//calculateAMRO calculateAMRO1(_dataSet, _paramArray[threadID], _thetas, cdevNum, gridNum, _dataSetLength, _phis);
			//residual = calculateAMRO1.returnvalue(_paramArray[threadID]);
			//std::cout << std::left << std::setfill(' ') << std::setw(10) << _paramArray[threadID][0] << "," << _paramArray[threadID][1] << "," << _paramArray[threadID][2] << "," << _paramArray[threadID][3] << "," << _paramArray[threadID][4] << "," << _paramArray[threadID][5] << "," << _paramArray[threadID][6] << "," << _paramArray[threadID][7] << "," << _paramArray[threadID][8] << std::endl;
			//std::cout << "valid dopping: " << area << std::setfill(' ')<<"residual: "<< residual <<std::endl;
			//std::cout << "valid P & doping" << std::endl;
		//}
		//else {
			//residual = std::numeric_limits<double>::infinity();
			//	std::cout << "invalid dopping: " << area <<std::endl;
		//}
	//	//std::cout << "dopping: " << area <<std::endl;
	//}
	//else {
	//	residual = std::numeric_limits<double>::infinity();
		//	std::cout << "invalid P" << std::endl;
			//std::cout << std::left << std::setfill(' ') << std::setw(10) << _paramArray[threadID][0] << "," << _paramArray[threadID][1] << "," << _paramArray[threadID][2] << "," << _paramArray[threadID][3] << "," << _paramArray[threadID][4] << "," << _paramArray[threadID][5] << "," << _paramArray[threadID][6] << "," << _paramArray[threadID][7] << "," << _paramArray[threadID][8] << std::endl;
	//}
	//double * frequencies = calculateFrequencies(_paramArray[threadID]);//<-directly change this one to conductivity 
		//need to look at this part ****************************************************************






	return area;
}




void GeneticAlgorithm::calculateNewGenerations(int nGenerations){
	double averageTime = 0;

	for(int i = 0; i < nGenerations; i++){
	 std::cout << "Generation" << i <<std::endl;

	 viRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, _nPopulation, ints1, 0, _nPopulation );
	 viRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, _nPopulation, ints2, 0, _nPopulation );
	 viRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, _nPopulation, ints3, 0, _nPopulation );
	 viRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, _nPopulation, shuffleIndex, 0, _nPopulation );
	
	
	 //Parameters::fitParameters  temp;
	 //for(int r = 0; r < _nPopulation; r++){
		//temp = _populationParametersOld[r];
		//_populationParametersOld[r]  = _populationParametersOld[shuffleIndex[r]];
		//_populationParametersOld[shuffleIndex[r]] = temp;
	 //}
	
	 				
	  for(int j = 0; j < _nPopulation; j++){

		
			double * pointerToOldVariable = &_populationParametersOld[j].h1;
			double * pointerToNewVariable = &_populationParametersNew[j].h1;
			double * pointerTog1Variable = &_populationParametersOld[ints1[j]].h1;
			double * pointerTog2Variable = &_populationParametersOld[ints2[j]].h1;
			double * pointerTog3Variable = &_populationParametersOld[ints3[j]].h1;
				
			EvolveParameter(pointerToOldVariable, pointerToNewVariable, pointerTog1Variable, pointerTog2Variable, pointerTog3Variable);
			//_populationParametersNew[j].area = return_area[_populationParametersNew[j]];
				//could be optimzed for vector arithmetic
			/*for(int k = 0; k < nVars; k++){
				double p = randomDouble(0,1);
				if(p > _crossingProbability){
					*pointerToNewVariable = *pointerToOldVariable;
				}
				else{
					*pointerToNewVariable = *pointerTog1Variable + _scaleFactor*(*pointerTog2Variable - *pointerTog3Variable);
				}
				pointerToNewVariable++;
				pointerToOldVariable++;
				pointerTog1Variable++;
				pointerTog2Variable++;
				pointerTog3Variable++;
			}	*/	
	   }	

		totalTime = 0;

		HANDLE threadEvents[nThreads];
		Parameters::arrayBounds threadBounds[nThreads];
		threadContents threadContents[nThreads];

	
		

		for(int m = 0; m<nThreads; m++){		

			 threadEvents[m] = CreateEvent(NULL, FALSE, FALSE, NULL);
			int nPopulationPerThread =  _nPopulation/nThreads;
			if(m != (nThreads-1)){
				threadBounds[m].start = m*nPopulationPerThread;
				threadBounds[m].end = (m+1)*nPopulationPerThread - 1;
			}
			else{
				threadBounds[m].start = m*nPopulationPerThread;
				threadBounds[m].end = (m+1)*nPopulationPerThread - 1 + _nPopulation%nThreads;
			}		
		
			threadBounds[m].handle = threadEvents[m];
			threadBounds[m].time = 0;
			threadBounds[m].threadID = m;
			threadContents[m].arrayBounds = threadBounds[m];
			threadContents[m].pThis = this;
	
			AfxBeginThread(startResidualThread, (LPVOID) &threadContents[m]);		
			
		}
		
	//	std::cout<<s2.lTotalCount<<std::endl;
		WaitForMultipleObjects(nThreads,threadEvents,TRUE,INFINITE);	
		std::string filename = "";
		filename = "";
		filename.append("generation");
	    filename.append(std::to_string(i+10));
		filename.append(".dat");
		std::ofstream out;
	    out.open(filename);
		out.precision(15);
		for (int genera = 0; genera < _nPopulation; ++genera) {
		
			out << _populationParametersNew[genera].h1 << '\t' << _populationParametersNew[genera].h2 << '\t' << _populationParametersNew[genera].h3 << '\t' << _populationParametersNew[genera].h4 << '\t' << _populationParametersNew[genera].h5 << '\t' << _populationParametersNew[genera].h6 << '\t' << _populationParametersNew[genera].h7 << '\t' << _populationParametersNew[genera].h8 << '\t' << _populationParametersNew[genera].h9 << '\t' << _populationParametersNew[genera].area <<'\t' <<_populationParametersNew[genera].chiSq << std::endl;
		}
		out.close();
		
		for(int timerIndex = 0; timerIndex < nThreads; timerIndex++){
		totalTime += threadContents[timerIndex].arrayBounds.time;
		}
	
		averageTime +=totalTime;
		
		mkl_free_buffers();
	
	/*	calculateMinimum();
		exportChiSq();*/
	
	}
	

	std::cout<<"Average time per generation for a thread: "<<averageTime/(nThreads*nGenerations)<<"ms"<<std::endl<<std::endl;
		
		
}

UINT GeneticAlgorithm::startResidualThread(LPVOID param){
	threadContents * contents = (threadContents*) param;
	contents->pThis->residualCalculatingThread(&(contents->arrayBounds));
	return 0;
}

void GeneticAlgorithm::residualCalculatingThread(Parameters::arrayBounds * arrayBounds){
	
	

	LARGE_INTEGER time1,time2;


	
	for(int i = arrayBounds->start; i <= arrayBounds->end; i++){
	
			QueryPerformanceCounter(&time1);

		_populationParametersNew[i].area = calculateArea(&_populationParametersNew[i],arrayBounds->threadID);
		_populationParametersNew[i].chiSq = calculateResidual(&_populationParametersNew[i], arrayBounds->threadID);
		QueryPerformanceCounter(&time2);
		arrayBounds->time += 1000*(double)(time2.QuadPart-time1.QuadPart)/(freq.QuadPart);
			if(_populationParametersNew[i].chiSq < _populationParametersOld[i].chiSq){
				_populationParametersOld[i] = _populationParametersNew[i];
			}
	}
	SetEvent(arrayBounds->handle);

	return;
}

//double ** GeneticAlgorithm::derivatives(){
//
//
//
//}

void GeneticAlgorithm::exportChiSq(){
    	std::ofstream out;
		out.open("chiSq3.dat",std::ios_base::app);
		out.precision(15);
		out<<_minimumParameters.chiSq<<"\t";
		out.close();
}

void GeneticAlgorithm::printMinimumParameters(){
	
		//double rms = 0;
		double _printmin[9] = { _minimumParameters.h1 ,_minimumParameters.h2 ,_minimumParameters.h3 ,_minimumParameters.h4 ,_minimumParameters.h5 ,_minimumParameters.h6 ,_minimumParameters.h7 ,_minimumParameters.h8 ,_minimumParameters.h9 };


	std::cout<<std::left<<std::setfill(' ')<<std::setw(10)<<"h1"<<std::setw(10)<<"h2"<<std::setw(10)<<"h3"<<std::setw(10)<<"h4"<<std::setw(10)<<"h5"<<std::setw(10)<<"h6"<<std::setw(10)<<"h7"<<std::setw(10)<<"h8"<<std::setw(10)<<"h9"<<std::endl;
	std::cout<<std::left<<std::setfill(' ')<<std::setw(10)<<_minimumParameters.h1<<std::setw(10)<<_minimumParameters.h2<<std::setw(10)<<_minimumParameters.h3<<std::setw(10)<<_minimumParameters.h4<<std::setw(10)<<_minimumParameters.h5<<std::setw(10)<<_minimumParameters.h6<<std::setw(10)<<_minimumParameters.h7<<std::setw(10)<<_minimumParameters.h8<<std::setw(10)<<_minimumParameters.h9<<std::endl<<std::endl<<"Residual: "<< _minimumParameters.chiSq <<" %"<<std::endl<<std::endl;

	//std::cout<<"c11: "<<_minimumParameters.c11<<" "<<"c22: "<<_minimumParameters.c22<<" "<<"c33: "<<_minimumParameters.c33<<std::endl<<"c44: "<<_minimumParameters.c44<<" "<<"c55: "<<_minimumParameters.c55<<" "<<"c66: "<<_minimumParameters.c66<<"c12: "<<_minimumParameters.c12<<" "<<"c13: "<<_minimumParameters.c13<<" "<<"c23: "<<_minimumParameters.c23<<std::endl<<"Residual: "<<_minimumParameters.chiSq<<std::endl<<std::endl;
	
		std::ofstream out;
		out.open("parameter.dat");
		out.precision(15);
		out<<_minimumParameters.h1<<'\t'<<_minimumParameters.h2<<'\t'<<_minimumParameters.h3<<'\t'<<_minimumParameters.h4<<'\t'<<_minimumParameters.h5<<'\t'<<_minimumParameters.h6<<'\t'<<_minimumParameters.h7<<'\t'<<_minimumParameters.h8<<'\t'<<_minimumParameters.h9;
		out.close();
		

		calculateAMRO minimum(_dataSet, _printmin, _thetas, cdevNum, gridNum, _dataSetLength,_phis);
		minimum.writefile(_printmin);
			//std::cout<<std::left<<std::setfill(' ')<< std::setw(10) << rms <<" %" << std::endl << std::endl;
			
			
			
		
	

	

}
Ipp64f GeneticAlgorithm::func(double *params, Ipp64f kx) {
	Ipp64f energy;
	Ipp64f t = 17582.4;
	Ipp64f a = 3.74767;
	energy = params[1]* t -2 *t *params[2]*(cos(kx *a) + 1) -4 *t *params[ 3]*cos(kx *a) -2 *t* params[ 4] * (cos(2 *(kx *a)) + 1) -2 *t *params[ 5] * (cos(kx *a) -1)*(cos(kx *a) - 1) * cos((kx *a) / 2)  - 2 * t*params[6] ;
	

	return energy;
}
bool GeneticAlgorithm::validparameterQ(double * parameters) {
	Ipp64f productenergy;
	Ipp64f a = 3.74767;
	productenergy = func(parameters, 0)* func(parameters, 3.1415926 / a);
	
	if (productenergy >= 0 || parameters[0]<=0) { 
		return false; 
	}
	else { 
		return true;
	}
}
int  GeneticAlgorithm::func_cal(double *params, Ipp64f * argkz, Ipp64f * argCos, Ipp64f * argSin, Ipp64f *r, int length, Ipp64f *temp, Ipp64f *out) {
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

}
Ipp64f  GeneticAlgorithm::return_area(double *params,double kzval) {
	int fineN = 1000;
	Ipp64f *kz = new Ipp64f[fineN];//for inverting
	ippsSet_64f(kzval, kz, fineN);
	Ipp64f *gtheta=new Ipp64f[fineN];
	Ipp64f *r = new Ipp64f[fineN];

	Ipp64f *temp1 = new Ipp64f[fineN*15];;
	Ipp64f *argCos = new Ipp64f[fineN];
	Ipp64f *argSin = new Ipp64f[fineN];


	Ipp64f area;
	
	Ipp64f *funcval = new Ipp64f[fineN];
	Ipp64f *dfunc1 = new Ipp64f[fineN];
	Ipp64f *dfunc2 = new Ipp64f[fineN];
	Ipp64f *dfunc = new Ipp64f[fineN];
	Ipp64f *rtemp = new Ipp64f[fineN];
	Ipp64f *kx = new Ipp64f[fineN];
	Ipp64f *ky = new Ipp64f[fineN];
	for (int j = 0; j < fineN; j++)
	{
		/*cout << starts[2 * i] << " " << starts[2 * i + 1] << endl;*/
		gtheta[ j] = 2 * 3.1415926 / (fineN - 1) * j;
		r[j] = 0.3;
	}
	vdSin(fineN, gtheta, argSin);
	vdCos(fineN, gtheta, argCos);
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
		func_cal(params, kz, argCos, argSin, r, fineN, temp1, funcval);
		//func(nPoints, temp1,funcval);
		//cout << argCos[1] << endl;
		//ippsAbs_64f(funcval, temp1, fineN);

		//cout << i << endl;
		//cout << funcval[0] << endl;
		//cout << subMaxR[0] << endl;

		ippsAddC_64f(r, 1E-6, rtemp, fineN);
		func_cal(params, kz, argCos, argSin, rtemp, fineN, temp1, dfunc1);
		//func(nPoints, temp1,dfunc1);
		ippsAddC_64f(r, -1E-6, rtemp, fineN);
		func_cal(params, kz, argCos, argSin, rtemp, fineN, temp1, dfunc2);
		//func(nPoints, temp1, dfunc2);
		ippsMulC_64f_I(-1, dfunc2, fineN);
		ippsAdd_64f(dfunc1, dfunc2, dfunc, fineN);
		ippsDivC_64f_I(2E-6, dfunc, fineN);

		ippsDiv_64f_I(dfunc, funcval, fineN);
		ippsSub_64f_I(funcval, r, fineN);

	}

	//for (int i = 0; i < nPoints; i++)
	//{
	//	cout << r[i] << endl;
	//}
	//ippsMul_64f(r, argCos, kx, nPoints);
	//ippsMul_64f(r, argSin, ky, nPoints);

	ippsMul_64f(r, argCos, kx, fineN);
	ippsMul_64f(r, argSin, ky, fineN);
	ippsMul_64f(kx, &ky[1], temp1, fineN - 1);
	ippsMul_64f(&kx[1], ky, &temp1[fineN - 1], fineN - 1);
	ippsSub_64f_I(&temp1[fineN - 1], temp1, fineN - 1);
	ippsSum_64f(temp1, fineN - 1, &area);
	
	delete temp1;
	delete argCos;
	delete argSin;
	delete gtheta;
	delete kz;
	delete funcval;
	delete dfunc1;
	delete dfunc2;
	delete dfunc;
	delete r;
	delete rtemp;
	delete kx;
	delete ky;

	return -1 + 2 *(1 - 0.177882 *area);
}
