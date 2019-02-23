#pragma once
#include "Parameters.h"


	
#define		timing1 QueryPerformanceCounter(&stime1)
#define		timing2 QueryPerformanceCounter(&stime2)
#define		timingOut	std::cout<<"Time is: "<<1000*(double)(stime2.QuadPart-stime1.QuadPart)/(freq.QuadPart)<<"ms"<<std::endl<<std::endl

class GeneticAlgorithm
{
public:
	GeneticAlgorithm(double* dataSet, int dataSetLength,  int nPopulation, double scaleFactor, double crossingProbability, int gridN,int cdev,Ipp64f* thetas, Ipp64f *phi);
	~GeneticAlgorithm(void);
	
	void calculateMinimum();
	void printMinimumParameters();
	void exportChiSq();
	void calculateNewGenerations(int nGenerations);
	bool validparameterQ(double * parameters);
	
private:
	int gridNum;
	int cdevNum;
	Ipp64f func(double *params, Ipp64f kx);
	int func_cal(double *params, Ipp64f * argkz, Ipp64f * argCos, Ipp64f * argSin, Ipp64f *r, int length, Ipp64f *temp, Ipp64f *out);
	Ipp64f return_area(double *params,double kzval);
	double totalTime;
	_CrtMemState s1,s2,s3;
	LARGE_INTEGER stime1,stime2,freq;  // stores times and CPU frequency for profiling

	static const int nVars = 9; // number of variables, F, dF, etc...
	static const int nParams = 11;  //total parameters including chiSq

	static const int nThreads = 4;

    VSLStreamStatePtr stream;
	int * ints1;
	int * ints2;
	int * ints3;
	int * shuffleIndex;




	int _nPopulation,_dataSetLength;
	double _scaleFactor, _crossingProbability;
	double * _dataSet;

	double * _emat;


	
	double * _gradientCalcs;

	double ** _residualArray;
	double ** _paramArray;

	Ipp64f * _thetas;
	Ipp64f * _phis;


	

	


	Parameters::fitParameters * _populationParametersOld, *_populationParametersNew;
	Parameters::fitParameters _minimumParameters;

	void initializeRandomNumberGenerators();
	void initializeParameters(double* dataSet, int dataSetLength, int nPopulation, double scaleFactor, double crossingProbability,Ipp64f *thetas,Ipp64f *phi);
	double randomDouble(double min, double max);
	double calculateResidual(Parameters::fitParameters * parameters,int threadID);
	double calculateArea(Parameters::fitParameters * parameters, int threadID);
	
	void resetParameters(int nPopulation, double scaleFactor, double crossingProbability);
	static	UINT startResidualThread(LPVOID param);
	void residualCalculatingThread(Parameters::arrayBounds * arrayBounds);
	
	void EvolveParameter(double * pOld, double * pNew, double * rand1, double * rand2, double * rand3);
	
	


	struct threadContents{
		Parameters::arrayBounds arrayBounds;
		GeneticAlgorithm* pThis;
	};

};
