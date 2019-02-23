#include "stdafx.h"
#include "GeneticAlgorithm.h"

void GeneticAlgorithm::initializeRandomNumberGenerators(){
	SYSTEMTIME t;
	GetLocalTime(&t);
	unsigned int max = _nPopulation - 1;

	vslNewStream( & stream, VSL_BRNG_SFMT19937, t.wMilliseconds );
	ints1 = new int[_nPopulation];
    ints2 = new int[_nPopulation];
    ints3 = new int[_nPopulation];
	shuffleIndex = new int[_nPopulation];
}



void GeneticAlgorithm::initializeParameters(double* dataSet, int dataSetLength, int nPopulation, double scaleFactor, double crossingProbability, Ipp64f* thetas,Ipp64f* phi){

	_scaleFactor = scaleFactor;
	_crossingProbability = crossingProbability;
	_dataSetLength = dataSetLength;
	_dataSet = dataSet;
	_thetas = thetas;
	_phis = phi;
	_residualArray = new double*[nThreads];
	_paramArray = new double*[nThreads];
	Ipp64f area_0 = 0;
	Ipp64f area_half = 0;
	double *paramPointer = new double[9];
	int totlen = _dataSetLength ;
	for(int i = 0; i < nThreads; i++){
	_residualArray[i] = new double[totlen];
	_paramArray[i] = new double[nParams];

	}
	
	_populationParametersOld = (Parameters::fitParameters *) mkl_malloc(sizeof(Parameters::fitParameters)*nPopulation,16);
	_populationParametersNew = (Parameters::fitParameters *) mkl_malloc(sizeof(Parameters::fitParameters)*nPopulation,16);
	
	for (int i = 0; i < _nPopulation; i++) {

		_populationParametersOld[i].h1 = (randomDouble(0.02, 0.18));

		//_populationParametersOld[i].h2 = (randomDouble(350, 750));
		_populationParametersOld[i].h2 = (randomDouble(190*0.8, 190*1.2));
		//_populationParametersOld[i].h2 = 181.98;
		//_populationParametersOld[i].h2 = 557.568650327878;

		//_populationParametersOld[i].h3 = (randomDouble(500, 550));
		_populationParametersOld[i].h3 = 190;

		_populationParametersOld[i].h4 = -1 * (randomDouble(190*0.09, 190*0.3));//(30,100)
		//_populationParametersOld[i].h4 = -30.4;
		_populationParametersOld[i].h5 = (randomDouble(190*0.02, 190*0.13));
		//_populationParametersOld[i].h5 = 13.397;

		_populationParametersOld[i].h6 = (randomDouble(190 * 0.02, 190 * 0.13));
		//_populationParametersOld[i].h6 = 13.397;
		_populationParametersOld[i].h7 = 0;//(randomDouble(0, 3)) - 1.5;
		//_populationParametersOld[i].h7 = -1.2;
		//_populationParametersOld[i].h5 = 16.1129375749168;
		//_populationParametersOld[i].h6 = 2.57290092103488;
		//_populationParametersOld[i].h7 = -1.28989597536799;

		_populationParametersOld[i].h8 = randomDouble(0.001, 0.019);
		//_populationParametersOld[i].h8 = 0;
		//_populationParametersOld[i].h9 = 2 * (ceil(randomDouble(0, 15)));
		_populationParametersOld[i].h9 = 12;
		paramPointer[0] = _populationParametersOld[i].h1;
		paramPointer[1] = _populationParametersOld[i].h2;
		paramPointer[2] = _populationParametersOld[i].h3;
		paramPointer[3] = _populationParametersOld[i].h4;
		paramPointer[4] = _populationParametersOld[i].h5;
		paramPointer[5] = _populationParametersOld[i].h6;
		paramPointer[6] = _populationParametersOld[i].h7;
		paramPointer[7] = _populationParametersOld[i].h8;
		paramPointer[8] = _populationParametersOld[i].h9;
		
		area_0 = return_area(paramPointer,0);//kz=0;
	    area_half = return_area(paramPointer, 0.237999);//kz=Pi/c
		
		if (area_0 > 0.2 && area_0 < 0.35 && validparameterQ(paramPointer)&& area_half > 0.2 && area_half < 0.35) {
			//std::cout << i << "       " << area << std::endl;
			_populationParametersOld[i].area = (area_0+ area_half)/2;
			_populationParametersOld[i].chiSq = calculateResidual(&_populationParametersOld[i], 0);
			std::cout << i << "       " << _populationParametersOld[i].area << std::endl;
			
		}
		else {
			--i;
		}

	}
	std::string filename = "";
	filename = "";
	filename.append("generation");
	filename.append(std::to_string(00));
	filename.append(".dat");
	std::ofstream out;
	out.open(filename);
	out.precision(15);
	for (int genera = 0; genera < _nPopulation; ++genera) {

		out << _populationParametersOld[genera].h1 << '\t' << _populationParametersOld[genera].h2 << '\t' << _populationParametersOld[genera].h3 << '\t' << _populationParametersOld[genera].h4 << '\t' << _populationParametersOld[genera].h5 << '\t' << _populationParametersOld[genera].h6 << '\t' << _populationParametersOld[genera].h7 << '\t' << _populationParametersOld[genera].h8 << '\t' << _populationParametersOld[genera].h9 << '\t' << _populationParametersOld[genera].area << '\t' << _populationParametersOld[genera].chiSq << std::endl;
	}
	out.close();
	delete paramPointer;
		_minimumParameters.h1 = 1;
		_minimumParameters.h2 = 1;
		_minimumParameters.h3 = 1;
		_minimumParameters.h4 = 1;
		_minimumParameters.h5 = 1;
		_minimumParameters.h6 = 1;
		_minimumParameters.h7 = 1;
		_minimumParameters.h8 = 1;
		_minimumParameters.h9 = 1;
		_minimumParameters.area = 1;
		_minimumParameters.chiSq = std::numeric_limits<double>::infinity();

		
}

void GeneticAlgorithm::resetParameters(int nPopulation, double scaleFactor, double crossingProbability){
	_scaleFactor = scaleFactor;
	_crossingProbability = crossingProbability;
	_nPopulation = nPopulation;
	Ipp64f area_0 = 0;
	Ipp64f area_half = 0;
	double * paramPointer=new double[9];
	delete [] _populationParametersOld;
	delete [] _populationParametersNew;

	_populationParametersOld = (Parameters::fitParameters *) mkl_malloc(sizeof(Parameters::fitParameters)*nPopulation,16);
	_populationParametersNew = (Parameters::fitParameters *) mkl_malloc(sizeof(Parameters::fitParameters)*nPopulation,16);

	for(int i  = 0; i < _nPopulation; i++){
	
		_populationParametersOld[i].h1 = (randomDouble(0.03, 0.06));

		//_populationParametersOld[i].h2 = (randomDouble(350, 750));
		//_populationParametersOld[i].h2 = (randomDouble(500, 600));
		_populationParametersOld[i].h2 = 569.868225464555;
		//_populationParametersOld[i].h2 = 557.568650327878;

		//_populationParametersOld[i].h3 = (randomDouble(500, 550));
		_populationParametersOld[i].h3 = 533.6164;

		//_populationParametersOld[i].h4 = -1 * (randomDouble(90, 120));//(30,100)
		_populationParametersOld[i].h4 = -113.561272;
		//_populationParametersOld[i].h5 = (randomDouble(10, 25));
		_populationParametersOld[i].h5 = 23.2192;

		//_populationParametersOld[i].h6 = (randomDouble(5, 11));
		_populationParametersOld[i].h6 = 8.7296719;
		//_populationParametersOld[i].h7 = (randomDouble(0, 3)) - 1.5;
		_populationParametersOld[i].h7 = -0.89335299;
		//_populationParametersOld[i].h5 = 16.1129375749168;
		//_populationParametersOld[i].h6 = 2.57290092103488;
		//_populationParametersOld[i].h7 = -1.28989597536799;

		_populationParametersOld[i].h8 = randomDouble(30, 70);
		//_populationParametersOld[i].h8 = 0;
		//_populationParametersOld[i].h9 = 2 * (ceil(randomDouble(0, 10)));
		_populationParametersOld[i].h9 = 2;
		paramPointer[0] = _populationParametersOld[i].h1;
		paramPointer[1] = _populationParametersOld[i].h2;
		paramPointer[2] = _populationParametersOld[i].h3;
		paramPointer[3] = _populationParametersOld[i].h4;
		paramPointer[4] = _populationParametersOld[i].h5;
		paramPointer[5] = _populationParametersOld[i].h6;
		paramPointer[6] = _populationParametersOld[i].h7;
		paramPointer[7] = _populationParametersOld[i].h8;
		paramPointer[8] = _populationParametersOld[i].h9;

		area_0 = return_area(paramPointer, 0);//kz=0;
		area_half = return_area(paramPointer, 0.237999);//kz=Pi/c
		if (area_0 > 0.2 && area_0 < 0.35 && validparameterQ(paramPointer) && area_half > 0.2 && area_half < 0.5) {
			//std::cout << i << "       " << area << std::endl;
			_populationParametersOld[i].area = (area_0 + area_half) / 2;
			_populationParametersOld[i].chiSq = calculateResidual(&_populationParametersOld[i], 0);
	     	
		}
		else {
			--i;
		}
		
	}

	//delete _integerDistribution;
	delete [] ints1;
	delete [] ints2;
	delete [] ints3;
	delete[] paramPointer;
	ints1 = new int[_nPopulation];
    ints2 = new int[_nPopulation];
    ints3 = new int[_nPopulation];
	
}

double GeneticAlgorithm::randomDouble(double min, double max){				//change this to return full aray: faster
	double randNum;
	 vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 1, &randNum, min, max);
	return randNum;
}


void GeneticAlgorithm::EvolveParameter(double * pOld, double * pNew, double * rand1, double * rand2, double * rand3) {
	double p = randomDouble(0, 1);
	
	if (p > _crossingProbability) {
		pNew[0] = pOld[0];
	}
	else {
		
		pNew[0] = rand1[0] + _scaleFactor * (rand2[0] - rand3[0]);
		
	}

	//p = randomDouble(0, 1);
	if (p > _crossingProbability) {
		pNew[1] = pOld[1];
	}
	else {
		pNew[1] = rand1[1] + _scaleFactor * (rand2[1] - rand3[1]);
	}

	p = randomDouble(0, 1);
	if (p > _crossingProbability) {
		pNew[2] = pOld[2];
	}
	else {
		pNew[2] = rand1[2] + _scaleFactor * (rand2[2] - rand3[2]);
	}

	p = randomDouble(0, 1);
	if (p > _crossingProbability) {
		pNew[3] = pOld[3];
	}
	else {
		pNew[3] = rand1[3] + _scaleFactor * (rand2[3] - rand3[3]);
	}

	p = randomDouble(0, 1);
	if (p > _crossingProbability) {
		pNew[4] = pOld[4];
	}
	else {
		pNew[4] = rand1[4] + _scaleFactor * (rand2[4] - rand3[4]);
	}

	p = randomDouble(0, 1);
	if (p > _crossingProbability) {
		pNew[5] = pOld[5];
	}
	else {
		pNew[5] = rand1[5] + _scaleFactor * (rand2[5] - rand3[5]);
	}

	p = randomDouble(0, 1);
	if (p > _crossingProbability) {
		pNew[6] = pOld[6];
	}
	else {
		pNew[6] = rand1[6] + _scaleFactor * (rand2[6] - rand3[6]);
	}
	
	//pNew[7] = pOld[7];
	p = randomDouble(0, 1);
	if (p > _crossingProbability) {
		pNew[7] = pOld[7];
	}
	else {
		pNew[7] = rand1[7] + _scaleFactor * (rand2[7] - rand3[7]);
	}
	//pNew[8] = pOld[8];
	p = randomDouble(0, 1);
	if (p > _crossingProbability) {
		pNew[8] = pOld[8];
	}
	else {
		//pNew[8] = ceil((rand1[8]/2 + _scaleFactor * (rand2[8] - rand3[8]))/2)*2;
		pNew[8] = pOld[8];
	}
	
}