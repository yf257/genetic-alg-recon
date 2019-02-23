// amro.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
using namespace std;

//static const double density = 8313; // density of the material in (Kg?) grams/meter^3. All units are SI


	/*		QueryPerformanceCounter(&time1);*/
	
	//	QueryPerformanceCounter(&time2);
	//std::cout<<"Time to memcpy: "<<1000*(double)(time2.QuadPart-time1.QuadPart)/(freq.QuadPart)<<"ms"<<std::endl<<std::endl;

int _tmain(int argc, _TCHAR* argv[]) //main function
{
	_CrtMemState s1,s2,s3;
	LARGE_INTEGER gtime1,gtime2,gfreq;  // stores times and CPU frequency for profiling
	QueryPerformanceFrequency(&gfreq);
	//double **** ctens = initElasticConstants(); // 4 dimensional elastic constant array. can probably be simpler (obviously)



	
		VSLStreamStatePtr stream;
		SYSTEMTIME t;
	    GetLocalTime(&t);
	    vslNewStream( & stream, VSL_BRNG_SFMT19937, t.wMilliseconds );
		DataExtractor extractor_theta("thetas.dat");
		DataExtractor extractor_phi("phis.dat");
		double * thetas = extractor_theta.getDataArray();
		double * phi = extractor_phi.getDataArray();
		DataExtractor extractor("data.dat");
		double * data = extractor.getDataArray();
		int nPoints = extractor.getNumberOfLines();
		
		int gridN, cdev;
		double scale, cross;
		cout << "number of gird points in z direction? ";
		cin >> cdev;
		cout  << endl;
		 
		cout << "number of actual inplane grid points";
		cin >> gridN;
		cout << endl; // output that dimension to the user
		
		cout << "Scale factor? ";
		cin >> scale;
		cout  << endl;

		cout << "Crossing Probability? ";
		cin >> cross;
		cout << endl;
		
	
		GeneticAlgorithm geneticAlgorithm(data, nPoints,16, scale, cross, gridN, cdev,thetas,phi);
		
		geneticAlgorithm.calculateMinimum();
		
		geneticAlgorithm.printMinimumParameters();	
	
	while(true){ // bad programming...

		int nGens;
		cout<<"Number of generations: ";
		cin>>nGens;
		cout<<endl;

			//CURRENTLY ONE THREAD CHECK
		//_CrtMemCheckpoint( &s1 );

		QueryPerformanceCounter(&gtime1);
	
		geneticAlgorithm.calculateNewGenerations(nGens);	

		QueryPerformanceCounter(&gtime2);
		std::cout<<"Total time per generation: "<<(1000*(double)(gtime2.QuadPart-gtime1.QuadPart)/(gfreq.QuadPart))/nGens<<"ms"<<std::endl<<std::endl;

	/*	_CrtMemCheckpoint( &s2 );
		if ( _CrtMemDifference( &s3, &s1, &s2 ) )
      _CrtMemDumpStatistics( &s3 );

	cout<<s3.lTotalCount<<endl;*/


		geneticAlgorithm.calculateMinimum();	

		geneticAlgorithm.printMinimumParameters();	
		
		
	

		}
	return 0;
}
