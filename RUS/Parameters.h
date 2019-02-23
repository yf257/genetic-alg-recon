#pragma once
namespace Parameters
{
	__declspec( align( 128) ) struct fitParameters
	{
		double h1; //tau
		double h2; //fermi energy
		double h3; //txpy
		double h4; //txy
		double h5; //t2xp2y
		double h6; //tz
		double h7; //tz2
		double h8; //phi_dep_tau
		double h9; //pow_tau
		double area;//dopping
		double chiSq; //
		
		
	};

    struct arrayBounds
	{
		double start;
		double end;
		HANDLE handle;
		double time;
		int threadID;
	};
}