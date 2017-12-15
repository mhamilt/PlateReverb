//
//  PlateClass.hpp
//  FDTD_C_Plate
//
//  Created by mhamilt on 15/10/2017.
//  Copyright © 2017 mhamilt. All rights reserved.
//
//
#ifndef PlateClass_hpp
#define PlateClass_hpp
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <algorithm>

class FD_Plate
{
public:
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Constructors/Assignments
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	FD_Plate()
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Flags
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		bcFlag_ = 0;		// set boundary condition(s);
		outFlag_ = 1;		// set output type 0: displacement, 1: velocity
		setupFlag_ = false;	// flag if Setup() has been run
		
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Physical Parameters
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		//         E	  nu	Rho
		// Steel : 2e11   0.30	8050
		// Alum  : 7e10   0.35	2700
		// Lead  : 1.6e10 0.44	11340
		// Wood  : 1e10   0.40	480
		
		// // wood
		E_ = 11e9;						// Young's modulus
		rho_ = 480;						// density (kg/m^3)
		nu_ = .5;						// Poisson Ratios (< .5)
		
//		E_ = 2e11;						// Young's modulus
//		rho_ = 8050;						// density (kg/m^3)
//		nu_ = .3;						// Poisson Ratios (< .5)
		
		
		H_ = .006;						// thickness (m)
		Lx_ = .7;							// x-axis plate length (m)
		Ly_ = .7;							// y-axis plate length (m)
		loss_[0] = 100; loss_[1] = 1;
		loss_[2] = 1000; loss_[3] = .2;    // loss [freq.(Hz), T60;...]
		
		// I/O Parameters
		rp_[0] = .45; rp_[1]=.65; rp_[2] = .45; rp_[3]= .15; // readout position as percentage.
		
		//Excitation
		ctr_[0] = .35; ctr_[1] = .45;	// centre point of excitation as percentage
		wid_ = .25;					// width (m)
		u0_ = 0; v0_ = 1;				// excitation displacement and velocity
		
	};
	
	~FD_Plate()
	{
//		delete u_;
//		delete u1_;
//		delete u2_;
//		//// cleanup is equally ugly:
//		for(int i = 0; i < interpOrder_; ++i)
//		{
//			delete[] interpLookTable_[i];
//		}
//		delete[] interpLookTable_;
	};
	
	// Copy Assignment
	//	FD_Plate& operator= (const FD_Plate&);
	
	// Move-Constructor
	//	FD_Plate (FD_Plate&&);
	
	// Move-Assignment
	//	FD_Plate& operator= (FD_Plate&&);
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Methods
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	void setup(double, bool);	//Semi-Finished
	void setLoss(double, double);
	void setCoefs(bool);
	void setGrid();
	void setInitialCondition();
	void setOutput(double, double);
	void setStereoOutput(double, double);
	void setOutType(bool);
	void setInterpOut(double, double);
	void updateScheme();
	void addForce(double);
	void addStrike();
	void processIO();
	void printInfo();
	void printCoefs();
	double getOutput(bool);
	double getInterpOut();
	double reverb(double);
	void getStereoOutput(bool, double&, double&);
	int sgn(double d){if(d<=0){return 0;}else{return 1;}};
	
private:

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Methods
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//	void updateInternal();
	//	void updateSides();
	//	void updateCorners();
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Interpolation
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double** getInterpLookTable();
	
	static constexpr const int interpOrder_ = 4;
	static constexpr const int interpRes_ = 1000;
	int interpZeroIndex_;
	double interpPointY_, interpAlphaY_;
	double interpPointX_, interpAlphaX_;
	double linearAlphas_[4];
	int xAlphaIndex_, yAlphaIndex_;
	int  linearInterpInds_[4];
	int* xInterpIndeces_ = new int[interpOrder_];
	int* yInterpIndeces_ = new int[interpOrder_];
	double** interpLookTable_ = getInterpLookTable();
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Constants
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	static constexpr double pi {3.14159265358979323846};
	static constexpr int max_grid_size {3000};		// real-time limit 3000 points approx.
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Flags
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	bool bcFlag_;		// set boundary condition(s);
	bool outFlag_;		// set output type 0: displacement, 1: velocity
	bool setupFlag_;		// flag if Setup() has been run
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Physical Parameters
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// // wood
	double E_;			// Young's modulus
	double rho_;		// density (kg/m^3)
	
	double H_;			// thickness (m)
	double Lx_;			// x-axis plate length (m)
	double Ly_;			// y-axis plate length (m)
	double loss_[4];	// loss [freq.(Hz), T60;...]
	double nu_;			// Poisson Ratios (< .5)
	
	// I/O Parameters
	double rp_[4]; // readout position as percentage.
	
	//Excitation
	double u0_ , v0_, wid_, ctr_[2]; // excitation displacement and velocity
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Derived Parameters
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double D_, kappa_, hmin_, h_, mu_, k_, SR_, readcoord_x_,readcoord_y_, readoutpos_;
	int Nx_, Ny_, ss_, li_, lo_, lo_l_, lo_r_;
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Allocate Memory
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double * u_,* u1_,* u2_;
	double *dummy_ptr_;
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Loss coefficients
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double sigma0_ ,sigma1_, z1_, z2_;
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Scheme Coefficient
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// coefficients are named based on position on the x and y axes.
	double A00_, B00_, B01_, B11_, B02_, BC1_, BC2_, C00_, C01_, d0_;
	
};


#endif /* PlateClass_hpp */
