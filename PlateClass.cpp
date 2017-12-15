//
//  PlateClass.cpp
//  FDTD_C_Plate
//
//  Created by mhamilt on 15/10/2017.
//  Copyright © 2017 mhamilt. All rights reserved.
//
//  Class file for a FDTD Plate

#include "PlateClass.hpp"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

void FD_Plate::setLoss(double lowt60, double hight60percent)
{
	
	loss_[1] = lowt60; loss_[3] = lowt60*hight60percent;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Loss coefficients
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	z1_ = 2*kappa_*(2*pi*loss_[0])/(2*pow(kappa_,2));
	z2_ = 2*kappa_*(2*pi*loss_[2])/(2*pow(kappa_,2));
	sigma0_ = 6*log(10)*(-z2_/loss_[1] + z1_/loss_[3])/(z1_-z2_);
	sigma1_ = 6*log(10)*(1/loss_[1] - 1/loss_[3])/(z1_-z2_);
	
	if (setupFlag_)
	{
		setGrid();
		setCoefs(bcFlag_);
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

void FD_Plate::setGrid()
{
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Scheme Spacing
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	const int max_x_grid = 30.;
	// stability condition
	if(Lx_/max_x_grid<(sqrt(4*k_*(sigma1_+sqrt(pow(sigma1_,2)+pow(kappa_,2))))))
	{
		hmin_ = (sqrt(4*k_*(sigma1_+sqrt(pow(sigma1_,2)+pow(kappa_,2)))));
	}
	else
	{
		hmin_ =	Lx_/max_x_grid;
	}
	Nx_ = floor(Lx_/hmin_);		// number of segments x-axis
	Ny_ = floor(Ly_/hmin_);		// number of segments y-axis
	
	h_ = sqrt(Lx_*Ly_/(Nx_*Ny_));;	// adjusted grid spacing x/y
	Nx_ = Nx_+1; Ny_ = Ny_+1;		// grid point number x and y
	mu_ = (kappa_ * k_)/pow(h_,2);	// scheme parameter
	ss_ = Nx_*Ny_;					// total grid size.
	
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

void FD_Plate::setCoefs(bool bcType = 0)
{
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Scheme Coefficients
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	//update flag
	bcFlag_ = bcType;
	
	// coefficients are named based on position on the x and y axes.
	A00_ = 1/(1+k_*sigma0_); // Central Loss Coeffient (INVERTED)
	
	//// Current time step (B) coeffients
	// There are six unique coefficients for B coefs
	B00_ = (-pow(mu_,2)*20 + (2*sigma1_*k_/pow(h_,2))*-4 + 2) * A00_; // center
	B01_ = (-pow(mu_,2)*-8 + (2*sigma1_*k_/pow(h_,2))) * A00_;		  // 1-off
	B11_ = (-pow(mu_,2)*2) * A00_;									  // diag
	B02_ = (-pow(mu_,2)*1) * A00_;									  // 2-off
	
	if(bcType){ // Clamped Boundary Coefficients
		BC1_ = (-pow(mu_,2)*21 + (2*sigma1_*k_/pow(h_,2))*-4 + 2) * A00_; // Side
		BC2_ = (-pow(mu_,2)*22 + (2*sigma1_*k_/pow(h_,2))*-4 + 2) * A00_; // Corner
	}
	else { // Simply Supported Boundary Coefficients
		BC1_ = (-pow(mu_,2)*19 + (2*sigma1_*k_/pow(h_,2))*-4 + 2) * A00_; // Side
		BC2_ = (-pow(mu_,2)*18 + (2*sigma1_*k_/pow(h_,2))*-4 + 2) * A00_; // Corner
	}
	
	//// Previous time step (C) coeffients
	C00_ = (-(2*sigma1_*k_/pow(h_,2))*-4 - (1-sigma0_*k_))  * A00_;
	C01_ = -(2*sigma1_*k_/pow(h_,2))  * A00_;
	
	//input force coefficient
	d0_ = pow(k_,2)/(rho_*H_*pow(h_,2))*(1/(1+k_*sigma0_))*A00_;
	
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\


void FD_Plate::setup(double samprate, bool bctype)
{
	
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Motion Coefficients
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	D_ = (E_*(pow(H_, 3)))/(12*(1-pow(nu_,2)));
	kappa_ = sqrt(D_ / (rho_*  H_) );
	
	SR_ = samprate;				// internal class sampling rate
	k_ = 1/SR_;					// time step
	
	setLoss(8,.75);
	setGrid();
	setCoefs(bctype);
	
	//Set Input and Output Indeces
	li_ = (Ny_*(ctr_[1]*Nx_)) + (ctr_[0]*Ny_);
	lo_ = (Ny_*(rp_[1]*Nx_)) +  (rp_[0]*Ny_);
	
	//	Update flags
	setupFlag_ = true;
	bcFlag_  = bctype;
	
	
	//	Allocate and Clear State Memory
	//	if (u_ || u1_ || u2_)
	//	{
	//		delete[] u_;
	//		delete[] u1_;
	//		delete[] u2_;
	//	}
	
	u_ = new double[ss_];
	u1_ = new double[ss_];
	u2_ = new double[ss_];
	std::fill(u_,  u_+ss_,  0);
	std::fill(u1_, u1_+ss_, 0);
	std::fill(u2_, u2_+ss_, 0);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Print Scheme Info
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void FD_Plate::printCoefs()
{
	printf("--- Coefficient Info --- \n\n");
	printf("Loss A		: %.4fm \n", A00_);
	printf("Centre B    : %.4fm \n", B00_);
	printf("1-Grid B    : %.4fm \n", B01_);
	printf("2-Grid B	: %.4fm \n", B02_);
	printf("Diagonal B  : %.4fm \n", B11_);
	printf("Centre C	: %.4fm \n", C00_);
	printf("1-Grid C    : %.4fm \n", C01_);
	printf("Side Bound	: %.4fm \n", BC1_);
	printf("Cornr Bound : %.4fm \n\n", BC2_);
}

void FD_Plate::printInfo()
{
	printf("--- Scheme Info --- \n\n");
	printf("Size			: %.1f m2 \n", Nx_*h_*Ny_*h_);
	printf("Thickness(mm)	: %.0f mm \n", H_*1e3);
	printf("Grid X-Ax		: %d \n", Nx_);
	printf("Grid Y-Ax		: %d \n", Ny_);
	printf("Total Ps		: %d \n", ss_);
	printf("In_cell			: %d\n", li_);
	printf("Out_cell		: %d\n", lo_);
	printf("TimeStep		: %.2e\n", k_);
	printf("SampRate		: %.2e\n", SR_);
	printf("Youngs			: %.2e\n", E_);
	printf("Sigma 0			: %f\n", sigma0_);
	printf("Sigma 1			: %f\n", sigma1_);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Update Plate State
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// TODO: update should not be run before Setup()
// Do a check to ensure

void FD_Plate::updateScheme()
{
	int cp;
	// Internal Gride Points
	for(int xi=Nx_-4; xi--; )
	{
		for(int yi=Ny_-4;yi--; )
		{
			cp = (yi+2)+((xi+2) * Ny_); // current point
			
			u_[cp] = B00_*u1_[cp] +
			B01_*( u1_[cp-1] + u1_[cp+1] + u1_[cp-Ny_] + u1_[cp+Ny_] ) +
			B02_*( u1_[cp-2] + u1_[cp+2] +u1_[cp-(2*Ny_)] + u1_[cp+(2*Ny_)] ) +
			B11_*( u1_[cp-1-Ny_] + u1_[cp+1-Ny_] +u1_[cp+1+Ny_] + u1_[cp-1+Ny_] ) +
			C00_*u2_[cp] +
			C01_*( u2_[cp-1] + u2_[cp+1] + u2_[cp-Ny_] + u2_[cp+Ny_] );
		}
	}
	
	// Update Side Boundaries
	//X-Axis
	
	for(int xi=Nx_-4; xi--; )
	{
		//North
		cp = 1+((xi+2) * Ny_); // current point
		u_[cp]  = BC1_*u1_[cp] +
		B01_*( u1_[cp+1] + u1_[cp-Ny_] + u1_[cp+Ny_] ) +
		B02_*( u1_[cp-2] + u1_[cp+2] + u1_[cp-(2*Ny_)] + u1_[cp+(2*Ny_)] ) +
		B11_*( u1_[cp+1-Ny_] + u1_[cp+1+Ny_] ) +
		C00_*u2_[cp] +
		C01_*( u2_[cp+1] + u2_[cp-Ny_] + u2_[cp+Ny_] );
		
		//South
		cp = Ny_-2 +((xi+2) * Ny_); // current point
		u_[cp]  = BC1_*u1_[cp] +
		B01_*( u1_[cp-1] + u1_[cp-Ny_] + u1_[cp+Ny_] ) +
		B02_*( u1_[cp-2] + u1_[cp-(2*Ny_)] + u1_[cp+(2*Ny_)] ) +
		B11_*( u1_[cp-1-Ny_] + u1_[cp-1+Ny_] ) +
		C00_*u2_[cp] +
		C01_*( u2_[cp-1] + u2_[cp-Ny_] + u2_[cp+Ny_] );
	}
	
	// Y-Axis
	
	for(int yi=Ny_-4;yi--; )
	{
		//West
		cp = yi+Ny_+2; // current point
		u_[cp]  = BC1_*u1_[cp] +
		B01_*( u1_[cp-1] + u1_[cp+1] + u1_[cp+Ny_] ) +
		B02_*( u1_[cp-2] + u1_[cp+2] + u1_[cp+(2*Ny_)] ) +
		B11_*( u1_[cp+1+Ny_] + u1_[cp-1+Ny_] ) +
		C00_*u2_[cp] +
		C01_*( u2_[cp-1] + u2_[cp+1] + u2_[cp+Ny_] );
		
		//East
		cp = (yi+2) + Ny_*(Nx_-2); // current point
		u_[cp]  = BC1_*u1_[cp] +
		B01_*( u1_[cp-1] + u1_[cp+1] + u1_[cp-Ny_] ) +
		B02_*( u1_[cp-2] + u1_[cp+2] +u1_[cp-(2*Ny_)] ) +
		B11_*( u1_[cp-1-Ny_] + u1_[cp+1-Ny_] ) +
		C00_*u2_[cp] +
		C01_*( u2_[cp-1] + u2_[cp+1] + u2_[cp-Ny_] );
	}
	
	// Corner Boundaries
	
	cp = Ny_+1;
	u_[cp] = BC2_*u1_[cp] +
	B01_*( u1_[cp-1] + u1_[cp+1] + u1_[cp-Ny_] + u1_[cp+Ny_] ) +
	B02_*( u1_[cp+2] + u1_[cp+(2*Ny_)] ) +
	B11_*( u1_[cp-1-Ny_] + u1_[cp+1-Ny_] +u1_[cp+1+Ny_] + u1_[cp-1+Ny_] ) +
	C00_*u2_[cp] +
	C01_*( u2_[cp-1] + u2_[cp+1] + u2_[cp-Ny_] + u2_[cp+Ny_] );
	
	cp = 2*(Ny_-1);
	u_[cp] = BC2_*u1_[cp] +
	B01_*( u1_[cp-1] + u1_[cp+1] + u1_[cp-Ny_] + u1_[cp+Ny_] ) +
	B02_*( u1_[cp-2] + u1_[cp+(2*Ny_)] ) +
	B11_*( u1_[cp-1-Ny_] + u1_[cp+1-Ny_] +u1_[cp+1+Ny_] + u1_[cp-1+Ny_] ) +
	C00_*u2_[cp] +
	C01_*( u2_[cp-1] + u2_[cp+1] + u2_[cp-Ny_] + u2_[cp+Ny_] );
	
	cp = Ny_*(Nx_-2)+1;
	u_[cp] = BC2_*u1_[cp] +
	B01_*( u1_[cp-1] + u1_[cp+1] + u1_[cp-Ny_] + u1_[cp+Ny_] ) +
	B02_*( u1_[cp+2] + u1_[cp-(2*Ny_)] ) +
	B11_*( u1_[cp-1-Ny_] + u1_[cp+1-Ny_] +u1_[cp+1+Ny_] + u1_[cp-1+Ny_] ) +
	C00_*u2_[cp] +
	C01_*( u2_[cp-1] + u2_[cp+1] + u2_[cp-Ny_] + u2_[cp+Ny_] );
	
	cp = Ny_*(Nx_-1) - 2;
	u_[cp] = BC2_*u1_[cp] +
	B01_*( u1_[cp-1] + u1_[cp+1] + u1_[cp-Ny_] + u1_[cp+Ny_] ) +
	B02_*( u1_[cp-2] + u1_[cp-(2*Ny_)] ) +
	B11_*( u1_[cp-1-Ny_] + u1_[cp+1-Ny_] +u1_[cp+1+Ny_] + u1_[cp-1+Ny_] ) +
	C00_*u2_[cp] +
	C01_*( u2_[cp-1] + u2_[cp+1] + u2_[cp-Ny_] + u2_[cp+Ny_] );
	
	// swap pointers
	dummy_ptr_ = u2_; u2_ = u1_; u1_ = u_; u_ = dummy_ptr_;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

// get the output from the plate based on the readout positino and type.
// Will need to implement interpolation, especially f this is meant to be a
// free moving read-out.

double FD_Plate::getOutput(bool outType)
{
	if(outType){
		return (u1_[lo_]- u2_[lo_])*SR_; // Velocity out
	}
	else{
		return u1_[lo_]; // Amplitude out
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void FD_Plate::getStereoOutput(bool outType, double &Left, double &Right)
{
	if(outType){// Velocity out
		Left = (u1_[lo_l_]- u2_[lo_l_])*SR_;
		Right = (u1_[lo_r_]- u2_[lo_r_])*SR_;
	}
	else{
		Left = u1_[lo_l_]; // Amplitude out
		Right = u1_[lo_r_];
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
// Set the readout position on the plate

void FD_Plate::setOutput(double xcoord, double ycoord)
{
	int readoutpos = floor((Ny_*(xcoord*Nx_)) +  (ycoord*Ny_));
	
	// Ensure that read out position is avalid grid point
	if(readoutpos % Ny_ == 0){readoutpos++;}
	else if(readoutpos % (Ny_) == (Ny_-1)){readoutpos--;}
	if(readoutpos - Ny_ < 0){readoutpos+=Ny_;}
	else if(readoutpos - (ss_-Ny_) > 0){readoutpos-=Ny_;}
	
	lo_ = readoutpos;
}

void FD_Plate::setStereoOutput(double lxcoord, double lycoord)
{
	// For simplicity, this mirrors left and right output
	// down the middle of the y-axis
	double rxcoord = lxcoord; double rycoord = 1-lycoord;
	
	//TODO: CHECK IF ON A VALID GRID POINT
	if((lxcoord*Nx_)-1 < 1 || (lxcoord*Nx_)-1 > Nx_-2){
		//Do something to ensure it is not on a zero point
	}
	if((rxcoord*Nx_)-1 < 1 || (rxcoord*Nx_)-1 > Nx_-2){
		//Do something to ensure it is not on a zero point
	}
	if((lycoord*Ny_)-1 < 1 || (lycoord*Ny_)-1 > Nx_-2){
		//Do something to ensure it is not on a zero point
	}
	
	if((rycoord*Ny_)-1 < 1 || (rycoord*Ny_)-1 > Nx_-2){
		//Do something to ensure it is not on a zero point
	}
	
	lo_l_ = (Ny_*(lxcoord*Nx_)) +  (lycoord*Ny_);
	lo_r_ = (Ny_*(rxcoord*Nx_)) +  (rycoord*Ny_);
	
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
// Method sets the output type, either velocity or, amplitude. Can probably be
// intergrated into the get output method.

void FD_Plate::setOutType(bool outtype)
{
	outFlag_ = outtype;  // set output to velocity amplitude
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
// Method will set the profile of the input force, for use in JUCE.
// Currently not interpolated, will need to look into that.
// Array values will them selves be multiplied by a rasied cosine/half-cosine

void FD_Plate::setInitialCondition()
{
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Excitation Force
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// raised cosine in 2D
	for(int xi=1;xi<Nx_-1;xi++)
	{
		const double X = xi*h_;
		
		for(int yi=1;yi<Ny_-1;yi++)
		{
			const int cp = yi+(xi * Ny_);
			const double Y = yi*h_;
			const double dist = sqrt(pow(X-(ctr_[0]*Lx_),2) + pow(Y-(ctr_[1]*Ly_),2));
			const double ind = sgn((wid_*0.5)-dist);			// displacement (logical)
			const double rc = .5*ind*(1+cos(2*pi*dist/wid_)); // displacement
			u2_[cp] = u0_*rc;
			u1_[cp] = v0_*k_*rc;
		}
	}
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void FD_Plate::addForce(double force)
{
	u1_[li_] += d0_* force;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void FD_Plate::addStrike()
{
	//TODO: Function to add a 2D rasied cosine multiplied by a gain rc in
	// displacement, dependent on strike velocity
	// know that strike has begun or ended (insert flag)
	// calculate gain rc based on velocity.
	// prepare for multiple strikes, potentially overlapping.
	
	// Need to work out in advance which indeces will actually be affected by a strike.
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

double FD_Plate::reverb(double force)
{
	updateScheme();
	addForce(force);
	//	return getOutput(true);
	return getInterpOut();
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

double FD_Plate::getInterpOut()
{
	const int order = interpOrder_;
	double interpOut = 0;

	for(int xi = 0; xi<order;xi++)
	{
		for(int yi = 0; yi<order;yi++)
		{
			// out of bound test to mark point as zero
			if(yInterpIndeces_[yi] < 0 || yInterpIndeces_[yi] > Ny_ || xInterpIndeces_[xi] < 0 || xInterpIndeces_[xi] > Nx_)
			{
//				hiResValue += 0;
			}
			else
			{
				int cp = yInterpIndeces_[yi] + ((xInterpIndeces_[xi]) * (Ny_-1));
				interpOut += (interpLookTable_[xi][xAlphaIndex_] * interpLookTable_[yi][yAlphaIndex_]) *
							 (u1_[cp]-u2_[cp]) * SR_;
			}
		}
	}
	return interpOut;
}


//======================================================================
//	LINEAR INTERP (COMMENT OUT)
//======================================================================

//double FD_Plate::getInterpOut()
//{
//	
//		const int cp = interpZeroIndex_;
//	
//		const double interpOut = (((1-interpAlphaX_)* (1-interpAlphaY_) * u1_[cp])	    +
//						   ((1-interpAlphaX_)* (interpAlphaY_)   * u1_[cp+1])	+
//						   ((interpAlphaX_)  * (1-interpAlphaY_) * u1_[cp+Ny_]) +
//						   ((interpAlphaX_)  * (interpAlphaY_)   * u1_[cp+1+Ny_])
//	
//						   -(((1-interpAlphaX_)* (1-interpAlphaY_) * u2_[cp])	+
//						   ((1-interpAlphaX_)* (interpAlphaY_)   * u2_[cp+1])	+
//						   ((interpAlphaX_)  * (1-interpAlphaY_) * u2_[cp+Ny_]) +
//						   ((interpAlphaX_)  * (interpAlphaY_)   * u2_[cp+1+Ny_])))*SR_;
//	
//	return interpOut;
//}
//======================================================================
//	LINEAR INTERP (COMMENT OUT)
//======================================================================


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//Split this so that it is arbitrary of coordinate, so input is xCoord*Nx_

void FD_Plate::setInterpOut(const double xCoord, const double yCoord)
{
//	lo_ = (Ny_*(xCoord*Nx_)) +  (yCoord*Ny_);

	const int order = interpOrder_;
	for(int i = 0; i<order;i++)
	{
		xInterpIndeces_[i] = (i+1 - order/2.) + floor(xCoord*(Nx_));
		yInterpIndeces_[i] = (i+1 - order/2.) + floor(yCoord*(Ny_-1));
//		printf("Y: %d \t X: %d \t I: %d\n",yInterpIndeces_[i],xInterpIndeces_[i],int(i+1 - order/2.));
	}
//	printf("\n");
//	printf("Y: %d \t X: %d\n",yInterpIndeces_[i],xInterpIndeces_[i]);
	
	
	interpPointY_ = (yCoord*(Ny_-1));
	interpPointX_ = (xCoord*Nx_);

//	printf("Y: %.2f \t X: %.2f\n",interpPointY_,interpPointX_);
//	printf("lo: %d \t Y: %f \t Alpha: %.2f \n",yInterpIndeces_[1] +  xInterpIndeces_[1]*Ny_,yCoord,interpPointY_-floor(interpPointY_));
	
	
//		const int order = interpOrder_;
	const int res = interpRes_;
	xAlphaIndex_ = floor((interpPointX_-floor(interpPointX_))*res);
	yAlphaIndex_ = floor((interpPointY_-floor(interpPointY_))*res);
}


//======================================================================
//	LINEAR INTERP (COMMENT OUT)
//======================================================================
//void FD_Plate::setInterpOut(const double xCoord, const double yCoord)
//{
//	interpAlphaY_ = (yCoord*(Ny_)) - floor(yCoord*(Ny_));
//	interpAlphaX_ = (xCoord*(Nx_)) - floor(xCoord*(Nx_));
//	interpZeroIndex_ = int( floor(yCoord*(Ny_-1)) + floor(xCoord*(Nx_))* (Ny_) );
//	
//	if (interpZeroIndex_ > ((Nx_*Ny_)-1-Ny_))
//	{
//		interpZeroIndex_ = int( floor(yCoord*(Ny_-2)) + (Nx_-2)* (Ny_) );
//		interpAlphaY_ = 0;
//		interpAlphaX_ = 0;
//	}
////	printf("Zero Index: %d\n",interpZeroIndex_);
//}
//======================================================================
//	LINEAR INTERP (COMMENT OUT)
//======================================================================




//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

double **FD_Plate::getInterpLookTable()
{
	//	const int order = 4;
	//	const int res = 100;
	const int order = interpOrder_;
	const int res = interpRes_;
	double **alphaTable = new double*[order];
	
	for(int i=0;i<order;i++)
	{
		alphaTable[i] = new double[res];
		std::fill(alphaTable[i], alphaTable[i]+res-1, 1);
	}
	
	double *polynomial_normaliser = new double [order];
	std::fill(polynomial_normaliser, polynomial_normaliser+order, 1);
	double *alphas = new double [res];
	
	for(int i = 0; i<res;i++)
	{
		alphas[i] = (i/float(res)) - 0.5;
	}
	
	double *anchors = new double [order];
	
	if ((order % 2)==0)
	{
		for(int i = 0; i<order;i++)
		{
			anchors[i] = -(order - 1)*0.5 + i;
			
		}
	}
	else
	{
		for(int i = 0; i<order;i++)
		{
			anchors[i] = (-(order)*0.5) + i;
		}
	}
	
	for (int q = 0; q<res;q++) // loop for every value of alpha
	{
		for (int j = 0; j<order;j++) // loop for sub polynomial
		{
			for (int m = 0; m < order; m++) //loop for each point in subpoly
			{
				if (m != j)
				{
					if (q == 0)
					{
						polynomial_normaliser[j] = polynomial_normaliser[j]*(anchors[j]-anchors[m]);
					}
					alphaTable[j][q] *= (alphas[q]-anchors[m]);
					
				}
			}
			alphaTable[j][q] /= polynomial_normaliser[j];
		}
	}
	delete[] polynomial_normaliser;
	delete[] alphas;
	delete[] anchors;
	return alphaTable;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
