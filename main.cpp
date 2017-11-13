//
//  PLATE REVERB
//
//
//  Created by mhamilt on 20/10/2017.
//  Copyright © 2017 mhamilt. All rights reserved.
//
//	Use Plate Class as a plate reverb unit

#include <iostream>
#include <cstring>
#include "PlateClass.cpp"
#include "AudioOut.h"	 // For .wav header, write and playback

int main (int argc, const char *argv[]) {

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Parameters
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// START EDIT HERE

	// High and low Frequency Loss
	// High freqeuncy in range  0 < x < 1
	const double lt60 = 8; const double ht60 = .1;

	// Set to 0 for Command Line Build
	// Or set 1 and change default file names below.
	const int fndebug = 1;

	char definfile[] = "/Users/admin/Documents/Masters/PBMMI/Audio_Examples/AfricanDrumming.wav";
	char defoutfile[] = "/Users/admin/Downloads/PlateReverb.wav";

	// STOP EDIT HERE

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Input
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	int inNf; int SR; int filenamelength;
	char *inputfname;

	switch (fndebug) { //This is just to allow debuggging of for arguements via command line

  case 0:	//CommandLine

			if(!argv[1]){
				printf("NO FILE DETECTED\n");
				return -1;
			}

			filenamelength = int(strlen(argv[1]));
			inputfname = new char[filenamelength+1]();
			strncpy(inputfname, argv[1], filenamelength);
			break;

  case 1:	// debug
			inputfname = definfile;
			break;

  default:  inputfname = definfile;
	}

	double *inwav = readWav(inputfname, &inNf, &SR);
	printWavHeader(inputfname);

	if(!inwav){
		printf("NO FILE DETECTED\n");
		return -1;
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Set Output File Name
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	int Nf = int(lt60*ht60) + inNf;
	double *out;
	out = new double[Nf];
	char *outputfname;


	switch (fndebug) {

  case 0:	//CommandLine
			if(!argv[2]){
				const char *homedir = getenv("HOME");
				if(!homedir){
					printf("Couldn't find Home Directory. Enter a filename\n");
					return 1;
				}
				printf("NO FILE NAME ENTERED\nUSING DEFAULT FILE DESTINATION\n");
				const char *def_fname = "/Downloads/PlateReverb.wav";
				filenamelength = int(strlen(homedir)) + int(strlen(def_fname));
				outputfname = new char[filenamelength+1]();
				strncpy(outputfname,homedir, int(strlen(homedir)));
				strcat(outputfname, def_fname);
			}
			else{
				filenamelength = int(strlen(argv[2]));
				outputfname = new char[filenamelength+1]();
				strncpy(outputfname, argv[2], filenamelength);
			}
			break;


  case 1:	// debug
			outputfname = defoutfile;
			break;

  default: outputfname = defoutfile;
	}


	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Plate Setup
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	FD_Plate PlateTest;

	//input is SR, and boundary condition type
	PlateTest.Setup(SR,1);

	PlateTest.setLoss(lt60,ht60);

	//	Print info
	PlateTest.printCoefs();
	PlateTest.printInfo();

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Process Loop
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// This is where the magic happens
	//	PlateTest.setInitialCondition();

	// START CLOCK
	std::clock_t start;
	double duration;
	start = std::clock();

	int n;

	for(n=0;n<Nf;n++){
		PlateTest.UpdateScheme();
		if(n<inNf){
			PlateTest.addForce(inwav[n]);
		}
		out[n] = PlateTest.getOutput(1);
	}

	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	std::cout<<"Time Elapsed (seconds): "<< duration <<'\n';

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Output
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	writeWavMS(out, outputfname, Nf, SR);
	delete out;
	playWavMS(outputfname);
	printf("\nComplete...\n");

	std::cout << "SUCCESS" << '\n';
	return 0;

}
