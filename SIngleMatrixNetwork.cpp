// SIngleMatrixNetwork.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Network.h"
#include <math.h>
#include <stdio.h>

int main(int argc, char* argv[])
{
	Network jack(1, 2, 2, "jacktest.txt");

	int i = 0;
	double input[6];
	while (i < 1000) {


		input[0] = 0.0;			// Set all inputs off
		if (i > 300 && i < 700)	input[0] = 1.0;		// Command jaw close 
		if (i > 600 && i < 800) input[0] = 1.0;		// command LEFT jaw open
		if (i > 650 && i < 850) input[0] = 1.0;		// command RIGHT jaw open
		jack.setNetworkInput(input);

		jack.cycleNetwork();
		jack.printNetworkOutputState();

		jack.writeNetworkOutputStateToFile("jacktest.txt");
	}
		
	//double input[6];
	//int i;
	//int j;

//	Network ted(3,10,2,"blank_network.txt");

//	Network fred;
	//Network fred("antJaw_05.txt");

	//fred.PrintNetworkState();


	// Initialize the weights to produce a test recurrent network
/*
	fred.setNetworkWeights(0.0);					// initialize all connections to zero
	fred.setNetworkWeightsRectangle(0.3,3,7,0,3);	// set feedforward connections from inputs to interneurons 
	fred.setNetworkWeightsRectangle(0.5,7,9,3,7);	// set feedforward connections from interneurons to output neurons 
	fred.setNetworkWeightsRectangle(-0.25,3,7,3,7);	// set recurrent connections from interneurons to interneurons
	fred.setNetworkWeightsDiagonalRange(0.1,0,9);	// set network autapses
  */

/*	
	//--- Identity Matrix ----
	fred.setNetworkWeights(0.0);
	fred.setNetworkWeightsDiagonalRange(1,0,8);
*/

	/*fred.writeNetworkOutputStateToFile( "testOutput10.txt" );

	i = 0;
	while( i < 1000 ){

		
		input[0] = input[1] = input[2]	= 0.0;			// Set all inputs off
		if( i> 300 && i < 700)	input[0] = 1.0;		// Command jaw close 
		if( i> 600 && i < 800) input[1] = 1.0;		// command LEFT jaw open
		if( i> 650 && i < 850) input[2] = 1.0;		// command RIGHT jaw open
		// Vary duration of command 
		// Get offset to work - correct resting state
		// close and then sensory input
		// approach close, then before it's closed get sensory input to see if it stops



/*
		//-----------------------------------------------------------
		for(j = 0; j < 6; ++j ){	// temporary input function (STEP RESPONSE)
			input[j]  =  sin((i + j*3)*0.01*i);
//			input[j]  = 1.0;
		}
		//---------------------------------------------------------------------
*/
/*		jack.setNetworkInput( input );

		jack.cycleNetwork();
		jack.printNetworkOutputState( );

		jack.writeNetworkOutputStateToFile( "jacktest.txt" );

//		fred.writeNetworkToFile("test.txt");
	//	fred.writeNetworkWeightsToFile("weights.txt");
	//	++i;
	//}

	//Test the read-write fucntions -------------
//	fred.writeNetworkToFile("temp.txt");
//	fred.readNetworkFromFile("temp.txt");
*/
	return 0;

}
