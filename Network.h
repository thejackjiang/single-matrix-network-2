// Network.h: interface for the Network class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_NETWORK_H__8E7C932B_D833_4E1F_9EDC_ED09AFCF876A__INCLUDED_)
#define AFX_NETWORK_H__8E7C932B_D833_4E1F_9EDC_ED09AFCF876A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#define MAX_NET_DIMENSION 100
#define MAX_NET_INPUTS 10		// Cannot exceed the size of the net
#define MAX_NET_OUTPUTS 10		// Cannot exceed the size of the net
#define MAX_DUMMY_STRING_LENGTH 30

class Network  
{
public:
	Network();								// default network constructor.
	Network::Network( char * file_name );	// Construct the network from a stored file.
	Network::Network(int inputs, int interneurons, int outputs, char * file_name ); // construct a blank network to user specified size and write it to a file for later editing
	virtual ~Network();

	// Members---------------------
	// Additions to the members should be added to the read and write netork functions

private:
	int numberOfInputs;
	int numberOfOutputs;
	int numberOfInterNeurons;
	int networkDimension;

	double neuronActivation[MAX_NET_DIMENSION];		// Individual neurons have individual activation levels that are tracked through timesteps an are not visibile as output ( transformed to output by a function)
	double neuronOutput[MAX_NET_DIMENSION];			// The output of individual neurons, used as inputs to the other neurons in the network throught the connection weights matrix
	double neuronThresholds[MAX_NET_DIMENSION];		// Individual Neurons each have a speficifed activation threshold
	double neuronLearningRate[MAX_NET_DIMENSION];	// Individual Neurons each have a speficifed learning rate -- rate of change of connection strength per time step
	short int neuronRefractoryState[MAX_NET_DIMENSION];	// Individual Neurons each have a speficifed period during which output is blocked -- should be 0 or greater.
	double neuronWeightTotal[MAX_NET_DIMENSION];	// Individual Neurons each have a speficifed total weight strength in their input connections.
	double networkWeights[MAX_NET_DIMENSION*MAX_NET_DIMENSION];
	double networkInputs[MAX_NET_INPUTS];
	double networkOutputs[MAX_NET_OUTPUTS];
	short int  plasticWeightsMask[MAX_NET_DIMENSION*MAX_NET_DIMENSION]; // a filter. Plastic weights are = 1, fixed = 0. THis allows for the specification of some fixed and some plastic weights in the same neuron. This could be a binary array ( type bool) to save space. 

	// Functions -------------------------

	void Network::instantiateDefaultNetwork( );
	void Network::setNetworkOuput( );
	void Network::copyNeuronActivationsToNeuronOutputs( );
	void Network::copyNetworkInputsToInputNeuronOutputs( );
	void Network::thresholdNeuronOutputs( );
	void Network::setNeuronOutput( double value );
	void Network::setNeuronThresholds( double value );
	void Network::setNeuronLearningRate( double value );
	int Network::setNeuronRefractoryState( int value );
	void Network::setPlasticWeightsMask( short int value ); // in general it is good to set this to 1 and let the learning rate determine plasticity.  This is to be used for special cases
	void Network::setNeuronActivation( double value );
	void Network::setNetworkOutputs( double value );
	void Network::networkActivation( );
	void Network::hebbianWeightUpdate( );
	void Network::hebbianExcitatoryWeightUpdate( );
	void Network::hebbianInhibitoryWeightUpdate( );
	void Network::normalizeNeuronWeights( );			// Update weight totals to neuron-specific values
	void Network::normalizeNeuronWeights( double value );	// Uptdate weight totals to specificed values
	void Network::normalizeNonDiagonalNeuronWeights( );
	void Network::normalizeNonDiagonalInhibitoryNeuronWeights( );
	void Network::normalizeNonDiagonalExcitatoryNeuronWeights( );
	void Network::setNeuronWeightTotal( double value);
	int Network::computeWeightIndex( int source_neuron_number, int target_neuron_number );

public:
	void Network::cycleNetwork( );
	void Network::cycleNetworkNormalizeHebbianLearning( );
	void Network::printNetworkOuput( );
	void Network::printNetworkOutputState( );
	void Network::setNetworkWeightsDiagonalRange( double value, int start_row_col, int end_row_col );
	void Network::setNetworkWeightsUpperLowerTriangleAndDiagonal( double diagonal_value, double upper_triangle_value, double lower_triangle_value);
	void Network::setNetworkWeightsRectangle( double value, int start_row, int end_row, int start_column, int end_column );
	void Network::setNetworkWeightsUpperTriangle( double value, int start_row, int end_row, int start_column, int end_column );
	void Network::setNetworkWeightsLowerTriangle( double value, int start_row, int end_row, int start_column, int end_column );
	void Network::writeNetworkOutputStateToFile( char * file_name );
	void Network::writeNetworkActivationStateToFile( char * file_name );
	void Network::writeNetworkWeightsToFile( char * file_name );
	void Network::setNetworkInput( double *vector);
	void Network::getNetworkOuput( double * vector );
	int Network::readNetworkFromFile( char * file_name );
	int Network::writeNetworkToFile( char * file_name );
	void Network::setNetworkWeights( double value );
	void Network::PrintNetworkState( );

};

#endif // !defined(AFX_NETWORK_H__8E7C932B_D833_4E1F_9EDC_ED09AFCF876A__INCLUDED_)
