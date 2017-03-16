#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"

#include <iostream>

const unsigned int numberOfEvents = 10000;

double getGaussianData( double mean, double sigma )
{
	return gRandom->Gaus( mean, sigma );
}

double getBoxData( double mean, double width )
{
	double z;
	do
	{
		z = gRandom->Exp(600);
	}
	while( z > mean + 2 * width );
	return (mean - width * 0.5) + z;
}

unsigned int getIntegerData( unsigned int from, unsigned int count )
{
	return from + gRandom->Integer( from + count );
}

void generateExampleFile()
{
	TFile* outFile = TFile::Open("example.root", "RECREATE");
	outFile->cd();

	unsigned int runNumber;

	double vertexX;
	double vertexY;
	double vertexZ;

	double beamGradientdXdZ;
	double beamGradientdYdZ;

	double recoilMomentumX;
	double recoilMomentumY;
	double recoilMomentumZ;
	double recoilEnergy;

	double scatteredMomentumX1;
	double scatteredMomentumY1;
	double scatteredMomentumZ1;

	double scatteredMomentumX2;
	double scatteredMomentumY2;
	double scatteredMomentumZ2;

	double scatteredMomentumX3;
	double scatteredMomentumY3;
	double scatteredMomentumZ3;

	TTree* outTree = new TTree("exampleEvents", "exampleEvents");
	outTree->Branch("RunNumber", &runNumber , "RunNumber/I" );

	outTree->Branch("VertexX", &vertexX, "VertexX/D" );
	outTree->Branch("VertexY", &vertexY, "VertexY/D" );
	outTree->Branch("VertexZ", &vertexZ, "VertexZ/D" );

	outTree->Branch("BeamGradientdXdZ", &beamGradientdXdZ, "BeamGradientdXdZ/D" );
	outTree->Branch("BeamGradientdYdZ", &beamGradientdYdZ, "BeamGradientdYdZ/D" );

	outTree->Branch("RecoilMomentumX", &recoilMomentumX, "RecoilMomentumX/D" );
	outTree->Branch("RecoilMomentumY", &recoilMomentumY, "RecoilMomentumY/D" );
	outTree->Branch("RecoilMomentumZ", &recoilMomentumZ, "RecoilMomentumZ/D" );
	outTree->Branch("RecoilEnergy"   , &recoilEnergy   , "RecoilEnergy/D"    );

	outTree->Branch("ScatteredMomentumX_1", &scatteredMomentumX1, "ScatteredMomentumX_1/D" );
	outTree->Branch("ScatteredMomentumY_1", &scatteredMomentumY1, "ScatteredMomentumY_1/D" );
	outTree->Branch("ScatteredMomentumZ_1", &scatteredMomentumZ1, "ScatteredMomentumZ_1/D" );

	outTree->Branch("ScatteredMomentumX_2", &scatteredMomentumX2, "ScatteredMomentumX_2/D" );
	outTree->Branch("ScatteredMomentumY_2", &scatteredMomentumY2, "ScatteredMomentumY_2/D" );
	outTree->Branch("ScatteredMomentumZ_2", &scatteredMomentumZ2, "ScatteredMomentumZ_2/D" );

	outTree->Branch("ScatteredMomentumX_3", &scatteredMomentumX3, "ScatteredMomentumX_3/D" );
	outTree->Branch("ScatteredMomentumY_3", &scatteredMomentumY3, "ScatteredMomentumY_3/D" );
	outTree->Branch("ScatteredMomentumZ_3", &scatteredMomentumZ3, "ScatteredMomentumZ_3/D" );

	for( unsigned int i = 0; i < numberOfEvents; ++i )
	{
		runNumber = getIntegerData( 10000, numberOfEvents );

		vertexX = getGaussianData(  0.1,  3.32 );
		vertexY = getGaussianData( -0.2,  2.82 );
		vertexZ = getBoxData     (    0, 22.42 );

		beamGradientdXdZ = getGaussianData( 0.001, 0.005 );
		beamGradientdYdZ = getGaussianData( 0.002, 0.007 );

		recoilMomentumX = getGaussianData( 0.1, 0.1 );
		recoilMomentumY = getGaussianData( 0.2, 0.3 );
		recoilMomentumZ = getGaussianData( 0.3, 0.2 );
		recoilEnergy    = TMath::Abs(getGaussianData(   5, 0.1 ));

		scatteredMomentumX1 = getGaussianData( 0.2, 0.4 );
		scatteredMomentumY1 = getGaussianData( 0.4, 0.7 );
		scatteredMomentumZ1 = getGaussianData( 75, 3.1  );

		scatteredMomentumX2 = getGaussianData( 0.2, 0.4 );
		scatteredMomentumY2 = getGaussianData( 0.4, 0.7 );
		scatteredMomentumZ2 = getGaussianData( 80, 3.1  );

		scatteredMomentumX3 = getGaussianData( 0.2, 0.4 );
		scatteredMomentumY3 = getGaussianData( 0.6, 0.7 );
		scatteredMomentumZ3 = getGaussianData(  55, 3.1 );

		outTree->Fill();
	}
	outTree->Write();

	outFile->Close();
}

int main()
{
	generateExampleFile();
}
