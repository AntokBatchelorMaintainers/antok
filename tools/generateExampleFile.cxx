#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"

#include <iostream>

const unsigned int numberOfEvents = 10000;

std::vector<double> getGaussianData( double mean, double sigma )
{
	std::vector<double> data;
	data.resize(numberOfEvents);
	for( unsigned int i = 0; i < numberOfEvents; i++ )
	{
		data[i] = gRandom->Gaus( mean, sigma );
	}
	return data;
}

std::vector<double> getBoxData( double mean, double width )
{
	std::vector<double> data;
	data.resize(numberOfEvents);
	for( unsigned int i = 0; i < numberOfEvents; i++ )
	{
		data[i]= gRandom->Uniform( mean - width, mean +width );
	}
	return data;
}

void generateExampleFile()
{
	TFile* outFile = TFile::Open("example.root", "RECREATE");
	outFile->cd();

	std::vector<double> vertexX;
	std::vector<double> vertexY;
	std::vector<double> vertexZ;

	std::vector<double> beamMomentumX;
	std::vector<double> beamMomentumY;
	std::vector<double> beamMomentumZ;

	std::vector<double> scatteredMomentumX;
	std::vector<double> scatteredMomentumY;
	std::vector<double> scatteredMomentumZ;

	TTree* outTree = new TTree("exampleEvents", "exampleEvents");
	outTree->Branch("VertexX", &vertexX );
	outTree->Branch("VertexY", &vertexY );
	outTree->Branch("VertexZ", &vertexZ );

	outTree->Branch("BeamMomentumX", &beamMomentumX,"beamMomentumX/D");
	outTree->Branch("BeamMomentumY", &beamMomentumY,"beamMomentumY/D");
	outTree->Branch("BeamMomentumZ", &beamMomentumZ,"beamMomentumZ/D");

	outTree->Branch("ScatteredMomentumX", &scatteredMomentumX,"scatteredMomentumX/D");
	outTree->Branch("ScatteredMomentumY", &scatteredMomentumY,"scatteredMomentumY/D");
	outTree->Branch("ScatteredMomentumZ", &scatteredMomentumZ,"scatteredMomentumZ/D");

	vertexX = getGaussianData( 0, 3.2    );
	vertexY = getGaussianData( 0, 2.8    );
	vertexZ = getBoxData     ( 0, 22.42  );

	beamMomentumX = getGaussianData( 0.1, 0.5 );
	beamMomentumY = getGaussianData( 0.2, 0.7 );
	beamMomentumZ = getGaussianData( 200, 2.4 );

	scatteredMomentumX = getGaussianData( 0.2, 0.4 );
	scatteredMomentumY = getGaussianData( 0.4, 0.7 );
	scatteredMomentumZ = getGaussianData( 198, 3.1 );

	outTree->Fill();
	outTree->Write();
	outFile->Close();
}

int main()
{
	generateExampleFile();
}
