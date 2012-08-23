
#include<assert.h>

#include<plotter.h>

hlib::Plotter::Plotter() {

	std::vector<int> cutmasks;
	std::vector<int> standard_cutmasks;

	standard_cutmasks.push_back(0);
	standard_cutmasks.push_back(256);
	standard_cutmasks.push_back(384);
	standard_cutmasks.push_back(448);
	standard_cutmasks.push_back(480);
	standard_cutmasks.push_back(496);
	standard_cutmasks.push_back(504);
	standard_cutmasks.push_back(508);
	standard_cutmasks.push_back(510);
	standard_cutmasks.push_back(511);

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	_plots.push_back(hlib::Plot(cutmasks, new TH1D("mass_5pi", "mass_5Pi", 500, 0, 7), &XMass));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(128);
	_plots.push_back(hlib::Plot(cutmasks, new TH1D("mom_5pi", "mom_5Pi", 500, 0, 250), &XMom));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(128);
	_plots.push_back(hlib::Plot(cutmasks, new TH1D("calc_beam_E", "calc_beam_E", 500, 0, 250), &CalcBeamE));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(8);
	cutmasks.push_back(264);
	cutmasks.push_back(136);
	cutmasks.push_back(392);
	_plots.push_back(hlib::Plot(cutmasks, new TH1D("rpd_mult", "rpd_mult", 10, 0, 10), &RPDMult));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(16);
	_plots.push_back(hlib::Plot(cutmasks, new TH1D("rpd_P_mass", "rpd_P_mass", 1000, 0, 10), &ProtonMass));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(4);
	_plots.push_back(hlib::Plot(cutmasks, new TH2D("vtx_pos", "vtx_pos", 1000, -5, 5, 1000, -5, 5), &PrimVX, &PrimVY));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(2);
	_plots.push_back(hlib::Plot(cutmasks, new TH1D("vtx_z", "vtx_z", 2000, -200, 200), &PrimVZ));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(64);
	_plots.push_back(hlib::Plot(cutmasks, new TH1D("t_prime", "t_prime", 1000, -5, 5), &TPrim));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(128);
	_plots.push_back(hlib::Plot(cutmasks, new TH1D("delta_phi", "delta_phi", 500, -7, 7), &RPDDeltaPhi));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(1);
	_plots.push_back(hlib::Plot(cutmasks, new TH1D("trigger_mask", "trigger_mask", 15, 0, 15), &TrigMask));
	cutmasks.clear();

};

void hlib::Plotter::fill(const hlib::Event& event, int cutmask) {

	XMass = event.get_pSum().M();
	XMom = event.get_pSum().Energy();
	CalcBeamE = event.get_pBeam().E();
	RPDMult = event.rawData->nbrRPDTracks;
	PrimVX = event.rawData->X_primV;
	PrimVY = event.rawData->Y_primV;
	PrimVZ = event.rawData->Z_primV;
	ProtonMass = event.get_pProton().M();
	TPrim = event.get_tPrime();
	RPDDeltaPhi = event.get_RpdDeltaPhi();
	TrigMask = event.rawData->TrigMask;

	for(unsigned int i = 0; i < _plots.size(); ++i) {
		_plots.at(i).fill(cutmask);
	}

}

void hlib::Plotter::save(TDirectory* dir) {

	dir->cd();
	for(unsigned int i = 0; i < _plots.size(); ++i) {

		hlib::Plot plot = _plots.at(i);
		const char* name = plot.get_template()->GetName();
		TDirectory* new_dir = dir->mkdir(name);
		assert(new_dir != NULL);
		new_dir->cd();
		std::vector<TH1*> histograms = plot.get_histograms();
		for(unsigned int j = 0; j < histograms.size(); ++j) {
			TH1* hist = histograms.at(j);
			hist->Write();
		}

	}

}
