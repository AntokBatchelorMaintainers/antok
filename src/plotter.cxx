
#include<assert.h>

#include<object_manager.h>
#include<plotter.h>

antok::Plotter* antok::Plotter::_plotter = 0;

antok::Plotter* antok::Plotter::instance() {
	if(_plotter == 0) {
		_plotter = new antok::Plotter();
	}
	return _plotter;
}

antok::Plotter::Plotter() {

/*
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

	standard_cutmasks.push_back(48);
	standard_cutmasks.push_back(56);
	standard_cutmasks.push_back(112);
	standard_cutmasks.push_back(120);
	standard_cutmasks.push_back(176);
	standard_cutmasks.push_back(184);
	standard_cutmasks.push_back(240);
	standard_cutmasks.push_back(248);

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	_plots.push_back(antok::Plot(cutmasks, new TH1D("mass_5pi", "mass_5Pi", 500, 0, 7), &XMass));
	antok::Data& data = antok::ObjectManager::instance()->getData();
	_plots.push_back(antok::Plot(cutmasks, new TH1D("mass_5pi_newWay", "mass_5Pi_newWay", 500, 0, 7), data.getDoubleAddr("XMass")));
	cutmasks.clear();



	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(128);
	_plots.push_back(antok::Plot(cutmasks, new TH1D("mom_5pi", "mom_5Pi", 500, 0, 250), &XMom));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(128);
	_plots.push_back(antok::Plot(cutmasks, new TH1D("calc_beam_E", "calc_beam_E", 500, 0, 250), &CalcBeamE));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(8);
	cutmasks.push_back(264);
	cutmasks.push_back(136);
	cutmasks.push_back(392);
	_plots.push_back(antok::Plot(cutmasks, new TH1D("rpd_mult", "rpd_mult", 10, 0, 10), &RPDMult));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(16);
	_plots.push_back(antok::Plot(cutmasks, new TH1D("rpd_P_mass", "rpd_P_mass", 1000, 0, 10), &ProtonMass));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(4);
	_plots.push_back(antok::Plot(cutmasks, new TH2D("vtx_pos", "vtx_pos", 1000, -5, 5, 1000, -5, 5), &PrimVX, &PrimVY));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(2);
	_plots.push_back(antok::Plot(cutmasks, new TH1D("vtx_z", "vtx_z", 2000, -200, 200), &PrimVZ));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(64);
	_plots.push_back(antok::Plot(cutmasks, new TH1D("t_prime", "t_prime", 1000, -5, 5), &TPrim));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(128);
	_plots.push_back(antok::Plot(cutmasks, new TH1D("delta_phi", "delta_phi", 500, -7, 7), &RPDDeltaPhi));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	cutmasks.push_back(1);
	_plots.push_back(antok::Plot(cutmasks, new TH1D("trigger_mask", "trigger_mask", 15, 0, 15), &TrigMask));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	_plots.push_back(antok::Plot(cutmasks, new TH1D("beam_time", "beam_time", 100, -10, 10), &beam_time));
	cutmasks.clear();

	cutmasks.insert(cutmasks.begin(), standard_cutmasks.begin(), standard_cutmasks.end());
	_plots.push_back(antok::Plot(cutmasks, new TH2D("cedar_theta", "cedar_theta", 3000, -300, 300, 3000, -300, 300), &cedarTheta_X, &cedarTheta_Y));
	cutmasks.clear();

	cutmasks.push_back(384);
	_plots.push_back(antok::Plot(cutmasks, new TH1D("delta_phi", "delta_phi", 500, -7, 7), &RPDDeltaPhi));
	_plots.push_back(antok::Plot(cutmasks, new TH1D("delta_phi_fhaas", "delta_phi", 500, -7, 7), &RPDDeltaPhi_fhaas));
	_plots.push_back(antok::Plot(cutmasks, new TH2D("delta_phi_comparison", "delta_phi_comparison", 1000, -7, 7, 1000, -7, 7), &RPDDeltaPhi, &RPDDeltaPhi_fhaas));
	_plots.push_back(antok::Plot(cutmasks, new TH1D("delta_phi_abs", "delta_abs_phi", 500, -7, 7), &RPDDeltaPhiAbs));
	_plots.push_back(antok::Plot(cutmasks, new TH1D("delta_phi_abs_fhaas", "delta_abs_phi", 500, -7, 7), &RPDDeltaPhiAbs_fhaas));
	_plots.push_back(antok::Plot(cutmasks, new TH2D("delta_phi_abs_comparison", "delta_phi_abs_comparison", 1000, -7, 7, 1000, -7, 7), &RPDDeltaPhiAbs, &RPDDeltaPhiAbs_fhaas));
	_plots.push_back(antok::Plot(cutmasks, new TH1D("phi_res", "phi_res", 500, -7, 7), &RPDPhiRes));
	_plots.push_back(antok::Plot(cutmasks, new TH1D("phi_res_fhaas", "phi_res_fhaas", 500, -7, 7), &RPDPhiRes_fhaas));
	_plots.push_back(antok::Plot(cutmasks, new TH2D("phi_res_comparison", "phi_res_comparison", 500, -7, 7, 500, -7, 7), &RPDPhiRes, &RPDPhiRes_fhaas));
*/
};

void antok::Plotter::fill(const antok::Event& event, int cutmask) {

/*
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
	beam_time = event.rawData->beam_time;
	cedarTheta_X = event.rawData->cedarTheta_X;
	cedarTheta_Y = event.rawData->cedarTheta_Y;

	RPDPhiRes = event.get_RpdPhiRes();
	RPDDeltaPhi_fhaas = event.get_RpdDeltaPhi_fhaas();
	RPDPhiRes_fhaas = event.get_RpdPhiRes_fhaas();
	RPDDeltaPhiAbs = std::fabs(RPDDeltaPhi);
	RPDDeltaPhiAbs_fhaas = std::fabs(RPDDeltaPhi_fhaas);
*/


	for(unsigned int i = 0; i < _plots.size(); ++i) {
		_plots.at(i).fill(cutmask);
	}

}

void antok::Plotter::save(TDirectory* dir) {

	dir->cd();
	for(unsigned int i = 0; i < _plots.size(); ++i) {

		antok::Plot plot = _plots.at(i);
		const char* name = plot.get_template()->GetName();
		TDirectory* new_dir = dir->mkdir(name);
		assert(new_dir != 0);
		new_dir->cd();
		std::vector<TH1*> histograms = plot.get_histograms();
		for(unsigned int j = 0; j < histograms.size(); ++j) {
			TH1* hist = histograms.at(j);
			hist->Write();
		}

	}

}
