
#include<plotter.h>

hlib::Plotter::Plotter() {

	TH1D* mass_5pi = new TH1D("mass_5pi", "mass_5Pi", 500, 0, 7);

	TH1D* mom_5pi = new TH1D("mom_5pi", "mom_5Pi", 500, 0, 250);

	TH1D* mom_5pi_raw = new TH1D("mom_5pi_raw", "mom_5Pi_raw", 500, 0, 250);

	TH1D* calc_beam_E = new TH1D("calc_beam_E", "calc_beam_E", 500, 0, 250);

	TH1D* rpd_mult = new TH1D("rpd_mult", "rpd_mult", 10, 0, 10);

	TH1D* rpd_Pxh = new TH1D("rpd_Px", "rpd_Px", 1000, 0, 10);

	TH2D* vtx_pos = new TH2D("vtx_pos", "vtx_pos", 1000, -5, 5, 1000, -5, 5);

	TH1D* vtx_zh = new TH1D("vtx_z", "vtx_z", 2000, -200, 200);

	TH1D* t_primh = new TH1D("t_prime", "t_prime", 1000, -5, 5);

	TH1D* delta_phih = new TH1D("delta_phi", "delta_phi", 500, -7, 7);

	TH1D* trig_maskh = new TH1D("trigger_mask", "trigger_mask", 15, 0, 15);



};
