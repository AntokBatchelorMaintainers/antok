
#include<assert.h>
#include<sstream>

#include<cutter.h>
#include<plot.h>

hlib::Plot::Plot(std::vector<int> cutmasks, TH1* hist_template, int dim) :
		_cutmasks(cutmasks),
		_dim(dim)
{

	assert(hist_template != NULL);
	if(_cutmasks.size() == 0) {
		_cutmasks.push_back(0);
	}
	hlib::Cutter* cutter = hlib::Cutter::instance();
	unsigned int no_cuts = cutter->get_no_cuts();
	std::stringstream sstr;
	bool found_zero = false;
	for(unsigned int i = 0; i < _cutmasks.size(); ++i) {
		if(_cutmasks.at(i) == 0) {
			found_zero = true;
			break;
		}
	}
	if(!found_zero) {
		_cutmasks.push_back(0);
	}
	for(unsigned int i = 0; i < _cutmasks.size(); ++i) {
		int cutmask = _cutmasks.at(i);
		sstr<<hist_template->GetName()<<"_";
		for(unsigned int j = no_cuts-1; j > 0; --j) {
			if((cutmask>>j)&0x1) {
				sstr<<"1";
			} else {
				sstr<<"0";
			}
		}
		sstr.str("");
		TH1* new_hist = dynamic_cast<TH1*>(hist_template->Clone(sstr.str().c_str()));
		assert(new_hist != NULL);
		sstr<<hist_template->GetTitle();
		sstr<<" "<<cutter->get_abbreviations(cutmask);
		new_hist->SetTitle(sstr.str().c_str());
		sstr.str("");
		_histograms.push_back(new_hist);
	}

};

void hlib::Plot::fill(int cutmask, double data1, double data2) {

	for(unsigned int i = 0; i < _cutmasks.size(); ++i) {
		if((~_cutmasks.at(i)&~cutmask) == _cutmasks.at(i)) {
			if(_dim == 1) {
				assert(data2 == -999999999.);
				_histograms.at(i)->Fill(data1);
			} else {
				assert(data2 != -999999999.);
				_histograms.at(i)->Fill(data1, data2);
			}
		}
	}

}

std::vector<TH1*> hlib::Plot::get_histograms() {

	std::vector<TH1*> hists;

	for(unsigned int i = 0; i < _cutmasks.size(); ++i) {

	};

	return hists;

};
