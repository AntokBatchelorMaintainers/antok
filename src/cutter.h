#ifndef ANTOK_CUTTER_H
#define ANTOK_CUTTER_H

#include<string>
#include<vector>

#include<TH1D.h>

#include<cut.hpp>
#include<event.h>

namespace antok {

	class Cutter {

		friend class Initializer;

	  public:

		static Cutter* instance();

		int get_cutmask();

		std::string get_abbreviations(int bitmask);

		unsigned int get_no_cuts() { return _cuts.size(); };

		bool set_stats_histogram(TH1D* stats);
		TH1D* get_stats_histogram() { return _statsHist; };

	  private:

		Cutter();
		~Cutter();

		static Cutter* _cutter;

		std::vector<antok::Cut*> _cuts;

		std::vector<std::vector<antok::Cut*> > _cutTrains;

		TH1D* _statsHist;

	};

}

#endif

