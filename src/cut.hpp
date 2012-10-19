#ifndef ANTOK_CUT_H
#define ANTOK_CUT_H

#include<string>
#include<sstream>

#include<event.h>

namespace antok {

	class Cut {

	  public:

		Cut(const std::string& shortname, const std::string& longname, const std::string& abbreviation, bool* outAddr)
			: _shortname(shortname),
			  _longname(longname),
			  _abbreviation(abbreviation),
			  _outAddr(outAddr) { };

		virtual ~Cut() { };

		virtual bool operator() () = 0;

		std::string get_shortname() const { return _shortname; };
		std::string get_longname() const { return _longname; };
		std::string get_abbreviation() const { return _abbreviation; };

	  protected:

		std::string _shortname;
		std::string _longname;
		std::string _abbreviation;

		bool* _outAddr;

	};

	namespace cuts {

		class RangeCut : public Cut {

		  public:

			RangeCut(const std::string& shortname,
			         const std::string& longname,
			         const std::string& abbreviation,
			         bool* outAddr,
			         double* lowerBoundAddr,
			         double* upperBoundAddr,
			         double* valueAddr,
			         int mode)
				: Cut(shortname, longname, abbreviation, outAddr),
				  _lowerBoundAddr(lowerBoundAddr),
				  _upperBoundAddr(upperBoundAddr),
				  _valueAddr(valueAddr),
				  _mode(mode) { }

			bool operator() () {
				switch(_mode)
				{
					case 0:
						// range exclusive
						(*_outAddr) = ((*_valueAddr) > (*_lowerBoundAddr)) and ((*_valueAddr) < (*_upperBoundAddr));
						return true;
					case 1:
						// range inclusive
						(*_outAddr) = ((*_valueAddr) >= (*_lowerBoundAddr)) and ((*_valueAddr) <= (*_upperBoundAddr));
						return true;
					case 2:
						// range open low exclusive
						(*_outAddr) = (*_valueAddr) < (*_upperBoundAddr);
						return true;
					case 3:
						// range open low inclusive
						(*_outAddr) = (*_valueAddr) <= (*_upperBoundAddr);
						return true;
					case 4:
						// range open high exclusive
						(*_outAddr) = (*_valueAddr) > (*_lowerBoundAddr);
						return true;
					case 5:
						// range open high inclusive
						(*_outAddr) = (*_valueAddr) >= (*_lowerBoundAddr);
						return true;
				}
				return false;
			};

		  private:

			double* _lowerBoundAddr;
			double* _upperBoundAddr;
			double* _valueAddr;
			int _mode;

		};

		template<typename T>
		class EqualityCut : public Cut {

		  public:

			EqualityCut(const std::string& shortname,
			            const std::string& longname,
			            const std::string& abbreviation,
			            bool* outAddr,
			            T* leftAddr,
			            T* rightAddr,
			            int mode)
				: Cut(shortname, longname, abbreviation, outAddr),
				  _leftAddr(leftAddr),
				  _rightAddr(rightAddr),
				  _mode(mode) { }

			bool operator() () {
				switch(_mode) {
					case 0:
						// ==
						(*_outAddr) = (*_leftAddr) == (*_rightAddr);
						return true;
					case 1:
						// !=
						(*_outAddr) = (*_leftAddr) != (*_rightAddr);
						return true;
				}
				return false;
			}

		  private:

			T* _leftAddr;
			T* _rightAddr;
			int _mode;

		};

		class TriggerMaskCut: public Cut {

		  public:

			TriggerMaskCut(const std::string& shortname,
			               const std::string& longname,
			               const std::string& abbreviation,
			               bool* outAddr,
			               int* maskAddr,
			               int* triggerAddr,
			               int mode)
				: Cut(shortname, longname, abbreviation, outAddr),
				  _maskAddr(maskAddr),
				  _triggerAddr(triggerAddr),
				  _mode(mode) { }

			bool operator() () {
				switch(_mode) {
					case 0:
						(*_outAddr) = (*_maskAddr)&(*_triggerAddr);
						return true;
				}
				return false;
			}

		  private:

			int* _maskAddr;
			int* _triggerAddr;
			int _mode;

		};

		class CutGroup : public Cut {

		  public:

			CutGroup(const std::string& shortname,
			         const std::string& longname,
			         const std::string& abbreviation,
			         bool* outAddr,
			         std::vector<antok::Cut*> cuts,
			         std::vector<bool*> results)
				: Cut(shortname, longname, abbreviation, outAddr),
				  _cuts(cuts),
				  _results(results) { }

			bool operator () () {
				bool retval = true;
				for(unsigned int i = 0; i < _cuts.size(); ++i) {
					(*_cuts[i])();
					retval = retval and (*_results[i]);
				}
				(*_outAddr) = retval;
				return true;
			}

		  private:

			std::vector<antok::Cut*> _cuts;
			std::vector<bool*> _results;

		};

	}
/*
// ---------------------------------- OLD STUFF

	// Cut on the trigger mask.
	class TrigMask : public Cut {

	  public:
		TrigMask(int trigmask, int num = 0);
		bool operator() (const antok::Event& event) const { return (!((event.rawData->TrigMask)&0x1)); };

	  private:
		int _trigmask;

	};
	inline TrigMask::TrigMask(int trigmask, int num) : _trigmask(trigmask) {
		std::stringstream sstr;
		sstr<<"trigmask_&_"<<std::hex<<_trigmask;
		shortname = sstr.str();
		sstr.str("");
		sstr<<"Trigger Mask & "<<std::hex<<_trigmask;
		longname = sstr.str();
		abbreviation = "tm";
		if(num != 0) {
			sstr.str("");
			sstr<<abbreviation<<num;
			abbreviation = sstr.str();
		}
	}

	// Cut on the vertex Z position.
	class VrtxZ : public Cut {

	  public:
		VrtxZ(double zmin, double zmax, int num = 0);
		bool operator() (const antok::Event& event) const { return ((event.rawData->Z_primV > _zmin) ||
		                                                            (event.rawData->Z_primV < _zmax)); };

	  private:
		double _zmin;
		double _zmax;

	};
	inline VrtxZ::VrtxZ(double zmin, double zmax, int num) : _zmin(zmin), _zmax(zmax) {
		std::stringstream sstr;
		sstr<<_zmin<<"<vtx_z<"<<_zmax;
		shortname = sstr.str();
		sstr.str("");
		sstr<<"Vertex Z in ]"<<_zmin<<","<<_zmax<<"[";
		longname = sstr.str();
		abbreviation = "vz";
		if(num != 0) {
			sstr.str("");
			sstr<<abbreviation<<num;
			abbreviation = sstr.str();
		}
	};

	// Cut on the vertex R.
	class VrtxR : public Cut {

	  public:
		VrtxR(double rmax, int num = 0);
		bool operator() (const antok::Event& event) const { return (std::pow(event.rawData->X_primV, 2) +
		                                                            std::pow(event.rawData->Y_primV, 2) > _rmax2); };

	  private:
		double _rmax;
		double _rmax2;

	};
	inline VrtxR::VrtxR(double rmax, int num) : _rmax(rmax), _rmax2(rmax * rmax) {
		std::stringstream sstr;
		sstr<<"vtx_R<"<<rmax;
		shortname = sstr.str();
		sstr.str("");
		sstr<<"Vertex.R() < "<<rmax;
		longname = sstr.str();
		abbreviation = "vr";
		if(num != 0) {
			sstr.str("");
			sstr<<abbreviation<<num;
			abbreviation = sstr.str();
		}
	};

	// Cut on the number of RPD tracks.
	class nRPDTracks : public Cut {

	  public:
		nRPDTracks(int nTracks, int num = 0);
		bool operator() (const antok::Event& event) const { return (event.rawData->nbrRPDTracks != _nTracks); };

	  private:
		int _nTracks;

	};
	inline nRPDTracks::nRPDTracks(int nTracks, int num) : _nTracks(nTracks) {
		std::stringstream sstr;
		sstr<<"nRPDTrks="<<_nTracks;
		shortname = sstr.str();
		sstr.str("");
		sstr<<_nTracks<<" RPD track(s)";
		longname = sstr.str();
		abbreviation = "Rt";
		if(num != 0) {
			sstr.str("");
			sstr<<abbreviation<<num;
			abbreviation = sstr.str();
		}
	};

	// Cut on the RPD's proton mass.
	class RPDProtMass : public Cut {

	  public:
		RPDProtMass(double minMass, int num = 0);
		bool operator() (const antok::Event& event) const { return (event.get_pProton().M() < _minMass); };

	  private:
		double _minMass;

	};
	inline RPDProtMass::RPDProtMass(double minMass, int num) : _minMass(minMass) {
		std::stringstream sstr;
		sstr<<"p+_mass>"<<_minMass;
		shortname = sstr.str();
		sstr.str("");
		sstr<<"RPD's proton mass > "<<_minMass;
		longname = sstr.str();
		abbreviation = "Rm";
		if(num != 0) {
			sstr.str("");
			sstr<<abbreviation<<num;
			abbreviation = sstr.str();
		}
	};

	// CEDAR kaon cut.
	class CedarKaon : public Cut {

	  public:
		CedarKaon();
		bool operator() (const antok::Event& event) const { return (event.rawData->cedarID_bayes != 0); };

	};
	inline CedarKaon::CedarKaon() {
		shortname = "cedarID_bayes=0";
		longname = "cedarID_bayes = 0";
		abbreviation = "ik";
	};

	// T-Prime cut.
	class TPrime : public Cut {

	  public:
		TPrime(double tmin, double tmax, int num = 0);
		inline bool operator() (const antok::Event& event) const {
			if (_tmin >= 0. && _tmax >= 0.) {
				return ((event.get_tPrime() <= _tmin) || (event.get_tPrime() >= _tmax));
			} else if (_tmin < 0. && _tmax >= 0.) {
				return (event.get_tPrime() >= _tmax);
			} else if (_tmin >= 0. && _tmax < 0.) {
				return (event.get_tPrime() <= _tmin);
			} else {
				return false;
			}
		};

	  private:
		double _tmin;
		double _tmax;

	};
	inline TPrime::TPrime(double tmin, double tmax, int num) : _tmin(tmin), _tmax(tmax) {
		std::stringstream sstr;
		sstr<<_tmin<<"<t_prim<"<<_tmax;
		shortname = sstr.str();
		sstr.str("");
		sstr<<"t' in ]"<<_tmin<<","<<_tmax<<"[";
		longname = sstr.str();
		abbreviation = "tp";
		if(num != 0) {
			sstr.str("");
			sstr<<abbreviation<<num;
			abbreviation = sstr.str();
		}
	};

	// RPD planarity cut.
	class RPDPlanarity : public Cut {

	  public:
		RPDPlanarity();
		bool operator() (const antok::Event& event) const { return (std::fabs(event.get_RpdDeltaPhi()) > event.get_RpdPhiRes()); };

	};
	inline RPDPlanarity::RPDPlanarity() {
		shortname = "rpd_plan";
		longname = "RPD planarity cut";
		abbreviation = "Rp";
	};

	// Exclusivity cut.
	class Exclusivity : public Cut {

	  public:
		Exclusivity(double mean, double win, int num = 0);
		bool operator() (const antok::Event& event) const { return (std::fabs(event.get_pBeam().Energy()-_mean) >= _win); };

	  private:
		double _mean;
		double _win;

	};
	inline Exclusivity::Exclusivity(double mean, double win, int num) : _mean(mean), _win(win) {
		std::stringstream sstr;
		sstr<<(_mean-_win)<<"<E_X<"<<(_mean+_win);
		shortname = sstr.str();
		sstr.str("");
		sstr<<"Total Energy in ]"<<(_mean-_win)<<","<<(_mean+_win)<<"[";
		longname = sstr.str();
		abbreviation = "ex";
		if(num != 0) {
			sstr.str("");
			sstr<<abbreviation<<num;
			abbreviation = sstr.str();
		}
	};
*/
}

#endif

