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

		std::string getShortName() const { return _shortname; };
		std::string getLongName() const { return _longname; };
		std::string getAbbreviation() const { return _abbreviation; };

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
			         std::vector<bool*> results,
					 int mode)
				: Cut(shortname, longname, abbreviation, outAddr),
				  _cuts(cuts),
				  _results(results),
				  _mode(mode) { }

			bool operator () () {
				bool retval;
				switch(_mode) {
					case 0:
						// and
						retval = true;
						for(unsigned int i = 0; i < _cuts.size(); ++i) {
							(*_cuts[i])();
							retval = retval and (*_results[i]);
						}
						(*_outAddr) = retval;
						return true;
					case 1:
						// or
						retval = false;
						for(unsigned int i = 0; i < _cuts.size(); ++i) {
							(*_cuts[i])();
							retval = retval or (*_results[i]);
						}
						(*_outAddr) = retval;
						return true;
				}
				return false;
			}

		  private:

			std::vector<antok::Cut*> _cuts;
			std::vector<bool*> _results;
			int _mode;

		};

		class NoCut : public Cut {

		public:

			NoCut(const std::string& shortname,
			      const std::string& longname,
			      const std::string& abbreviation,
			      bool* outAddr)
				: Cut(shortname, longname, abbreviation, outAddr) { }

			bool operator () () {
				(*_outAddr) = true;
				return true;
			}

		};

	}

}

#endif

