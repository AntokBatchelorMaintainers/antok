#ifndef ANTOK_CUT_H
#define ANTOK_CUT_H

#include<string>
#include<sstream>
#include<iostream>

#include<event.h>

namespace antok {

	class Cut {

	  public:

		Cut(const std::string& shortname, const std::string& longname, const std::string& abbreviation, bool* outAddr)
			: _shortname(shortname),
			  _longname(longname),
			  _abbreviation(abbreviation),
			  _outAddr(outAddr) { };

		virtual ~Cut() { delete _outAddr; }

		virtual bool operator() () = 0;
		virtual bool operator==(const Cut& rhs) = 0;

		std::string getShortName() const { return _shortname; }
		std::string getLongName() const { return _longname; }
		std::string getAbbreviation() const { return _abbreviation; }

	  protected:

		std::string _shortname;
		std::string _longname;
		std::string _abbreviation;

		bool* _outAddr;

	};

	namespace cuts {

		template< typename T>
		class RangeCut : public Cut {

		  public:

			RangeCut(const std::string& shortname,
			         const std::string& longname,
			         const std::string& abbreviation,
			         bool* outAddr,
			         T* lowerBoundOrCentralValueAddr,
			         T* upperBoundOrRangeAddr,
			         T* valueAddr,
			         int mode)
				: Cut(shortname, longname, abbreviation, outAddr),
				  _lowerBoundOrCentralValueAddr(lowerBoundOrCentralValueAddr),
				  _upperBoundOrRangeAddr(upperBoundOrRangeAddr),
				  _valueAddr(valueAddr),
				  _mode(mode) { }

			bool operator() () {
				switch(_mode)
				{
					case 0:
						// range exclusive
						(*_outAddr) = ((*_valueAddr) > (*_lowerBoundOrCentralValueAddr)) and ((*_valueAddr) < (*_upperBoundOrRangeAddr));
						return true;
					case 1:
						// range inclusive
						(*_outAddr) = ((*_valueAddr) >= (*_lowerBoundOrCentralValueAddr)) and ((*_valueAddr) <= (*_upperBoundOrRangeAddr));
						return true;
					case 2:
						// range open low exclusive
						(*_outAddr) = (*_valueAddr) < (*_upperBoundOrRangeAddr);
						return true;
					case 3:
						// range open low inclusive
						(*_outAddr) = (*_valueAddr) <= (*_upperBoundOrRangeAddr);
						return true;
					case 4:
						// range open high exclusive
						(*_outAddr) = (*_valueAddr) > (*_lowerBoundOrCentralValueAddr);
						return true;
					case 5:
						// range open high inclusive
						(*_outAddr) = (*_valueAddr) >= (*_lowerBoundOrCentralValueAddr);
						return true;
					case 6:
						// central value and range exclusive
						(*_outAddr) = ((*_valueAddr) > (*_lowerBoundOrCentralValueAddr - *_upperBoundOrRangeAddr)) and
						              ((*_valueAddr) < (*_lowerBoundOrCentralValueAddr + *_upperBoundOrRangeAddr));
						return true;
					case 7:
						// central value and range inclusive
						(*_outAddr) = ((*_valueAddr) >= (*_lowerBoundOrCentralValueAddr - *_upperBoundOrRangeAddr)) and
						              ((*_valueAddr) <= (*_lowerBoundOrCentralValueAddr + *_upperBoundOrRangeAddr));
						return true;
				}
				return false;
			};

			bool operator==(const Cut& arhs) {
				const RangeCut* rhs = dynamic_cast<const RangeCut*>(&arhs);
				if(not rhs) {
					return false;
				}
				if((this->_lowerBoundOrCentralValueAddr and not rhs->_lowerBoundOrCentralValueAddr) or (not this->_lowerBoundOrCentralValueAddr and rhs->_lowerBoundOrCentralValueAddr)) {
					return false;
				}
				if((this->_lowerBoundOrCentralValueAddr and rhs->_lowerBoundOrCentralValueAddr) and (*this->_lowerBoundOrCentralValueAddr != *rhs->_lowerBoundOrCentralValueAddr)) {
					return false;
				}
				if((this->_upperBoundOrRangeAddr and not rhs->_upperBoundOrRangeAddr) or (not this->_upperBoundOrRangeAddr and rhs->_upperBoundOrRangeAddr)) {
					return false;
				}
				if((this->_upperBoundOrRangeAddr and rhs->_upperBoundOrRangeAddr) and (*this->_upperBoundOrRangeAddr != *rhs->_upperBoundOrRangeAddr)) {
					return false;
				}
				return (*this->_valueAddr == *rhs->_valueAddr) and
				       (this->_mode == rhs->_mode);
			}

		  private:

			T* _lowerBoundOrCentralValueAddr;
			T* _upperBoundOrRangeAddr;
			T* _valueAddr;
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

			bool operator==(const Cut& arhs) {
				const EqualityCut* rhs = dynamic_cast<const EqualityCut*>(&arhs);
				if(not rhs) {
					return false;
				}
				return (*this->_leftAddr == *rhs->_leftAddr) and
				       (*this->_rightAddr == *rhs->_rightAddr) and
				       (this->_mode == rhs->_mode);
			}

		  private:

			T* _leftAddr;
			T* _rightAddr;
			int _mode;

		};

		class EllipticCut : public Cut {

		  public:

			EllipticCut(const std::string& shortname,
			            const std::string& longname,
			            const std::string& abbreviation,
			            bool* outAddr,
			            double* meanX,
			            double* meanY,
			            double* cutX,
			            double* cutY,
			            double* X,
			            double* Y,
			            int mode)
				: Cut(shortname, longname, abbreviation, outAddr),
				  _meanX(meanX),
				  _meanY(meanY),
				  _cutX(cutX),
				  _cutY(cutY),
				  _X(X),
				  _Y(Y),
				  _mode(mode) { }

			bool operator() () {
				switch(_mode) {
					case 0:
						// inclusive
						(*_outAddr) = (  ((*_X-*_meanX)*(*_X-*_meanX)/(*_cutX * *_cutX) +
						                  (*_Y-*_meanY)*(*_Y-*_meanY)/(*_cutY * *_cutY))  <= 1.  );
						return true;
					case 1:
						// exclusive
						(*_outAddr) = (  ((*_X-*_meanX)*(*_X-*_meanX)/(*_cutX * *_cutX) +
						                  (*_Y-*_meanY)*(*_Y-*_meanY)/(*_cutY * *_cutY))  < 1.  );
						return true;
				}
				return false;
			}

			bool operator==(const Cut& arhs) {
				const EllipticCut* rhs = dynamic_cast<const EllipticCut*>(&arhs);
				if(not rhs) {
					return false;
				}
				return (*this->_meanX == *rhs->_meanX) and
				       (*this->_meanY == *rhs->_meanY) and
				       (*this->_cutX == *rhs->_cutX) and
				       (*this->_cutY == *rhs->_cutY) and
				       (*this->_X == *rhs->_X) and
				       (*this->_Y == *rhs->_Y) and
				       (this->_mode == rhs->_mode);
			}

		  private:

			double* _meanX;
			double* _meanY;
			double* _cutX;
			double* _cutY;
			double* _X;
			double* _Y;
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

			bool operator==(const Cut& arhs) {
				const TriggerMaskCut* rhs = dynamic_cast<const TriggerMaskCut*>(&arhs);
				if(not rhs) {
					return false;
				}
				return (*this->_maskAddr == *rhs->_maskAddr) and
				       (*this->_triggerAddr == *rhs->_triggerAddr) and
				       (this->_mode == rhs->_mode);
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

			~CutGroup() {
				for(unsigned int i = 0; i < _cuts.size(); ++i) {
					delete _cuts[i];
				}
				_cuts.clear();
				_results.clear();
			}

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
					case 2:
						// nand
						retval = true;
						for(unsigned int i = 0; i < _cuts.size(); ++i) {
							(*_cuts[i])();
							retval = retval and (*_results[i]);
						}
						(*_outAddr) = not retval;
						return true;
				}
				return false;
			}

			bool operator==(const Cut& arhs) {
				const CutGroup* rhs = dynamic_cast<const CutGroup*>(&arhs);
				if(not rhs) {
					return false;
				}
				if((this->_cuts.size() != rhs->_cuts.size()) or (this->_mode != rhs->_mode)) {
					return false;
				}
				for(unsigned int i = 0; i < this->_cuts.size(); ++i) {
					if(not ( *(this->_cuts[i]) == *(rhs->_cuts[i]) )) {
						return false;
					}
				}
				return true;
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

			bool operator==(const Cut& arhs) {
				const NoCut* rhs = dynamic_cast<const NoCut*>(&arhs);
				if(not rhs) {
					return false;
				}
				return true;
			}

		};


		template<typename T>
		class IsNotNANCut: public Cut {

		  public:

			IsNotNANCut(const std::string& shortname,
			            const std::string& longname,
			            const std::string& abbreviation,
			            bool* outAddr,
			            T* varAddr)
				: Cut(shortname, longname, abbreviation, outAddr),
				  _var(*varAddr)
		  {}

			bool operator() () {
				(*_outAddr) = not std::isnan(_var);
				return true;
			}

			bool operator==(const Cut& arhs) {
				const IsNotNANCut* rhs = dynamic_cast<const IsNotNANCut*>(&arhs);
				if(not rhs) return false;
				return (this->_var== rhs->_var);
			}

		  private:
			T& _var;
		};

	}

}

#endif

