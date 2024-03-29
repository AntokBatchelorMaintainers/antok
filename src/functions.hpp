#ifndef ANTOK_FUNCTIONS_HPP
#define ANTOK_FUNCTIONS_HPP

#include <iostream>
#include <typeinfo>
#include <vector>

#include <boost/core/demangle.hpp>

#include "TVector3.h"
#include "TLorentzVector.h"

#include "basic_calcs.h"
#include "constants.h"
#include "charged_fit.h"

namespace antok {

	class Function
	{

	public:

		virtual bool operator() () = 0;
		virtual ~Function() { }

		std::string
		name() const
		{
			return boost::core::demangle(typeid(*this).name());
		}

	};

	namespace functions {

		template <typename T>
		class Sum : public Function
		{

		public:

			Sum(const std::vector<T*>& inputAddrsSummands,
			    const std::vector<T*>& inputAddrsSubtrahends,
			    T&                     out)
				: _inputAddrsSummands   (inputAddrsSummands),
				  _inputAddrsSubtrahends(inputAddrsSubtrahends),
				  _out                  (out)
			{
				if ((_inputAddrsSummands.size() < 1) and (_inputAddrsSubtrahends.size() < 1)) {
					std::cerr << "Got empty address vectors as input for a sum." << std::endl;
					throw 1;
				}
			}

			virtual ~Sum() { }

			bool
			operator() ()
			{
				_out = T();
				for (size_t i = 0; i < _inputAddrsSummands.size(); ++i) {
					_out += *(_inputAddrsSummands[i]);
				}
				for (size_t i = 0; i < _inputAddrsSubtrahends.size(); ++i) {
					_out -= *(_inputAddrsSubtrahends[i]);
				}
				return true;
			}


		private:

			const std::vector<T*> _inputAddrsSummands;     // need to copy input pointers, because vector in constructor is a temporary variable
			const std::vector<T*> _inputAddrsSubtrahends;  // need to copy input pointers, because vector in constructor is a temporary variable
			T&                    _out;

		};


		class Mass: public Function
		{

		public:

			Mass(const TLorentzVector& inputLV,
			     double&               out)
				: _inputLV(inputLV),
				  _out    (out)
			{ }

			virtual ~Mass() { }

			bool
			operator() ()
			{
				_out = _inputLV.M();
				return true;
			}

		private:

			const TLorentzVector& _inputLV;
			double&               _out;

		};


		class Masses: public Function
		{

		public:

			Masses(const std::vector<TLorentzVector>& inputLVs,
			     std::vector<double>&               out)
				: _inputLVs(inputLVs),
				  _out     (out)
			{ }

			virtual ~Masses() { }

			bool
			operator() ()
			{
				_out.clear();
				_out.reserve(_inputLVs.size());
				for (auto& entry : _inputLVs) {
					_out.push_back(entry.M());
				}
				return true;
			}

		private:

			const std::vector<TLorentzVector>& _inputLVs;
			std::vector<double>&               _out;

		};

		class Mass2: public Function
		{

		public:

			Mass2(const TLorentzVector& inputLV,
			      double&               out)
				: _inputLV(inputLV),
				  _out    (out)
			{ }

			virtual ~Mass2() { }

			bool
			operator() ()
			{
				_out = _inputLV.M2();
				return true;
			}

		private:

			const TLorentzVector& _inputLV;
			double&               _out;

		};


		class GetLorentzVec: public Function
		{

		public:

			enum lorentzVecDefType {XYZM = 0, PxPyPzE = 1, Vec3M = 2, Vec3E = 3};

			GetLorentzVec(const double&   x,
			              const double&   y,
			              const double&   z,
			              const double&   fourthComp,  // m or E
			              TLorentzVector& out,
			              const int       defType)
				: _xAddr         (&x),
				  _yAddr         (&y),
				  _zAddr         (&z),
				  _vec3Addr      (nullptr),
				  _fourthCompAddr(&fourthComp),
				  _out           (out),
				  _defType       (defType)
			{ }

			GetLorentzVec(const TVector3& vec3,
			              const double&   fourthComp,  // m or E
			              TLorentzVector& out,
			              const int       defType)
				: _xAddr         (nullptr),
				  _yAddr         (nullptr),
				  _zAddr         (nullptr),
				  _vec3Addr      (&vec3),
				  _fourthCompAddr(&fourthComp),
				  _out           (out),
				  _defType       (defType)
			{ }

			virtual ~GetLorentzVec() { }

			bool
			operator() ()
			{
				switch(_defType) {
					case XYZM:
						_out.SetXYZM(*_xAddr, *_yAddr, *_zAddr, *_fourthCompAddr);
						break;
					case PxPyPzE:
						_out.SetPxPyPzE(*_xAddr, *_yAddr, *_zAddr, *_fourthCompAddr);
						break;
					case Vec3M:
						_out.SetVectM(*_vec3Addr, *_fourthCompAddr);
						break;
					case Vec3E:
						_out.SetPxPyPzE(_vec3Addr->X(), _vec3Addr->Y(), _vec3Addr->Z(), *_fourthCompAddr);
						break;
				}
				return true;
			}

		private:

			const double*   _xAddr;
			const double*   _yAddr;
			const double*   _zAddr;
			const TVector3* _vec3Addr;
			const double*   _fourthCompAddr;
			TLorentzVector& _out;
			const int       _defType;

		};


		class GetBeamLorentzVec : public Function
		{

		public:

			GetBeamLorentzVec(const double&         gradX,
			                  const double&         gradY,
			                  const TLorentzVector& xLorentzVec,
			                  TLorentzVector&       out,
			                  const double*         massBeamAddr,
			                  const double*         massTargetAddr)
				: _gradX         (gradX),
				  _gradY         (gradY),
				  _xLorentzVec   (xLorentzVec),
				  _out           (out),
				  _massBeamAddr  (massBeamAddr),
				  _massTargetAddr(massTargetAddr)
			{
				if (_massTargetAddr == nullptr) {
					std::cout << "INFO: No target mass is given to reconstruct the beam energy -> using proton mass" << std::endl;
					_massTargetAddr = &antok::Constants::protonMass();
				}
				if (_massBeamAddr == nullptr) {
					std::cout << "INFO: No beam mass is given to reconstruct the beam energy -> using pion mass" << std::endl;
					_massBeamAddr= &antok::Constants::chargedPionMass();
				}
			}

			virtual ~GetBeamLorentzVec() { }

			bool
			operator() ()
			{
				const TVector3 p3Beam(_gradX, _gradY, 1.0);
				_out = antok::getBeamEnergy(p3Beam, _xLorentzVec, *_massBeamAddr, *_massTargetAddr);
				return true;
			}

		private:

			const double&         _gradX;
			const double&         _gradY;
			const TLorentzVector& _xLorentzVec;
			TLorentzVector&       _out;
			const double*         _massBeamAddr;
			const double*         _massTargetAddr;

		};


		class GetTs : public Function
		{

		public:

			GetTs(const TLorentzVector& xLorentzVec,
			      const TLorentzVector& beamLorentzVec,
			      double&               t,
			      double&               tMin,
			      double&               tPrime,
			      const double*         targetMassAddr)
				: _xLorentzVec   (xLorentzVec),
				  _beamLorentzVec(beamLorentzVec),
				  _t             (t),
				  _tMin          (tMin),
				  _tPrime        (tPrime),
				  _targetMassAddr(targetMassAddr)
			{
				if (_targetMassAddr == nullptr) {
					std::cout << "INFO: No target mass given to calculate t'. Using approximation." << std::endl;
				}
			}

			virtual ~GetTs() { }

			bool
			operator() ()
			{
				_t = (_beamLorentzVec - _xLorentzVec).Mag2();
				if (_targetMassAddr != nullptr) {
					const double EB_Lab = _beamLorentzVec.E();
					// calculate center-of-mass frame quantities to calculate tmin
					const double mT2   = (*_targetMassAddr) * (*_targetMassAddr);  // target mass squared
					const double mX2   = _xLorentzVec.M2();  // X-state mass squared
					const double mB2   = _beamLorentzVec.M2();  // beam mass squared
					const double s     = mB2 + mT2 + 2.0 * EB_Lab * std::sqrt(mT2);
					const double EB_CM = (s - mT2 + mB2) / (2.0 * std::sqrt(s));
					const double EX_CM = (s + mX2 - mT2) / (2.0 * std::sqrt(s));
					const double pB_CM = std::sqrt(EB_CM * EB_CM - mB2);
					const double pX_CM = std::sqrt(EX_CM * EX_CM - mX2);
					_tMin = -mB2 - mX2 + 2.0 * (EB_CM * EX_CM - pB_CM * pX_CM);
				} else {
					// use approximation
					_tMin = std::fabs((std::pow(_xLorentzVec.M2() - _beamLorentzVec.M2(), 2)) / (4.0 * _beamLorentzVec.Vect().Mag2()));
				}
				_tPrime = std::fabs(_t) - _tMin;
				return true;
			}

		private:

			const TLorentzVector& _xLorentzVec;
			const TLorentzVector& _beamLorentzVec;
			double&               _t;
			double&               _tMin;
			double&               _tPrime;
			const double*         _targetMassAddr;

		};


		class Sum2 : public Function
		{

		public:

			Sum2(const std::vector<double*>& in,
			     double&                     out)
				: _in (in),
				  _out(out)
			{ }

			virtual ~Sum2() { }

			bool
			operator() ()
			{
				_out = 0.0;
				for (size_t i = 0; i < _in.size(); ++i) {
					const double val = *(_in[i]);
					_out += val * val;
					;
				}
				_out = std::sqrt(_out);
				return true;
			}

		private:

			const std::vector<double*> _in;  // need to copy input pointers, because vector in constructor is a temporary variable
			double&                    _out;

		};


		template <typename T>
		class Diff : public Function
		{

		public:

			Diff(const T& in1,
			     const T& in2,
			     T&       out)
				: _in1(in1),
				  _in2(in2),
				  _out(out)
			{ }

			virtual ~Diff() { }

			bool
			operator() ()
			{
				_out = _in1 - _in2;
				return true;
			}

		private:

			const T& _in1;
			const T& _in2;
			T&       _out;
		};


		template <>
		class Diff<std::vector<double>> : public Function
		{

		public:

			Diff(const std::vector<double>& in1,
			     const std::vector<double>& in2,
			     std::vector<double>&       out)
				: _in1(in1),
				  _in2(in2),
				  _out(out)
			{ }

			virtual ~Diff() { }

			bool
			operator() ()
			{
				const size_t sizeVec = _in1.size();
				if (_in2.size() != sizeVec) {
					return false;
				}
				_out.resize(sizeVec);
				for (size_t i = 0; i < sizeVec; ++i) {
					_out[i] = _in1[i] - _in2[i];
				}
				return true;
			}

		private:

			const std::vector<double>& _in1;
			const std::vector<double>& _in2;
			std::vector<double>&       _out;

		};


		template <typename T>
		class Quotient: public Function
		{

		public:

			Quotient(const T& in1,
			         const T& in2,
			         T&       out)
				: _in1(in1),
				  _in2(in2),
				  _out(out)
			{ }

			virtual ~Quotient() { }

			bool
			operator() ()
			{
				_out = _in1 / _in2;
				return true;
			}

		private:

			const T& _in1;
			const T& _in2;
			T&       _out;

		};


		template <>
		class Quotient<std::vector<double>>: public Function
		{

		public:

			Quotient(const std::vector<double>& in1,
			         const std::vector<double>& in2,
			         std::vector<double>&       out)
				: _in1(in1),
				  _in2(in2),
				  _out(out)
			{ }

			virtual ~Quotient() { }

			bool
			operator() ()
			{
				const size_t sizeVec = _in1.size();
				if (_in2.size() != sizeVec) {
					return false;
				}
				_out.resize(sizeVec);
				for (size_t i = 0; i < sizeVec; ++i) {
					_out[i] = _in1[i] / _in2[i];
				}
				return true;
			}

		private:

			const std::vector<double>& _in1;
			const std::vector<double>& _in2;
			std::vector<double>&       _out;

		};


		template <typename T>
		class Mul: public Function
		{

		public:

			Mul(const T& in1,
			    const T& in2,
			    T&       out)
				: _in1(in1),
				  _in2(in2),
				  _out(out)
			{ }

			virtual ~Mul() { }

			bool
			operator() ()
			{
				_out = _in1 * _in2;
				return true;
			}

		private:

			const T& _in1;
			const T& _in2;
			T&       _out;

		};


		template <>
		class Mul<std::vector<double>>: public Function
		{

		public:

			Mul(const std::vector<double>& in1,
			    const std::vector<double>& in2,
			    std::vector<double>&       out)
				: _in1(in1),
				  _in2(in2),
				  _out(out)
			{ }

			virtual ~Mul() { }

			bool
			operator() ()
			{
				const size_t sizeVec = _in1.size();
				if (_in2.size() != sizeVec) {
					return false;
				}
				_out.resize(sizeVec);
				for (size_t i = 0; i < sizeVec; ++i) {
					_out[i] = _in1[i] * _in2[i];
				}
				return true;
			}

		private:

			const std::vector<double>& _in1;
			const std::vector<double>& _in2;
			std::vector<double>&       _out;

		};


		template <typename T>
		class Abs : public Function
		{

		public:

			Abs(const T& in,
			    double&  out)
				: _in (in),
				  _out(out)
			{ }

			virtual ~Abs() { }

			bool
			operator() ()
			{
				_out = __abs(_in);
				return true;
			}

		private:

			const T& _in;
			double&  _out;

			double __abs (const TVector2& x) { return x.Mod(); }
			double __abs (const TVector3& x) { return x.Mag(); }
			template <typename Q>
			double __abs(const Q& x) { return std::fabs(x); }

		};


		template <typename T>
		class Abs<std::vector<T>>: public Function
		{

		public:

			Abs(const std::vector<T>& in,
			    std::vector<double>&  out)
				: _in(in),
				  _out(out)
			{ }

			virtual ~Abs() { }

			bool
			operator() ()
			{
				const size_t sizeVec = _in.size();
				_out.resize(sizeVec);
				for (size_t i = 0; i < sizeVec; ++i) {
				    _out[i] = __abs(_in[i]);
				}
				return true;
			}

		private:

			const std::vector<T>& _in;
			std::vector<double>&  _out;

			double __abs (const TVector2& x) { return x.Mod(); }
			double __abs (const TVector3& x) { return x.Mag(); }
			template <typename Q>
			double __abs(const Q& x) { return std::fabs(x); }

		};


		template <typename T>
		class Log : public Function
		{

		public:

			Log(const T&      in,
			    const double& base,
			    double&       out)
				: _in  (in),
				  _base(base),
				  _out (out)
			{ }

			virtual ~Log() { }

			bool
			operator() ()
			{
				if (std::isnan(_base)) {
					// default: natural logarithm
					_out = std::log(_in);
					return true;
				} else {
					_out = std::log(_in) / std::log(_base);
					return true;
				}
			}

		private:

			const T&     _in;
			const double _base;  // constant parameter, needs to be copied; if nan, ln is calculated
			double&      _out;

		};


		template <>
		class Log<std::vector<double>>: public Function
		{

		public:

			Log(const std::vector<double>& in,
			    const double&              base,
			    std::vector<double>&       out)
				: _in  (in),
				  _base(base),
				  _out (out)
			{ }

			virtual ~Log() { }

			bool
			operator() ()
			{
				const size_t sizeVec = _in.size();
				_out.resize(sizeVec);
				for (size_t i = 0; i < sizeVec; ++i) {
					if (std::isnan(_base)) {
						// default: natural logarithm
						_out[i] = std::log(_in[i]);
					} else {
						_out[i] = std::log(_in[i]) / std::log(_base);
					}
				}
				return true;
			}

		private:

			const std::vector<double>& _in;
			const double               _base;  // constant parameter, needs to be copied
			std::vector<double>&       _out;

		};


		template <typename T>
		class Sqrt : public Function
		{

		public:

			Sqrt(const T& in,
			     double&  out)
				: _in (in),
				  _out(out)
			{ }

			virtual ~Sqrt() { }

			bool
			operator() ()
			{
				_out = std::sqrt(_in);
				return true;
			}

		private:

			const T& _in;
			double&  _out;

		};


		template <>
		class Sqrt<std::vector<double>>: public Function
		{

		public:

			Sqrt(const std::vector<double>& in,
			     std::vector<double>&       out)
				: _in(in),
				  _out(out)
			{ }

			virtual ~Sqrt() { }

			bool
			operator() ()
			{
				const size_t sizeVec = _in.size();
				_out.resize(sizeVec);
				for (size_t i = 0; i < sizeVec; ++i) {
					_out[i] = std::sqrt(_in[i]);
				}
				return true;
			}

		private:

			const std::vector<double>& _in;
			std::vector<double>&       _out;

		};


		class Energy : public Function
		{

		public:

			Energy(const TLorentzVector& in,
			       double&               out)
				: _in (in),
				  _out(out)
			{ }

			virtual ~Energy() { }

			bool
			operator() ()
			{
				_out = _in.E();
				return true;
			}

		private:

			const TLorentzVector& _in;
			double&               _out;

		};


		class Energies : public Function
		{

		public:

			Energies(const std::vector<TLorentzVector>& in,
			         std::vector<double>&               out)
				: _in (in),
				  _out(out)
			{ }

			virtual ~Energies() { }

			bool
			operator() ()
			{
				_out.clear();
				_out.reserve(_in.size());
				for (auto& entry : _in) {
					_out.push_back(entry.E());
				}
				return true;
			}

		private:

			const std::vector<TLorentzVector>& _in;
			std::vector<double>&               _out;

		};


		class RadToDegree : public Function
		{

		public:

			RadToDegree(const double& in,
			            double&       out)
				: _in (in),
				  _out(out)
			{ }

			virtual ~RadToDegree() { }

			bool
			operator() ()
			{
				_out = _in * TMath::RadToDeg();
				return true;
			}

		private:

			const double& _in;
			double&       _out;

		};


		class ConvertIntToDouble : public Function
		{

		public:

			ConvertIntToDouble(const int& in,
			                   double&    out)
				: _in (in),
				  _out(out)
			{ }

			virtual ~ConvertIntToDouble() { }

			bool
			operator() ()
			{
				_out = _in;
				return true;
			}

		private:

			const int& _in;
			double&    _out;

		};


		class GetGradXGradY : public Function
		{

		public:

			GetGradXGradY(const TLorentzVector& lorentzVec,
			              double&               xGrad,
			              double&               yGrad)
				: _lorentzVec(lorentzVec),
				  _xGrad     (xGrad),
				  _yGrad     (yGrad)
			{ }

			virtual ~GetGradXGradY() { }

			bool
			operator() ()
			{
				const TVector3 vect = _lorentzVec.Vect();
				const double   z    = vect.Z();
				_xGrad = vect.X() / z;
				_yGrad = vect.Y() / z;
				return true;
			}

		private:

			const TLorentzVector& _lorentzVec;
			double&               _xGrad;
			double&               _yGrad;

		};


		class GetLorentzVectorAttributes : public Function
		{

		public:

			GetLorentzVectorAttributes(const TLorentzVector& lorentzVec,
			                           double&               x,
			                           double&               y,
			                           double&               z,
			                           double&               phi,
			                           double&               theta)
				: _lorentzVec(lorentzVec),
				  _x         (x),
				  _y         (y),
				  _z         (z),
				  _phi       (phi),
				  _theta     (theta)
			{ }

			virtual ~GetLorentzVectorAttributes() { }

			bool
			operator() ()
			{
				_x     = _lorentzVec.X();
				_y     = _lorentzVec.Y();
				_z     = _lorentzVec.Z();
				_phi   = _lorentzVec.Phi();
				_theta = _lorentzVec.Theta();
				return true;
			}

		private:

			const TLorentzVector& _lorentzVec;
			double&               _x;
			double&               _y;
			double&               _z;
			double&               _phi;
			double&               _theta;

		};


		class GetVectorTVector2 : public Function
		{

		public:

			GetVectorTVector2(const std::vector<double>& x,
			                  const std::vector<double>& y,
			                  std::vector<TVector2>&     v)
				: _x(x),
				  _y(y),
				  _v(v)
			{ }

			virtual ~GetVectorTVector2() { }

			bool
			operator() ()
			{
				if (_x.size() != _y.size()) {
					return false;
				}
				_v.resize(_x.size());
				for (size_t i = 0; i < _x.size(); ++i) {
					_v[i] = TVector2(_x[i], _y[i]);
				}
				return true;
			}

		private:

			const std::vector<double>& _x;
			const std::vector<double>& _y;
			std::vector<TVector2>&     _v;

		};


		class GetTVector2 : public Function
		{

		public:

			GetTVector2(const double& x,
			            const double& y,
			            TVector2&     out)
				: _x  (x),
				  _y  (y),
				  _out(out)
			{ }

			virtual ~GetTVector2() { }

			bool
			operator() ()
			{
				_out = TVector2(_x, _y);
				return true;
			}

		private:

			const double& _x;
			const double& _y;
			TVector2&     _out;

		};


		class GetVectorTVector3 : public Function
		{

		public:

			GetVectorTVector3(const std::vector<double>& x,
			                  const std::vector<double>& y,
			                  const std::vector<double>& z,
			                  std::vector<TVector3>&     v)
				: _x(x),
				  _y(y),
				  _z(z),
				  _v(v)
			{ }

			virtual ~GetVectorTVector3() { }

			bool
			operator() ()
			{
				if ((_x.size() != _y.size()) and (_x.size() != _z.size())) {
					return false;
				}
				_v.resize(_x.size());
				for (size_t i = 0; i < _x.size(); ++i) {
					_v[i] = TVector3(_x[i], _y[i], _z[i]);
				}
				return true;
			}

		private:

			const std::vector<double>& _x;
			const std::vector<double>& _y;
			const std::vector<double>& _z;
			std::vector<TVector3>&     _v;

		};


		class GetTVector3 : public Function
		{

		public:

			GetTVector3(const double& x,
			            const double& y,
			            const double& z,
			            TVector3&     out)
				: _x  (x),
				  _y  (y),
				  _z  (z),
				  _out(out)
			{ }

			virtual ~GetTVector3() { }

			bool
			operator() ()
			{
				_out.SetXYZ(_x, _y, _z);
				return true;
			}

		private:

			const double& _x;
			const double& _y;
			const double& _z;
			TVector3&     _out;

		};


		template <typename T>
		class GetVectorEntry : public Function
		{

		public:

			GetVectorEntry(const std::vector<T>& vector,
			               const int&            entry,
			               T&                    result)
				: _vector(vector),
				  _entry (entry),
				  _result(result)
			{
				if (_entry < 0) {
					std::cerr << "Got invalid entry position for a vector entry." << std::endl;
					throw 1;
				}
			};

			virtual ~GetVectorEntry() { }

			bool
			operator() ()
			{
				if (_vector.size() == 0) {
					_result = T();
					return true;
				}
				if (_entry >= (int)_vector.size()) {
					return false;
				}
				_result = _vector[_entry];
				return true;
			};

		private:

			const std::vector<T>& _vector;
			const int             _entry;  // constant parameter, needs to be copied
			T&                    _result;

		};


		template <typename T>
		class GetVectorSize : public Function
		{

		public:

			GetVectorSize(const std::vector<T>& vector,
			              int&                  result)
				: _vector(vector),
				  _result(result)
			{ };

			virtual ~GetVectorSize() { }

			bool
			operator() ()
			{
				_result = _vector.size();
				return true;
			};

		private:

			const std::vector<T>& _vector;
			int&                  _result;

		};


		class GetTVector3FromTLorenzVector : public Function
		{

		public:

			GetTVector3FromTLorenzVector(const TLorentzVector& lv,
			                             TVector3&             v)
				: _lv(lv),
				  _v (v)
			{ }

			virtual ~GetTVector3FromTLorenzVector() { }

			bool
			operator() ()
			{
				_v = _lv.Vect();
				return true;
			}

		private:

			const TLorentzVector& _lv;
			TVector3&             _v;

		};

		// Performs the kinematic fit of two charged particles:
		class GetChargedKinematicFitting : public Function
		{

		public:

			GetChargedKinematicFitting(const TVector3& Particle1Momentum, // momentum of track of particle 1, extrapolated to Z of SV
									   const TVector3& Particle2Momentum,
									   const double&   Particle1Energy, // sqrt(p_track*p_track + m*m)
									   const double&   Particle2Energy,
									   const std::vector<double>& Particle1MomentumCovariance,  // covariance of track 4-momentum as vector of len 16=(4x4)
									   const std::vector<double>& Particle2MomentumCovariance,
									   const double    Particle1Mass,
									   const double    Particle2Mass,
									   const double    Mass, // mass to which the invariant mass of particles 1 and 2 is fitted
									   const double    PrecisionGoal,             // Precision goal relative to Mass**2
									   TLorentzVector& ResultLorentzVector1,      // Lorentz vectors of particle 1 after fit
									   TLorentzVector& ResultLorentzVector2,      // Lorentz vectors of particle 2 after fit
									   double&         ResultChi2,                // chi^2 values of fit
									   double&         ResultPValue,              // P-values of fit
									   int&            ResultNmbIterations,       // number of iterations required to reach PrecisionGoal
									   int&            ResultSuccess,             // indicates whether fit was successful
									   double&         ResultPullsX0,             // pulls for x direction of first particle in pair
									   double&         ResultPullsY0,             // pulls for y direction of first particle in pair
									   double&         ResultPullsE0,             // pulls for energy of first particle in pair
									   double&         ResultPullsX1,             // pulls for x direction of second particle in pair
									   double&         ResultPullsY1,             // pulls for y direction of second particle in pair
									   double&         ResultPullsE1,             // pulls for energy of second particle in pair
									   std::vector<double>& ResultTransfCov)      // Transformed covariance matrix for x1,y1,E1,x2,y2,E2
				:   _particle1Momentum          (Particle1Momentum),
					_particle2Momentum          (Particle2Momentum),
					_particle1Energy            (Particle1Energy),
					_particle2Energy            (Particle2Energy),
					_particle1Mass              (Particle1Mass),
					_particle2Mass              (Particle2Mass),
					_particle1MomentumCovariance(Particle1MomentumCovariance),
					_particle2MomentumCovariance(Particle2MomentumCovariance),
					_Mass                       (Mass),
					_PrecisionGoal              (PrecisionGoal),
					_ResultLorentzVector1       (ResultLorentzVector1),
					_ResultLorentzVector2       (ResultLorentzVector2),
					_ResultChi2                 (ResultChi2),
					_ResultPValue               (ResultPValue),
					_ResultNmbIterations        (ResultNmbIterations),
					_ResultSuccess              (ResultSuccess),
					_ResultPullsX0              (ResultPullsX0),
					_ResultPullsY0              (ResultPullsY0),
					_ResultPullsE0              (ResultPullsE0),
					_ResultPullsX1              (ResultPullsX1),
					_ResultPullsY1              (ResultPullsY1),
					_ResultPullsE1              (ResultPullsE1),
					_ResultTransfCov            (ResultTransfCov)
			{ }

			virtual ~GetChargedKinematicFitting() { }

			bool
			operator() ()
			{

				_ResultLorentzVector1.SetXYZT(0.,0.,0.,0.);
				_ResultLorentzVector2.SetXYZT(0.,0.,0.,0.);
				_ResultChi2 = 0.;
				_ResultPValue = 0.;
				_ResultNmbIterations = 0;
				_ResultPullsX0 = 0.;
				_ResultPullsY0 = 0.;
				_ResultPullsE0 = 0.;
				_ResultPullsX1 = 0.;
				_ResultPullsY1 = 0.;
				_ResultPullsE1 = 0.;
				_ResultSuccess = 0;

				// Fit particle momenta to given _Mass
				bool success = false;

				antok::ChargedFit chargedFit(
					_particle1Momentum,
					_particle2Momentum,
					_particle1Energy,
					_particle2Energy,
					_particle1Mass,
					_particle2Mass,
					_particle1MomentumCovariance,
					_particle2MomentumCovariance,
					_Mass,
					_PrecisionGoal);
				success = chargedFit.doFit();
				if (success) {
					_ResultLorentzVector1 = chargedFit.getImprovedLV1();
					_ResultLorentzVector2 = chargedFit.getImprovedLV2();
					_ResultPullsX0        = chargedFit.pullValues()[0];
					_ResultPullsY0        = chargedFit.pullValues()[1];
					_ResultPullsE0        = chargedFit.pullValues()[2];
					_ResultPullsX1        = chargedFit.pullValues()[3];
					_ResultPullsY1        = chargedFit.pullValues()[4];
					_ResultPullsE1        = chargedFit.pullValues()[5];
					_ResultChi2           = chargedFit.chi2Value();
					_ResultPValue         = chargedFit.pValue();
					_ResultNmbIterations  = chargedFit.nmbIterations();
					_ResultTransfCov      = chargedFit.getTransfCov();
				}

				if (success) {
					_ResultSuccess = 1;
				}
				return true;
			}

		private:

			const TVector3& _particle1Momentum;
			const TVector3& _particle2Momentum;
			const double&   _particle1Energy;
			const double&   _particle2Energy;
			const double    _particle1Mass;
			const double    _particle2Mass;
			const std::vector<double>& _particle1MomentumCovariance;
			const std::vector<double>& _particle2MomentumCovariance;
			const double    _Mass;
			const double    _PrecisionGoal;
			TLorentzVector& _ResultLorentzVector1;      // Lorentz vectors of particle 1 after fit
			TLorentzVector& _ResultLorentzVector2;      // Lorentz vectors of particle 2 after fit
			double&         _ResultChi2;               // chi^2 values of fits
			double&         _ResultPValue;             // P-values of fits
			int&            _ResultNmbIterations;       // number of iterations required to reach PrecisionGoal
			int&            _ResultSuccess;             // indicates whether fit was successful
			double&         _ResultPullsX0;             // pulls for x direction of first particle in pairs
			double&         _ResultPullsY0;             // pulls for y direction of first particle in pairs
			double&         _ResultPullsE0;             // pulls for energy of first particle in pairs
			double&         _ResultPullsX1;             // pulls for x direction of second particle in pairs
			double&         _ResultPullsY1;             // pulls for y direction of second particle in pairs
			double&         _ResultPullsE1;
			std::vector<double>& _ResultTransfCov;


		};

	}  // functions namespace

}  // antok namespace

#endif  // ANTOK_FUNCTIONS_HPP
