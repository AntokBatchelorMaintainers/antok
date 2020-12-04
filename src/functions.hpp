#ifndef ANTOK_FUNCTIONS_HPP
#define ANTOK_FUNCTIONS_HPP

#include <iostream>
#include <typeinfo>
#include <vector>

#include <boost/core/demangle.hpp>

#include "TLorentzVector.h"

#include "basic_calcs.h"
#include "constants.h"

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

			double __abs (const TVector3& x) { return x.Mag(); }
			template <typename Q, class Dummy = int>
			double __abs(const Q& x) { return std::fabs(x); }

		};


		template <>
		class Abs<std::vector<double>>: public Function
		{

		public:

			Abs(const std::vector<double>& in,
			    std::vector<double>&       out)
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
				    _out[i] = std::fabs(_in[i]);
				}
				return true;
			}

		private:

			const std::vector<double>& _in;
			std::vector<double>&       _out;

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

	}  // functions namespace

}  // antok namespace

#endif  // ANTOK_FUNCTIONS_HPP
