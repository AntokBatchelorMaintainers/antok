#ifndef ANTOK_USER_JBECKERS_FUNCTIONS_HPP
#define ANTOK_USER_JBECKERS_FUNCTIONS_HPP

#include "TVector3.h"
#include "TLorentzVector.h"
#include "charged_fit.h"

namespace antok {
namespace user {
namespace jbeckers {
namespace functions {

	// Multiplies a variable by a constant scale (must be of type double):
	template <typename T>
	class Scale: public Function
	{

	public:

		Scale(const double& scaleFactor,
		      const T&      value,
		      T&            result)
			: _scaleFactor(scaleFactor),
			  _value      (value),
			  _result     (result)
		{ }

		virtual ~Scale() { }

		bool
		operator() ()
		{
			_result = _scaleFactor * _value;
			return true;
		}

	private:

		const double _scaleFactor;  // constant parameter, needs to be copied (-> no & !)
		const T&     _value;
		T&           _result;

	};


	template <>
	class Scale<std::vector<double>>: public Function
	{

	public:

		Scale(const double&               scaleFactor,
		      const std::vector<double>&  values,
		      std::vector<double>&        results)
			: _scaleFactor(scaleFactor),
			  _values     (values),
			  _results    (results)
		{ }

		virtual ~Scale() { }

		bool
		operator() ()
		{
			const size_t sizeVec = _values.size();
			_results.resize(sizeVec);
			for (size_t i = 0; i < sizeVec; ++i) {
				_results[i] = _scaleFactor * _values[i];
			}
			return true;
		}

	private:

		const double               _scaleFactor;  // constant parameter, needs to be copied (-> no & !)
		const std::vector<double>& _values;
		std::vector<double>&       _results;

	};

template <typename T>
	class Increase: public Function
	{

	public:

		Increase(const double& increment,
		         const T&      value,
		         T&            result)
			   : _increment(increment),
			     _value    (value),
			     _result   (result)
		{ }

		virtual ~Increase() { }

		bool
		operator() ()
		{
			_result = _increment + _value;
			return true;
		}

	private:

		const double _increment;  // constant parameter, needs to be copied (-> no & !)
		const T&     _value;
		T&           _result;

	};


	template <>
	class Increase<std::vector<double>>: public Function
	{

	public:

		Increase(const double&            increment,
		      const std::vector<double>&  values,
		      std::vector<double>&        results)
			: _increment(increment),
			  _values   (values),
			  _results  (results)
		{ }

		virtual ~Increase() { }

		bool
		operator() ()
		{
			const size_t sizeVec = _values.size();
			_results.resize(sizeVec);
			for (size_t i = 0; i < sizeVec; ++i) {
				_results[i] = _increment + _values[i];
			}
			return true;
		}

	private:

		const double               _increment;  // constant parameter, needs to be copied (-> no & !)
		const std::vector<double>& _values;
		std::vector<double>&       _results;

	};

	// Performs the kinematic fit of two charged particles:
	class GetChargedKinematicFitting : public Function
	{

	public:

		GetChargedKinematicFitting(const TVector3& Particle1Momentum, // from tracks, extrapolated to Z of SV!
                                   const TVector3& Particle2Momentum,
                                   const double&   Particle1Energy, // sqrt(p_track*p_track + m*m)
                                   const double&   Particle2Energy,
                                   const std::vector<double>& Particle1MomentumCovariance,
                                   const std::vector<double>& Particle2MomentumCovariance,
                                   const double    Particle1Mass,
                                   const double    Particle2Mass,
                                   const double    Mass,
                                   //const double    MassLowerLimit, // not used
                                   //const double    MassUpperLimit, // not used
                                   const double    PrecisionGoal,
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
				//_MassLowerLimit             (MassLowerLimit),
				//_MassUpperLimit             (MassUpperLimit),
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
				0, // MassLowerLimit for chargedFit::massIsInWindow(), which is not used
				0, // MassUpperLimit for chargedFit::massIsInWindow(), which is not used
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
			//std::cout << "Mass after Fit:" << (_ResultLorentzVectors[i]).M() << std::endl;
			//std::cout << "Mass difference:" << (_ResultLorentzVectors[i]).M()-_Mass << std::endl;

			//std::cout << "LV[" << i << "] = (" << _ResultLorentzVectors[i].X() << "," << _ResultLorentzVectors[i].Y() << "," << _ResultLorentzVectors[i].Z() << "," << _ResultLorentzVectors[i].E() << std::endl;

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
        //const double    _MassLowerLimit;
        //const double    _MassUpperLimit;
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

	// calculate collinearity angle, i.e. the angle between the geometrical line from the primary to the secondary vertex
	// and the momentum of the V0 particle thought to have decayed in the secondary vertex.
	// gives the value in rad

	class CalculateCollinearityAngle: public Function
	{

	public:

		CalculateCollinearityAngle(const TVector3& primVertexPos, const TVector3& secVertexPos, const TLorentzVector& V0LVector, double& result)
			: _primVertexPos(primVertexPos),
			  _secVertexPos (secVertexPos),
			  _V0LVector    (V0LVector),
			  _result       (result)
		{ }

		virtual ~CalculateCollinearityAngle() { }

		bool
		operator() ()
		{
			_result = -1.;
			const TVector3 lineToSV( _secVertexPos - _primVertexPos );
			_result = _V0LVector.Angle(lineToSV);
			return true;
		}

	private:

		const TVector3&       _primVertexPos;
		const TVector3&       _secVertexPos;
		const TLorentzVector& _V0LVector;
		double&               _result;

	};


	// Compare indices in vector to indices in vector of vectors
	// Check whether all vectors of indices in `vectorOfVectorsOfIndices` are fully overlapping with original index vector `vectorOfIndices`
	// if there are new indices in one or more vectors, return 1
	// else return 0

	class CompareIndices: public Function
	{

	public:

		CompareIndices(const std::vector<int>& vectorOfIndices, const std::vector<std::vector<int>>& vectorOfVectorsOfIndices, int& result)
			: _vectorOfIndices          (vectorOfIndices),
			  _vectorOfVectorsOfIndices (vectorOfVectorsOfIndices),
			  _result                   (result)
		{ }

		virtual ~CompareIndices() { }

		bool
		operator() ()
		{
			_result = 0;
			int NVectorNonOverlaps = 0;
			for (uint i = 0; i < _vectorOfVectorsOfIndices.size(); i++)
			{
				std::vector<int> vector = _vectorOfVectorsOfIndices[i];
				uint overlap = 0;
				for (uint j = 0; j < vector.size(); j++)
				{
					int currentIndex = vector[j];
					for (uint k = 0; k < _vectorOfIndices.size(); k++)
					{
						if (currentIndex == _vectorOfIndices[k]) { overlap++; break; }
					}
				}
				if (overlap != vector.size())
				{
					NVectorNonOverlaps++;
				}
			}
			if (NVectorNonOverlaps != 0)
			{
				_result = 1;
			}
			return true;
		}

	private:

		const std::vector<int>&               _vectorOfIndices;
		const std::vector<std::vector<int>>&  _vectorOfVectorsOfIndices;
		int&                                  _result;

	};


	// Helper class to print stuff

	class PrintEvent: public Function
	{

	public:

		PrintEvent(const int& runNmbr,const int& spillNmbr,const int& evtNmbr)
			: _runNmbr     (runNmbr),
			  _spillNmbr   (spillNmbr),
			  _evtNmbr     (evtNmbr)
		{ }

		virtual ~PrintEvent() { }

		bool
		operator() ()
		{
			std::cout << _runNmbr << ", " << _spillNmbr << ", " << _evtNmbr << std::endl;
			return true;
		}

	private:

		const int&     _runNmbr;
		const int&     _spillNmbr;
		const int&     _evtNmbr;

	};

}
}
}
}

#endif  // ANTOK_USER_JBECKERS_FUNCTIONS_HPP