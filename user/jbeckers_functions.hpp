#ifndef ANTOK_USER_JBECKERS_FUNCTIONS_HPP
#define ANTOK_USER_JBECKERS_FUNCTIONS_HPP

#include "TVector3.h"
#include "TLorentzVector.h"

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