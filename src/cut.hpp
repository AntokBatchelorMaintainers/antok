#ifndef HLIB_CUT_H
#define HLIB_CUT_H

#include<string>

#include<event.h>

namespace hlib {

	class Cut {
	
	  public:

		virtual bool operator() (const hlib::Event& event) const = 0;

		static const std::string shortname;
		static const std::string longname;
		static const std::string abbreviation;
	
	};

	// Cut on the trigger mask.
	class TrigMask : public Cut {

	  public:
		bool operator() (const hlib::Event& event) const { return (!((event.rawData->TrigMask)&0x1)); };
		static const std::string shortname;
		static const std::string longname;
		static const std::string abbreviation;
	};
	const std::string TrigMask::shortname("trigmask_&_0x1");
	const std::string TrigMask::longname("Trigger Mask & 0x1");
	const std::string TrigMask::abbreviation("tm");

	// Cut on the vertex Z position.
	class VrtxZ : public Cut {

	  public:
		bool operator() (const hlib::Event& event) const { return ((event.rawData->Z_primV > -28.4) ||
		                                                           (event.rawData->Z_primV < -68.4)); };
		static const std::string shortname;
		static const std::string longname;
		static const std::string abbreviation;

	};
	const std::string VrtxZ::shortname("-28.4<vtx_z<-68.4");
	const std::string VrtxZ::longname("Vertex Z in ]-28.4,-68.4[");
	const std::string VrtxZ::abbreviation("vz");

	// Cut on the vertex R.
	class VrtxR : public Cut {

	  public:
		VrtxR(double rmax) { _rmax = rmax; _rmax2 = _rmax*_rmax; };
		bool operator() (const hlib::Event& event) const { return (std::pow(event.rawData->X_primV, 2) +  
																   std::pow(event.rawData->Y_primV, 2) > _rmax2); };
		static const std::string shortname;
		static const std::string longname;
		static const std::string abbreviation;

	  private:
		double _rmax;
		double _rmax2;
	};
	const std::string VrtxR::shortname("vtx_R < 1.75");
	const std::string VrtxR::longname("Vertex.R() < 1.75");
	const std::string VrtxR::abbreviation("vr");


}

#endif

