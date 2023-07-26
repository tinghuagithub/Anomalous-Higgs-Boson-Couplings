#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include <map>
#include <sstream>

namespace Rivet {
  

  bool VBFRegion0( Jets jets){
        const FourMomentum& j1(jets[0].momentum());
        const FourMomentum& j2(jets[1].momentum());
          // require first two jets the have at least 25GeV of pT
        return (j1.pT()/GeV > 25.0 ) && (j2.pT()/GeV > 25.0 );
  }
  bool VBFRegion1Tight( Jets jets){
        const FourMomentum& j1(jets[0].momentum());
        const FourMomentum& j2(jets[1].momentum());
	double mjjval = (j1+j2).mass()/GeV;
          // invariant mass cut
        return (mjjval > 600.0);
  }

  bool VBFRegion1Loose( Jets jets){
        const FourMomentum& j1(jets[0].momentum());
        const FourMomentum& j2(jets[1].momentum());
	double mjjval = (j1+j2).mass()/GeV;
          // invariant mass cut
        return (mjjval > 200.0);
  }
  

  bool VBFRegion2Tight( Jets jets ){
        const FourMomentum& j1(jets[0].momentum());
        const FourMomentum& j2(jets[1].momentum());
	double dyjjval= fabs(j1.rapidity()-j2.rapidity());
	// rapidity gap
	return (dyjjval>4.5);
  }

  bool VBFRegion2Loose( Jets jets ){
        const FourMomentum& j1(jets[0].momentum());
        const FourMomentum& j2(jets[1].momentum());
	double dyjjval= fabs(j1.rapidity()-j2.rapidity());
	// rapidity gap
	return (dyjjval>1.0);
  }


  bool VBFRegion3( Jets jets){
	const FourMomentum& j1(jets[0].momentum());
        const FourMomentum& j2(jets[1].momentum());
	double yopp = j1.rapidity()*j2.rapidity();
	// opposite detectors 
	return (yopp < 0.0);
  }

  

}
