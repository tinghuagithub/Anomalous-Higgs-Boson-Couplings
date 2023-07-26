  // -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include "NLOAnalysis.h"
#include "VBFRegions.h"
#include <map>
#include <sstream>

namespace Rivet {
 
  
    /// @brief Higgs boson analysis
    ///
    /// @author Terrance Figy <terrance.figy@wichita.edu>
  

   struct Variables {

    Variables(const vector<const Jet*>& jets) {
      FourMomentum j1 = jets.at(0)->momentum();
      FourMomentum j2 = jets.at(1)->momentum();
      jet1pt = j1.pT();
      jet2pt = j2.pT();
      assert(jet1pt > jet2pt);

     

      deltay = fabs(j1.rapidity() - j2.rapidity());
      mjj = (j1 + j2).mass();
      ptbar = (j1.pT()+j2.pT())/2.0;

      FourMomentum gapjet(0., 0., 0., 0.);
      ngapjets = _getNumGapJets(jets, gapjet);

      gapjet1pt = gapjet.pT()/GeV;
      gapjet1y = gapjet.rapidity();
      //double ptbal_vec = (j1 + j2 + lep1->mom() + lep2->mom()).pT();
      //double ptbal_sc = j1.pT() + j2.pT() + lep1->pT() + lep2->pT();
      //ptbalance2 = ptbal_vec / ptbal_sc;

      //double ptbal3_vec = (j1 + j2 + gapjet + lep1->mom() + lep2->mom()).pT();
      //double ptbal3_sc = j1.pT() + j2.pT() + gapjet.pT() + lep1->pT() + lep2->pT();
      //ptbalance3 = ptbal3_vec / ptbal3_sc;

    
      pass_jetveto = gapjet.pT() < 25.0*GeV;
      
    }


    double jet1pt;
    double jet2pt;
    double ptbar;
    double gapjet1pt;
    double gapjet1y;

    double deltay;
    double mjj;
    
    //double deltaphijj;
    //double ptbalance2;
    //double ptbalance3;
    int ngapjets;

    //double dilepton_dr;

    bool pass_jetveto;
    //bool pass_ptbaleff;


  private:

    bool _isBetween(const Jet* probe, const Jet* boundary1, const Jet* boundary2) {
      double y_p = probe->rapidity();
      double y_b1 = boundary1->rapidity();
      double y_b2 = boundary2->rapidity();

      double y_min = std::min(y_b1, y_b2);
      double y_max = std::max(y_b1, y_b2);

      if (y_p > y_min && y_p < y_max) return true;
      else return false;
    }

    int _getNumGapJets(const vector<const Jet*>& jets, FourMomentum& thirdJet) {
      if (jets.size() < 2) return 0;
      // The vector of jets is already sorted by pT. So the boundary jets will be the first two.
      const Jet* bj1 = jets.at(0);
      const Jet* bj2 = jets.at(1);

      int n_between = 0;
      // Start loop at the 3rd hardest pT jet
      for (size_t i = 2; i < jets.size(); ++i) {
        const Jet* j = jets.at(i);
        // If this jet is between the boundary jets and is hard enough, increment counter
        if (_isBetween(j, bj1, bj2)) {
          if (n_between == 0) thirdJet = j->momentum();
          ++n_between;
        }
      }
      return n_between;
    }

  };
	
  class MC_H2JETS : public NLOAnalysis {
    
  protected:
    
    double _R,_jrap,_jpT;
    
    int _VBFDef;
    
  private:
    std::map<std::string,NLOHisto1DPtr> histos;
    
    std::vector<Histo1DPtr> _h_log10_d;
    std::vector<Scatter2DPtr> _h_log10_R;

    
    Profile1DPtr _prof_ngapjets_dy;
    Profile1DPtr _prof_ngapjets_mjj;
    Profile1DPtr _prof_ngapjets_pt;
    Profile1DPtr _prof_ngapjets_ht;
  public:
      /// Constructor
      // DEFAULT_RIVET_ANALYSIS_CTOR(MC_H2JETS);
    MC_H2JETS(string name="MC_H2JETS")
    : NLOAnalysis(name),_h_log10_d(4), _h_log10_R(5) 
    {
      _R=0.4;
      _jpT=25.0;
      _jrap=4.5;
      _VBFDef=0;
    }
    std::complex<double> EPSTENSOR(const FourMomentum &p1,
				   const FourMomentum &p2,
				   const FourMomentum &p3,
				   const FourMomentum &p4)
    {
      return -std::complex<double>(0.,1.)
             *(-p1.z()*p2.x()*p3.E()*p4.y()+p1.z()*p2.y()*p3.E()*p4.x()
	       +p1.E()*p2.x()*p3.z()*p4.y()-p1.E()*p2.y()*p3.z()*p4.x());
    }



      /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs;
      IdentifiedFinalState higgses(PID::HIGGS);
      IdentifiedFinalState photons(PID::PHOTON);
      VetoedFinalState rest(fs);
      rest.addVetoOnThisFinalState(higgses);
      rest.addVetoOnThisFinalState(photons);
      //rest.addVetoId(82);
      declare(fs, "FS");
      declare(higgses, "Higgses");
      declare(photons, "Photons");
      declare(rest, "Rest");
      declare(FastJets(rest, FastJets::ANTIKT, _R), "Jets");
      declare(FastJets(rest,FastJets::KT,0.6),"KTJets");
      defineHistos();
     
    }
    
    
      /// Perform the per-event analysis
    void analyze(const Event& e) {

      double _mH(125.0*GeV);
      double _mHdev(1.0*GeV);
      
      ParticleVector higgses = applyProjection<IdentifiedFinalState>(e, "Higgses").particles();
      ParticleVector photons = applyProjection<IdentifiedFinalState>(e, "Photons").particles();
      ParticleVector rest = applyProjection<VetoedFinalState>(e, "Rest").particles();
     

      const double weight = e.weight();

        // require either one stable Higgs or at least two photons
      if (higgses.size()>1) vetoEvent;
      FourMomentum hmom;
      size_t idph1(0),idph2(0);
      std::vector<FourMomentum> phs;
      if (higgses.size()==1) {
        hmom = higgses[0].momentum();
      }
      else if (photons.size()>1) {
          // reconstruct Higgs from photon pair with correct inv. mass
          // only take first one
        bool foundone(false);
        for (size_t i(0);i<photons.size();++i) {
          for (size_t j(i+1);j<photons.size();++j) {
            if (!foundone &&
                fabs((photons[i].momentum()
                      +photons[j].momentum()).mass()-_mH)<_mHdev) {
              idph1=i; idph2=j;
              hmom = photons[i].momentum()+photons[j].momentum();
              phs.push_back(photons[i].momentum());
              phs.push_back(photons[j].momentum());
              foundone=true;
              break;
            }
          }
          if (foundone) break;
        }
      }
      else vetoEvent;
      
        // check that found one Higgs
      if (higgses.size()==0 && phs.size()!=2) vetoEvent;
      
        // add remaining photons to the remaining final state
      for (size_t i(0); i<photons.size(); ++i) {
        if (idph1==idph2 || (i!=idph1 && i!=idph2)) rest.push_back(photons[i]);
      }
      
        // Get the jet candidates
        // PT ordered
        // Rapdity ordered
        // Other observable
        // Different methods of ordering jets
      auto _jetalgo=apply<FastJets>(e, "Jets");
      
      Jets jets = _jetalgo.jetsByPt(Cuts::pT > 25.*GeV && Cuts::absrap < 4.5);

      //Jets jets = apply<FastJets>(e, "Jets").jetsByPt(15.*GeV);
      
      Jets PTJets,RapJets,alljets;
      foreach (const Jet& jetcand, _jetalgo.pseudoJetsByPt(0.*GeV)) {
	     if (fabs(jetcand.momentum().rapidity()) < _jrap) {
	          alljets.push_back(jetcand);
	       }
      }
      foreach (const Jet& jetcand, _jetalgo.pseudoJetsByPt(_jpT)) {
	     if (fabs(jetcand.momentum().rapidity()) < _jrap) {
	          PTJets.push_back(jetcand);
	        }
      }
      foreach (const Jet& jetcand, _jetalgo.pseudojetsByRapidity(_jpT)) {
	    if (fabs(jetcand.momentum().rapidity()) < _jrap) {
	          RapJets.push_back(jetcand);
	        }
      } 

      vector<const Jet*> good_jets;
     // const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 15*GeV && Cuts::absrap < 4.5);
      foreach(const Jet& j, jets) {
          good_jets.push_back(&j);
      }
     
      // If we don't have at least 2 good jets, we won't use this event.
      if (good_jets.size() < 2) vetoEvent;
      
      Variables vars(good_jets);
      
      // cross check
      if (PTJets.size()!=RapJets.size()) abort();


      //if (jets.size()>1){ // need at least two jets
        const FourMomentum& j1(jets[0].momentum());
        const FourMomentum& j2(jets[1].momentum());  
        
        switch (_VBFDef) {
          case 0:
            if (!VBFRegion0( jets )) vetoEvent;
            break;
          case 1:
            if(!VBFRegion0( jets)) vetoEvent;
            if(!VBFRegion1Tight( jets)) vetoEvent;
            if(!VBFRegion2Tight( jets)) vetoEvent; 
            if(!VBFRegion3( jets)) vetoEvent;
            break;
          case 2:
            if(!VBFRegion0( jets)) vetoEvent;
            if(!VBFRegion1Loose( jets)) vetoEvent;
            if(!VBFRegion2Loose( jets)) vetoEvent;
            break;
          default:
            break;
        }
          // check if event passes cuts
        //double mjjval = (j1+j2).mass()/GeV;
        double dyjjval= fabs(j1.rapidity()-j2.rapidity());
        double y_mid = (j1.rapidity()+j2.rapidity())/(2.0);
          // at this stage we should have events with the first and second jets of at least 30GEV
          // resonance veto snippet to follow (is this robust?)
          //
        /*
         // W/Z veto
         bool leafloop=false;
         for ( Jets::const_iterator jet1=PTJets.begin(); jet1!=PTJets.end() && !leafloop; ++jet1){
         // construct masses
         double m_i    = jet1->momentum().mass()/GeV;
         // check masses
         if ( fabs(m_i   -_mZ)<_mZdev || fabs(m_i   -_mW)<_mWdev ) {
         leafloop=true;
         break;
         }
         // next multiplicity
         for ( Jets::const_iterator jet2=jet1; jet2!=PTJets.end() && !leafloop; ++jet2){
         // construct masses
         double m_ij   = ( jet1->momentum()+jet2->momentum() ).mass()/GeV;
         // check masses
         if ( fabs(m_ij  -_mZ)<_mZdev || fabs(m_ij  -_mW)<_mWdev ) {
         leafloop=true;
         break;
         }
         // next multiplicity
         for ( Jets::const_iterator jet3=jet2; jet3!=PTJets.end() && !leafloop; ++jet3){
         // construct masses
         double m_ijk  = ( jet1->momentum()+jet2->momentum()+jet3->momentum() ).mass()/GeV;
         // check masses
         if ( fabs(m_ijk -_mZ)<_mZdev || fabs(m_ijk -_mW)<_mWdev ) {
         leafloop=true;
         break;
         }
         // next multiplicity
         for ( Jets::const_iterator jet4=jet3; jet4!=PTJets.end() && !leafloop; ++jet4){
         // construct masses
         double m_ijkl = ( jet1->momentum()+jet2->momentum()+jet3->momentum()+jet4->momentum() ).mass()/GeV;
         // check masses
         if ( fabs(m_ijkl-_mZ)<_mZdev || fabs(m_ijkl-_mW)<_mWdev ) {
         leafloop=true;
         break;
         }
         }
         }
         }
         }
         if (leafloop) vetoEvent;
         */
       const FastJets& jetpro = apply<FastJets>(e, "KTJets");
       const auto seq = jetpro.clusterSeq();
       if (!seq) vetoEvent; //< the cseq is the whole point in this sort of analysis!!
      size_t m_njet(4);
    // Jet resolutions and integrated jet rates
       double previous_dij = 10.0;
       for (size_t i = 0; i < min(m_njet,(size_t)seq->n_particles()); ++i) {
       const double d_ij2 = seq->exclusive_dmerge_max(i);
       if (d_ij2 <= 0) continue; ///< @todo Is < 0 possible? Feels like no; I should check ;-)
      // Jet resolution i -> j
       const double d_ij = log10(sqrt(d_ij2));

      // Fill differential jet resolution
      _h_log10_d[i]->fill(d_ij, weight);

      // Fill integrated jet resolution
      for (size_t ibin = 0; ibin < _h_log10_R[i]->numPoints(); ++ibin) {
        Point2D& dp = _h_log10_R[i]->point(ibin);
        if (dp.x() > d_ij && dp.x() < previous_dij) {
          dp.setY(dp.y() + weight);
        }
      }
      previous_dij = d_ij;
    }
    // One remaining integrated jet resolution
    for (size_t ibin = 0; ibin<_h_log10_R[m_njet]->numPoints(); ++ibin) {
      Point2D & dp = _h_log10_R[m_njet]->point(ibin);
      if (dp.x() < previous_dij) {
        dp.setY(dp.y() + weight);
      }
    }

        
        histos["XS"]->fill(0.5,e);
//--------------------------------------------------------------------
 // Calculation of phi_2 from arXiv:1001.3822  
	// phi_2 = azimuthal angle between the vector sum of jets 
	// forward and jets backward of the Higgs boson

	      const FourMomentum& jb(RapJets.front().momentum());
      	const FourMomentum& jf(RapJets.back().momentum());
	
	      FourMomentum vsumf(0.,0.,0.,0.);
	      FourMomentum vsumb(0.,0.,0.,0.);
	      FourMomentum p1(1.,0.,0.,1.);
	      FourMomentum p2(1.,0.,0.,-1.);

	      bool f_nonzero(false),b_nonzero(false);
	      foreach (const Jet& jj, RapJets) {
	      if (jj.momentum().rapidity()>hmom.rapidity()) {
	          vsumf += jj.momentum();
	          f_nonzero = true;
	      }
	      else {
	           vsumb += jj.momentum();
	           b_nonzero = true;
	            }
	      }
      	double phi2(-10.);
	// Calculate phi_2
	      if (f_nonzero && b_nonzero) {
	        phi2 = acos((vsumb.x()*vsumf.x()+vsumb.y()*vsumf.y())/
		      (sqrt(vsumb.x()*vsumb.x()+vsumb.y()*vsumb.y())*
		       sqrt(vsumf.x()*vsumf.x()+vsumf.y()*vsumf.y()))); 
	        if (imag(EPSTENSOR(p1,vsumb,p2,vsumf))<0.) phi2 *= -1.;
	        }
	      else if (!f_nonzero) {
	        vsumb -= jf;
	        phi2 = acos((vsumb.x()*jf.x()+vsumb.y()*jf.y())/
		      (sqrt(vsumb.x()*vsumb.x()+vsumb.y()*vsumb.y())*
		       sqrt(jf.x()*jf.x()+jf.y()*jf.y())));
	        if (imag(EPSTENSOR(p1,vsumb,p2,jf))<0.) phi2 *= -1.;
	        }
	      else { 
	        vsumf -= jb;
	        phi2 = acos((jb.x()*vsumf.x()+jb.y()*vsumf.y())/
		      (sqrt(jb.x()*jb.x()+jb.y()*jb.y())*
		       sqrt(vsumf.x()*vsumf.x()+vsumf.y()*vsumf.y())));  
	        if (imag(EPSTENSOR(p1,jb,p2,vsumf))<0.)  phi2 *= -1.;
	       }
	
	      histos["deltaphi2"]->fill(phi2,e);

          // inclusive histograms
        histos["H_pT_incl"]              ->fill(hmom.pT()/GeV,e);
        histos["jet1_pT_incl"]           ->fill(j1.pT()/GeV,e);
        histos["jet2_pT_incl"]           ->fill(j2.pT()/GeV,e);
        
        histos["jet1_y_incl"]            ->fill(j1.rapidity(),e);
        histos["jet2_y_incl"]            ->fill(j2.rapidity(),e);
        histos["jet1_y_abs_incl"]        ->fill(fabs(j1.rapidity()),e);
        histos["jet2_y_abs_incl"]        ->fill(fabs(j2.rapidity()),e);
        
          // use y* definition here instead of rap() only
        histos["H_y"]                    ->fill((hmom.rapidity()-(j1.rapidity()+j2.rapidity())/(2.0))/dyjjval,e);
        // JB 1.5.2019: Either use y* or y, but this y and y_abs is misleading ???
        histos["H_y_abs"]                ->fill(fabs(hmom.rapidity()),e);
        histos["NJet_excl"]              ->fill(jets.size(),e);
        
        for (size_t i(0); i<15;++i){
          if (jets.size()>i) {
            histos["NJet_incl"]          ->fill(i,e);
            histos["H_j_pT_incl"]        ->fill(jets[i].pT()/GeV,e);
          }
        }
         

	// fill gap jet
        //
         histos["NGapJet_excl"]            ->fill(vars.ngapjets,e);

         for (int i(0); i<15;++i){
           if (vars.ngapjets>i){
            histos["NGapJet_incl"]       ->fill(i,e);
            }
          }
         if(vars.ngapjets>0){
       // cout << "here we go";
        _prof_ngapjets_dy->fill(vars.deltay, vars.ngapjets,weight);
        _prof_ngapjets_mjj->fill(vars.mjj, vars.ngapjets, weight);
	      _prof_ngapjets_pt->fill(vars.ptbar,vars.ngapjets,weight);
        }
        // \Delta R(H,j1), 
        histos["Hj_pT_incl"]             ->fill((hmom+j1).pT()/GeV,e);
        histos["jet1_mass"]              ->fill(j1.mass()/GeV,e);
        histos["jet2_mass"]              ->fill(j2.mass()/GeV,e);
        histos["deltaphi_jj_incl"]       ->fill(deltaPhi(j1,j2),e);
        histos["deltaR_H_jj_incl"]       ->fill(deltaR(hmom,j1+j2),e);
        histos["deltaR_H_j1_incl"]       ->fill(deltaR(hmom,j1),e);
        histos["deltaR_H_j2_incl"]       ->fill(deltaR(hmom,j2),e);
        histos["deltaR_jj_incl"]         ->fill(deltaR(j1,j2),e);
        histos["deltaphi_Hjj_incl"]      ->fill(deltaPhi(hmom,j1+j2),e);
        histos["deltaphi_Hj1_incl"]      ->fill(deltaPhi(hmom,j1),e);
        histos["deltaphi_Hj2_incl"]      ->fill(deltaPhi(hmom,j2),e);
        histos["Hjj_pT_incl"]            ->fill((hmom+j1+j2).pT()/GeV,e);
        // add the sign version for the deltaphi_jj
        double deltaPhiJfJb(0.0);
        if (j1.rapidity()>j2.rapidity()){
          deltaPhiJfJb = j1.phi()-j2.phi();
        }
        else {
          deltaPhiJfJb = j2.phi()-j1.phi();
        }

        if (deltaPhiJfJb > PI){
          deltaPhiJfJb = deltaPhiJfJb - 2.0*PI;
        }
        else if (deltaPhiJfJb < -PI){
          deltaPhiJfJb = deltaPhiJfJb + 2.0*PI;
        }
        histos["deltaphi_jfjb"]       ->fill(deltaPhiJfJb,e);
        // add the y_sep defined in the arXiv:1001.3822
        double y_sep(0.0);
        if (j1.rapidity() > hmom.rapidity() && j2.rapidity() < hmom.rapidity()){
              y_sep = min(fabs(hmom.rapidity()-j1.rapidity()),fabs(hmom.rapidity()-j2.rapidity()));
        }
        else if (j2.rapidity() > hmom.rapidity() && j1.rapidity() < hmom.rapidity()){
              y_sep = min(fabs(hmom.rapidity()-j1.rapidity()),fabs(hmom.rapidity()-j2.rapidity()));
        }
        histos["H_jj_sep"]              ->fill(y_sep,e);

        // JB 1.5.2019: Is this j1j2_pT, name is misleading (compared to deltay_H_jj)?
        histos["H_jj_pT_incl"]           ->fill((j1+j2).pT()/GeV,e);
        histos["deltay_jj"]              ->fill(dyjjval,e);
        histos["Hjj_y"]                  ->fill((hmom+j1+j2).rapidity(),e); 
        histos["deltay_H_j1"]            ->fill(abs(hmom.rapidity()-j1.rapidity()),e);
        histos["deltay_H_j2"]            ->fill(abs(hmom.rapidity()-j2.rapidity()),e);
        histos["dijet_mass"]             ->fill((j1+j2).mass(),e);
        histos["dijet_mass_fine"]        ->fill((j1+j2).mass(),e);
        histos["yj1_yj2"]                ->fill(j1.rapidity()*j2.rapidity(),e);
        // there could be two different definitions
        double pull_angle = fabs(CalculatePullAngle(jets[0], jets[1], 0));
        histos["j1j2_pullangle"]->fill(pull_angle / Rivet::PI,e);
        
        histos["H_dijet_mass"]           ->fill((hmom+j1+j2).mass(),e);
        histos["deltay_H_jj"]            ->fill(fabs((hmom.rapidity()-(j1+j2).rapidity())),e);
        histos["H3_y"]                   ->fill((hmom.rapidity()-(j1.rapidity()+j2.rapidity())/(2.0))/dyjjval,e);
          
          // njets == 2;
        if (jets.size()==2) {
          histos["deltaphi_jj_excl"]     ->fill(deltaPhi(j1,j2),e);
          histos["deltaphi_Hjj_excl"]    ->fill(deltaPhi(hmom,j1+j2),e);
          histos["Hjj_pT_excl"]          ->fill((hmom+j1+j2).pT()/GeV,e);
          // JB 1.5.2019: This seems wrong??? (include the jets!)
          histos["H_jj_pT_excl"]         ->fill(hmom.pT()/GeV,e);
        }
          // njets > 2;
        if (jets.size()>2) {
	  
          const FourMomentum& j3(jets[2].momentum());
          histos["j1j3_mass"]           ->fill((j1+j3).mass(),e);
	        histos["j1j3_mass_fine"]      ->fill((j1+j3).mass(),e);
          histos["j2j3_mass"]           ->fill((j2+j3).mass(),e);
          histos["j2j3_mass_fine"]      ->fill((j2+j3).mass(),e);
          histos["trijet_mass"]         ->fill((j1+j2+j3).mass(),e);
          histos["trijet_mass_fine"]    ->fill((j1+j2+j3).mass(),e);
          histos["H_jjj_pT_incl"]       ->fill((j1+j2+j3).pT()/GeV,e);
          histos["jet3_pT_incl"]        ->fill(j3.pT()/GeV,e);
          histos["deltay_H_j3"]         ->fill(abs(hmom.rapidity()-j3.rapidity()),e);
            // use y* definition here instead of rap() only
            //define z_j3 use the definition in HXWSG paper
          double z_j3=(j3.rapidity()-(j1.rapidity()+j2.rapidity())/(2.0))/dyjjval;
          histos["jet3_y"]              ->fill(z_j3,e);
          histos["jet3_mass"]           ->fill(j3.mass()/GeV,e);
           //add xj3 use the definitoin in HXWSG paper
          double x_j3=min(fabs(j1.rapidity()-j3.rapidity()),fabs(j2.rapidity()-j3.rapidity()));
           //if z_j3 is large, x_j3 is negative
          if (fabs(z_j3)>0.5) x_j3 *=-1;
          histos["jet3_x"]		        ->fill(x_j3,e);
          // add y_star variable
          double y_j3_star = j3.rapidity()-y_mid;
          histos["jet3_y_star"]           ->fill(y_j3_star,e);

	        double pull_angle13 = fabs(CalculatePullAngle(jets[0], jets[2], 0));
          histos["j1j3_pullangle"]->fill(pull_angle13 / Rivet::PI,e);

	        double pull_angle23 = fabs(CalculatePullAngle(jets[1], jets[2], 0));
          histos["j2j3_pullangle"]->fill(pull_angle23 / Rivet::PI,e);

           // Only fill gapjet if it exists.
          if (vars.ngapjets>0){
             histos["gapjet1_pT_incl"]       ->fill(vars.gapjet1pt,e);
	           histos["gapjet1_y_incl"]        ->fill(vars.gapjet1y,e);

             // HT of gap jets start at position 2
            double HT_gapjets(0.0);
            for(size_t i(2);i<good_jets.size();i++) HT_gapjets += good_jets[i]->momentum().pT()/GeV;
            histos["HT_gapjets"]             ->fill(HT_gapjets,e);

            _prof_ngapjets_ht->fill(HT_gapjets,vars.ngapjets,weight);
            
          }
          if (jets.size()==3) {
            histos["H_jjj_pT_excl"]     ->fill((j1+j2+j3).pT()/GeV,e);
            histos["jet3_pT_excl"]      ->fill(j3.pT()/GeV,e);
            histos["jet3_y_incl"]       ->fill(j3.rapidity(),e);
          }
          
        }
        if (jets.size()>3) {
          const FourMomentum& j3(jets[2].momentum());
          const FourMomentum& j4(jets[3].momentum());
          
          histos["j1j4_mass"]          ->fill((j1+j4).mass(),e);
          histos["j1j4_mass_fine"]     ->fill((j1+j4).mass(),e);
          histos["j2j4_mass"]          ->fill((j2+j4).mass(),e);
          histos["j2j4_mass_fine"]     ->fill((j2+j4).mass(),e);
          histos["j3j4_mass"]          ->fill((j3+j4).mass(),e);
          histos["j3j4_mass_fine"]     ->fill((j3+j4).mass(),e);
          histos["j1j2j4_mass"]        ->fill((j1+j2+j4).mass(),e);
          histos["j1j2j4_mass_fine"]   ->fill((j1+j2+j4).mass(),e);
          histos["j1j3j4_mass"]        ->fill((j1+j3+j4).mass(),e);
          histos["j1j3j4_mass_fine"]   ->fill((j1+j3+j4).mass(),e);
          histos["j2j3j4_mass"]        ->fill((j2+j3+j4).mass(),e);
          histos["j2j3j4_mass_fine"]   ->fill((j2+j3+j4).mass(),e);
          histos["j1j2j3j4_mass"]      ->fill((j1+j2+j3+j4).mass(),e);
          histos["j1j2j3j4_mass_fine"] ->fill((j1+j2+j3+j4).mass(),e);
          histos["H_jjj_pT_incl"]      ->fill((j1+j2+j3).pT()/GeV,e);
          histos["jet4_pT_incl"]       ->fill(j4.pT()/GeV,e);
            // use y* definition here instead of rap() only
          double z_j4= (j4.rapidity()-(j1.rapidity()+j2.rapidity())/(2.0))/dyjjval;
          histos["jet4_y"]             ->fill(z_j4,e);
          double y_j4_star = j4.rapidity()-y_mid;
          histos["jet4_y_star"]           ->fill(y_j4_star,e);
          // min function is needed
	  // xj4=min(fabs(j1.rap-j4.rap),fabs(j2.rap-j4.rap))
	  double x_j4=min(fabs(j1.rapidity()-j4.rapidity()),fabs(j2.rapidity()-j4.rapidity()));
	  if (fabs(z_j4)>0.5) x_j4 *=-1;
	  histos["jet4_x"]             ->fill(x_j4,e);
	  histos["jet4_mass"]          ->fill(j4.mass()/GeV,e);
          histos["jet4_y_incl"]        ->fill(j4.rapidity(),e);
          
        }
        double HT_jets(0.),HT_all(0.);
        for(size_t i(0);i<jets.size();i++) HT_jets += jets[i].momentum().pT();
       
        HT_all=HT_jets+hmom.Et();
        histos["HT_jets"]->fill(HT_jets,e);
        histos["HT_all"]->fill(HT_all,e);
        // Scalar transverse momentum sum
        double HT_jets_central(0.), HT_jets_mid(0.);
        for(size_t i(0);i<jets.size();i++){
          if (jets[i].rapidity()>=-0.5 && jets[i].rapidity()<=0.5){
            HT_jets_central += abs(jets[i].momentum().pT());
          }
          if ((jets[i].rapidity()-y_mid)>=-0.5 && (jets[i].rapidity()-y_mid)<=0.5){
             HT_jets_mid += abs(jets[i].momentum().pT());
          }
        }        
        histos["HT_jets_central"]->fill(HT_jets_central,e);
        histos["HT_jets_mid"]->fill(HT_jets_mid,e);
        
      //}
      
    }


    Vector3 CalculatePull(Jet& jet, bool &isCharged) {
      Vector3 pull(0.0, 0.0, 0.0);
      double PT = jet.pT();
      Particles& constituents = jet.particles();
      Particles charged_constituents;
      if (isCharged) {
        for (Particle p : constituents) {
          if (p.charge3() != 0)  charged_constituents += p;
        }
        constituents = charged_constituents;
      }
      // calculate axis
      FourMomentum axis;
      for (Particle p : constituents)  axis += p.momentum();
      Vector3 J(axis.rap(), axis.phi(MINUSPI_PLUSPI), 0.0);
      // calculate pull
      for (Particle p : constituents) {
        Vector3 ri = Vector3(p.rap(), p.phi(MINUSPI_PLUSPI), 0.0) - J;
        while (ri.y() >  Rivet::PI) ri.setY(ri.y() - Rivet::TWOPI);
        while (ri.y() < -Rivet::PI) ri.setY(ri.y() + Rivet::TWOPI);
        pull.setX(pull.x() + (ri.mod() * ri.x() * p.pT()) / PT);
        pull.setY(pull.y() + (ri.mod() * ri.y() * p.pT()) / PT);
      }
      return pull;
    }

    double CalculatePullAngle(Jet& jet1, Jet& axisjet, bool isCharged) {
      Vector3 pull_vector = CalculatePull(jet1, isCharged);
      pull_vector = Vector3(1000.*pull_vector.x(), 1000.*pull_vector.y(), 0.);
      double drap = axisjet.rap() - jet1.rap();
      double dphi = axisjet.phi(MINUSPI_PLUSPI) - jet1.phi(MINUSPI_PLUSPI);
      Vector3 j2_vector(drap, dphi, 0.0);
      return mapAngleMPiToPi(deltaPhi(pull_vector, j2_vector));
    }

      /// Normalise histograms etc., after the run
    void finalize() {

    
      for (std::map<std::string,NLOHisto1DPtr>::iterator
           hit=histos.begin(); hit!=histos.end();hit++)
        hit->second->finalize();
      
      double scalefactor(crossSection()/sumOfWeights());
      
      for (std::map<std::string,NLOHisto1DPtr>::iterator
           hit=histos.begin(); hit!=histos.end();hit++)
        scale(hit->second,scalefactor);

      const double xsec_unitw = crossSection()/picobarn/sumOfWeights();

      size_t m_njet(4);
      for (size_t i = 0; i < m_njet; ++i) {
       scale(_h_log10_d[i], xsec_unitw);
        for (size_t ibin = 0; ibin<_h_log10_R[i]->numPoints(); ++ibin) {
        Point2D& dp = _h_log10_R[i]->point(ibin);
        dp.setY(dp.y()*xsec_unitw);
      }
    }

      for (size_t ibin = 0; ibin < _h_log10_R[m_njet]->numPoints(); ++ibin) {
      Point2D& dp =_h_log10_R[m_njet]->point(ibin);
      dp.setY(dp.y()*xsec_unitw);
      }


      }
    
    
    void defineHistos(){
      histos["XS"] = bookNLOHisto1D("XS",1,0.,1.);
      histos["yj1_yj2"] = bookNLOHisto1D("y1j_yj2",50,-10.0,10.0);
        // fine binnings always end on ""
        // incl. and excl. jet multis
      histos["NJet_excl"] = bookNLOHisto1D("NJet_excl",20,-0.5,19.5);
      histos["NJet_incl"] = bookNLOHisto1D("NJet_incl",20,-0.5,19.5);
        // incl. and excl. gap jet multis
      histos["NGapJet_incl"] = bookNLOHisto1D("NGapJet_incl",20,-0.5,19.5);
      histos["NGapJet_excl"] = bookNLOHisto1D("NGapJet_excl",20,-0.5,19.5);
        // pT(H) in incl. and excl. jet bins
      histos["H_pT_incl"] = bookNLOHisto1D("H_pT_incl",50,0,500);
      histos["H_pT_excl"]= bookNLOHisto1D("H_pT_excl",50,0,500);
      
      histos["H_j_pT_incl"]= bookNLOHisto1D("H_j_pT_incl",50,0,500);
      
      histos["H_jj_pT_incl"] = bookNLOHisto1D("H_jj_pT_incl",50,0,500);
      histos["H_jj_pT_excl"] = bookNLOHisto1D("H_jj_pT_excl",50,0,500);
      
      histos["H_jjj_pT_incl"] = bookNLOHisto1D("H_jjj_pT_incl",50,0,500);
      histos["H_jjj_pT_excl"] = bookNLOHisto1D("H_jjj_pT_excl",50,0,500);
      
        // pT(H+nj) in incl. and excl. jet bins
      histos["Hj_pT_incl"]= bookNLOHisto1D("Hj_pT_incl",50,0,500);
      
      histos["Hjj_pT_incl"] = bookNLOHisto1D("Hjj_pT_incl",50,0,500);
      histos["Hjj_pT_excl"] = bookNLOHisto1D("Hjj_pT_excl",50,0,500);
      
        // pT(j) in incl. and excl. jet bins
      histos["jet1_pT_incl"] = bookNLOHisto1D("jet1_pT_incl",40,0,400);
      histos["jet2_pT_incl"] = bookNLOHisto1D("jet2_pT_incl",30,0,300);
      histos["jet3_pT_incl"] = bookNLOHisto1D("jet3_pT_incl",50,0,500);
      histos["jet4_pT_incl"] = bookNLOHisto1D("jet4_pT_incl",50,0,500);
      // pT(hardest gap jet)

      histos["gapjet1_pT_incl"] = bookNLOHisto1D("gapjet1_pt_incl",50,0,500);

        // jet rapidities
      histos["jet1_y_incl"] = bookNLOHisto1D("jet1_y_incl",50,-5.0,5.0);
      histos["jet2_y_incl"] = bookNLOHisto1D("jet2_y_incl",50,-5.0,5.0);
      histos["jet3_y_incl"] = bookNLOHisto1D("jet3_y_incl",50,-5.0,5.0);
      histos["jet4_y_incl"] = bookNLOHisto1D("jet4_y_incl",50,-5.0,5.0);
      histos["jet3_y_star"] = bookNLOHisto1D("jet3_y_star",50,-5.0,5.0);
      histos["jet4_y_star"] = bookNLOHisto1D("jet4_y_star",50,-5.0,5.0);
      histos["jet1_y_abs_incl"] = bookNLOHisto1D("jet1_y_abs_incl",50,0,8.0);
      histos["jet2_y_abs_incl"] = bookNLOHisto1D("jet2_y_abs_incl",50,0,8.0);
      
        //
      histos["jet3_pT_excl"] = bookNLOHisto1D("jet3_pT_excl",40,0,400);
      // rapidity of gap jet
      histos["gapjet1_y_incl"]=bookNLOHisto1D("gapjet1_y_incl",50,-5.0,5.0);

        // inclusive Higgs and jet rapidities
      histos["H_y"] = bookNLOHisto1D("H_y",25,-5,5);
      histos["H_y_abs"] = bookNLOHisto1D("H_y_abs",50,0.0,8.0);
      histos["Hjj_y"] = bookNLOHisto1D("Hjj_y",80,-8.0,8.0);
      histos["jet3_y"] = bookNLOHisto1D("jet3_y",25,-5,5);
      histos["jet4_y"] = bookNLOHisto1D("jet4_y",25,-5,5);
      histos["jet3_x"] = bookNLOHisto1D("jet3_x",25,-5,5);
      histos["jet4_x"] = bookNLOHisto1D("jet4_x",25,-5,5);  

      // y_sep as in arXiv:1001.3822
      histos["H_jj_sep"] = bookNLOHisto1D("H_jj_sep",50,0,10);


        // \Delta y(jj), \Delta y(H,j1),\Delta y(H,j2),\Delta y(H,j3)
      histos["deltay_jj"] = bookNLOHisto1D("deltay_jj",50,0,10);
      histos["deltay_H_j1"] = bookNLOHisto1D("deltay_H_j1",50,0,10);
      histos["deltay_H_j2"] = bookNLOHisto1D("deltay_H_j2",50,0,10);
      histos["deltay_H_j3"] = bookNLOHisto1D("deltay_H_j3",50,0,10);
        // Jet pull angle (J1,J2) 
      histos["j1j2_pullangle"] = bookNLOHisto1D("j1j2_pullangle",20,0,1);
       // Jet pull angle (J1,J3)
      histos["j1j3_pullangle"] = bookNLOHisto1D("j1j3_pullangle",20,0,1);   
       // Jet pull angle (J2,J3)
      histos["j2j3_pullangle"] = bookNLOHisto1D("j2j3_pullangle",20,0,1);         

      histos["H3_y"] = bookNLOHisto1D("H3_y",25,-3,3);
        //m(j)
      histos["jet1_mass"] = bookNLOHisto1D("jet1_mass",100,0,100);
      histos["jet2_mass"] = bookNLOHisto1D("jet2_mass",100,0,100);
      histos["jet3_mass"] = bookNLOHisto1D("jet3_mass",100,0,100);
      histos["jet4_mass"] = bookNLOHisto1D("jet4_mass",100,0,100);
      
        // m(jj)
        // 2 jet
      histos["dijet_mass"] = bookNLOHisto1D("dijet_mass",150,0,6000);
     
        // 3 jet
      histos["j1j3_mass"] = bookNLOHisto1D("j1j3_mass",150,0,6000);
      histos["j2j3_mass"] = bookNLOHisto1D("j2j3_mass",150,0,6000);
        // 4 jet
      histos["j1j4_mass"] = bookNLOHisto1D("j1j4_mass",150,0,6000);
      histos["j2j4_mass"] = bookNLOHisto1D("j2j4_mass",150,0,6000);
      histos["j3j4_mass"] = bookNLOHisto1D("j3j4_mass",150,0,6000);
      
      
        // 2 jet
      histos["dijet_mass_fine"] = bookNLOHisto1D("dijet_mass_fine",50,0,200);
        // 3 jet
      histos["j1j3_mass_fine"] = bookNLOHisto1D("j1j3_mass_fine",50,0,200);
      histos["j2j3_mass_fine"] = bookNLOHisto1D("j2j3_mass_fine",50,0,200);
        // 4 jet
      histos["j1j4_mass_fine"] = bookNLOHisto1D("j1j4_mass_fine",50,0,200);
      histos["j2j4_mass_fine"] = bookNLOHisto1D("j2j4_mass_fine",50,0,200);
      histos["j3j4_mass_fine"] = bookNLOHisto1D("j3j4_mass_fine",50,0,200);
      
      
        // m(jjj)
        // 3 jet
      histos["trijet_mass"] = bookNLOHisto1D("trijet_mass",50,0,1000);
      histos["trijet_mass_fine"] = bookNLOHisto1D("trijet_mass_fine",50,0,200);
        // 4 jet
      histos["j1j2j4_mass"] = bookNLOHisto1D("j1j2j4_mass",50,0,1000);
      histos["j1j2j4_mass_fine"] = bookNLOHisto1D("j1j2j4_mass_fine",50,0,200);
      
      histos["j1j3j4_mass"] = bookNLOHisto1D("j1j3j4_mass",50,0,1000);
      histos["j1j3j4_mass_fine"] = bookNLOHisto1D("j1j3j4_mass_fine",50,0,200);
      
      histos["j2j3j4_mass"] = bookNLOHisto1D("j2j3j4_mass",50,0,500);
      histos["j2j3j4_mass_fine"] = bookNLOHisto1D("j2j3j4_mass_fine",50,0,200);
      
        //m(jjjj)
      histos["j1j2j3j4_mass"] = bookNLOHisto1D("j1j2j3j4_mass",50,0,500);
      histos["j1j2j3j4_mass_fine"] = bookNLOHisto1D("j1j2j3j4_mass_fine",50,0,200);
      
        // m(Hjj)
      histos["H_dijet_mass"]= bookNLOHisto1D("H_dijet_mass",50,0,1000);
      
        // \Delta\phi(H,jj) incl. and excl. \Delta\phi(H,j1), \Delta\phi(H,j2)
      histos["deltaphi_jj_incl"] = bookNLOHisto1D("deltaphi_jj_incl",30,0.0,PI);
      histos["deltaphi_jj_excl"] = bookNLOHisto1D("deltaphi_jj_excl",30,0.0,PI);
      histos["deltaphi_jfjb"] = bookNLOHisto1D("deltaphi_jfjb",60,-PI,PI);
      histos["deltaphi_Hjj_incl"] = bookNLOHisto1D("deltaphi_Hjj_incl",30,0.0,PI);
      histos["deltaphi_Hjj_excl"] = bookNLOHisto1D("deltaphi_Hjj_excl",30,0.0,PI);
      histos["deltaphi_Hj1_incl"] = bookNLOHisto1D("deltaphi_Hj1_incl",30,0.0,PI);
      histos["deltaphi_Hj2_incl"] = bookNLOHisto1D("deltaphi_Hj2_incl",30,0.0,PI);
      // incl.
      histos["deltaphi2"] = bookNLOHisto1D("deltaphi2",60,-PI,PI);
      // \Delta R
      histos["deltaR_H_jj_incl"] = bookNLOHisto1D("deltaR_H_jj_incl",40,0.0,8);
      histos["deltaR_H_j1_incl"] = bookNLOHisto1D("deltaR_H_j1_incl",40,0.0,8);
      histos["deltaR_H_j2_incl"] = bookNLOHisto1D("deltaR_H_j2_incl",40,0.0,8);
      histos["deltaR_jj_incl"] = bookNLOHisto1D("deltaR_jj_incl",40,0.0,8);
        // \Delta y(H,jj)
      histos["deltay_H_jj"] = bookNLOHisto1D("deltay_H_jj",25,0,8);
      
        // HT
      histos["HT_all"] = bookNLOHisto1D("HT_all",150,0,6000);
      histos["HT_jets"] = bookNLOHisto1D("HT_jets",150,0,6000);
      histos["HT_gapjets"] = bookNLOHisto1D("HT_gapjets",100,0,1000);
      histos["HT_jets_central"] = bookNLOHisto1D("HT_jets_central",150,0,6000);
      histos["HT_jets_mid"] = bookNLOHisto1D("HT_jets_mid",150,0,6000);
      // add jet splittings (Is 4 jets enought?)

      const double sqrts = sqrtS() ? sqrtS() : 14000.*GeV;

      size_t m_njet(4);
      for (size_t i = 0; i < m_njet; ++i) {
      string dname = "log10_d_" + to_str(i) + to_str(i+1);
      _h_log10_d[i] = bookHisto1D(dname, 100, 0.2, log10(0.5*sqrts/GeV));
      string Rname = "log10_R_" + to_str(i);
      _h_log10_R[i] = bookScatter2D(Rname, 50, 0.2, log10(0.5*sqrts/GeV));
      }
      string Rname = "log10_R_" + to_str(m_njet);
      _h_log10_R[m_njet] = bookScatter2D(Rname, 50, 0.2, log10(0.5*sqrts/GeV));




      // define profiles
      _prof_ngapjets_mjj = bookProfile1D("avgNGapJets_mjj",20,0.0,3000.0,"avgNGapJets","mjj");
      _prof_ngapjets_dy = bookProfile1D("avgNGapJets_dy",20,0.0,8.0,"avgNGapJets","dyjj");
      _prof_ngapjets_pt = bookProfile1D("avgNGapJets_pt",50,0.0,500.0,"avgNGapJets","ptbar");
      _prof_ngapjets_ht = bookProfile1D("avgNGapJets_ht",50,0.0,500.0,"avgNGapJets","ht");
    }
    
    
    
    
  };
  
  
    // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_H2JETS);
  
  struct MC_H2JETS_01_0 : public MC_H2JETS {    MC_H2JETS_01_0() : MC_H2JETS("MC_H2JETS_01_INC") { _R = 0.1; _VBFDef=0;}  };
  struct MC_H2JETS_02_0 : public MC_H2JETS {    MC_H2JETS_02_0() : MC_H2JETS("MC_H2JETS_02_INC") { _R = 0.2; _VBFDef=0;}  };
  struct MC_H2JETS_03_0 : public MC_H2JETS {    MC_H2JETS_03_0() : MC_H2JETS("MC_H2JETS_03_INC") { _R = 0.3; _VBFDef=0;}  };
  struct MC_H2JETS_04_0 : public MC_H2JETS {    MC_H2JETS_04_0() : MC_H2JETS("MC_H2JETS_04_INC") { _R = 0.4; _VBFDef=0;}  };
  struct MC_H2JETS_05_0 : public MC_H2JETS {    MC_H2JETS_05_0() : MC_H2JETS("MC_H2JETS_05_INC") { _R = 0.5; _VBFDef=0;}  };
  struct MC_H2JETS_06_0 : public MC_H2JETS {    MC_H2JETS_06_0() : MC_H2JETS("MC_H2JETS_06_INC") { _R = 0.6; _VBFDef=0;}  };
  struct MC_H2JETS_07_0 : public MC_H2JETS {    MC_H2JETS_07_0() : MC_H2JETS("MC_H2JETS_07_INC") { _R = 0.7; _VBFDef=0;}  };
  struct MC_H2JETS_08_0 : public MC_H2JETS {    MC_H2JETS_08_0() : MC_H2JETS("MC_H2JETS_08_INC") { _R = 0.8; _VBFDef=0;}  };
  struct MC_H2JETS_09_0 : public MC_H2JETS {    MC_H2JETS_09_0() : MC_H2JETS("MC_H2JETS_09_INC") { _R = 0.9; _VBFDef=0;}  };
  struct MC_H2JETS_10_0 : public MC_H2JETS {    MC_H2JETS_10_0() : MC_H2JETS("MC_H2JETS_10_INC") { _R = 1.0; _VBFDef=0;}  };
  struct MC_H2JETS_11_0 : public MC_H2JETS {    MC_H2JETS_11_0() : MC_H2JETS("MC_H2JETS_11_INC") { _R = 1.1; _VBFDef=0; }  };
  struct MC_H2JETS_12_0 : public MC_H2JETS {    MC_H2JETS_12_0() : MC_H2JETS("MC_H2JETS_12_INC") { _R = 1.2; _VBFDef=0; }  };
  struct MC_H2JETS_13_0 : public MC_H2JETS {    MC_H2JETS_13_0() : MC_H2JETS("MC_H2JETS_13_INC") { _R = 1.3; _VBFDef=0; }  };
  struct MC_H2JETS_14_0 : public MC_H2JETS {    MC_H2JETS_14_0() : MC_H2JETS("MC_H2JETS_14_INC") { _R = 1.4; _VBFDef=0; }  };
  struct MC_H2JETS_15_0 : public MC_H2JETS {    MC_H2JETS_15_0() : MC_H2JETS("MC_H2JETS_15_INC") { _R = 1.5; _VBFDef=0; }  };


  struct MC_H2JETS_01_1 : public MC_H2JETS {    MC_H2JETS_01_1() : MC_H2JETS("MC_H2JETS_01_TIGHT") { _R = 0.1; _VBFDef=1;}  };
  struct MC_H2JETS_02_1 : public MC_H2JETS {    MC_H2JETS_02_1() : MC_H2JETS("MC_H2JETS_02_TIGHT") { _R = 0.2; _VBFDef=1;}  };
  struct MC_H2JETS_03_1 : public MC_H2JETS {    MC_H2JETS_03_1() : MC_H2JETS("MC_H2JETS_03_TIGHT") { _R = 0.3; _VBFDef=1;}  };
  struct MC_H2JETS_04_1 : public MC_H2JETS {    MC_H2JETS_04_1() : MC_H2JETS("MC_H2JETS_04_TIGHT") { _R = 0.4; _VBFDef=1;}  };
  struct MC_H2JETS_05_1 : public MC_H2JETS {    MC_H2JETS_05_1() : MC_H2JETS("MC_H2JETS_05_TIGHT") { _R = 0.5; _VBFDef=1;}  };
  struct MC_H2JETS_06_1 : public MC_H2JETS {    MC_H2JETS_06_1() : MC_H2JETS("MC_H2JETS_06_TIGHT") { _R = 0.6; _VBFDef=1;}  };
  struct MC_H2JETS_07_1 : public MC_H2JETS {    MC_H2JETS_07_1() : MC_H2JETS("MC_H2JETS_07_TIGHT") { _R = 0.7; _VBFDef=1;}  };
  struct MC_H2JETS_08_1 : public MC_H2JETS {    MC_H2JETS_08_1() : MC_H2JETS("MC_H2JETS_08_TIGHT") { _R = 0.8; _VBFDef=1;}  };
  struct MC_H2JETS_09_1 : public MC_H2JETS {    MC_H2JETS_09_1() : MC_H2JETS("MC_H2JETS_09_TIGHT") { _R = 0.9; _VBFDef=1;}  };
  struct MC_H2JETS_10_1 : public MC_H2JETS {    MC_H2JETS_10_1() : MC_H2JETS("MC_H2JETS_10_TIGHT") { _R = 1.0; _VBFDef=1;}  };
  struct MC_H2JETS_11_1 : public MC_H2JETS {    MC_H2JETS_11_1() : MC_H2JETS("MC_H2JETS_11_TIGHT") { _R = 1.1; _VBFDef=1; }  };
  struct MC_H2JETS_12_1 : public MC_H2JETS {    MC_H2JETS_12_1() : MC_H2JETS("MC_H2JETS_12_TIGHT") { _R = 1.2; _VBFDef=1; }  };
  struct MC_H2JETS_13_1 : public MC_H2JETS {    MC_H2JETS_13_1() : MC_H2JETS("MC_H2JETS_13_TIGHT") { _R = 1.3; _VBFDef=1; }  };
  struct MC_H2JETS_14_1 : public MC_H2JETS {    MC_H2JETS_14_1() : MC_H2JETS("MC_H2JETS_14_TIGHT") { _R = 1.4; _VBFDef=1; }  };
  struct MC_H2JETS_15_1 : public MC_H2JETS {    MC_H2JETS_15_1() : MC_H2JETS("MC_H2JETS_15_TIGHT") { _R = 1.5; _VBFDef=1; }  };


  
  struct MC_H2JETS_01_2 : public MC_H2JETS {    MC_H2JETS_01_2() : MC_H2JETS("MC_H2JETS_01_LOOSE") { _R = 0.1; _VBFDef=2;}  };
  struct MC_H2JETS_02_2 : public MC_H2JETS {    MC_H2JETS_02_2() : MC_H2JETS("MC_H2JETS_02_LOOSE") { _R = 0.2; _VBFDef=2;}  };
  struct MC_H2JETS_03_2 : public MC_H2JETS {    MC_H2JETS_03_2() : MC_H2JETS("MC_H2JETS_03_LOOSE") { _R = 0.3; _VBFDef=2;}  };
  struct MC_H2JETS_04_2 : public MC_H2JETS {    MC_H2JETS_04_2() : MC_H2JETS("MC_H2JETS_04_LOOSE") { _R = 0.4; _VBFDef=2;}  };
  struct MC_H2JETS_05_2 : public MC_H2JETS {    MC_H2JETS_05_2() : MC_H2JETS("MC_H2JETS_05_LOOSE") { _R = 0.5; _VBFDef=2;}  };
  struct MC_H2JETS_06_2 : public MC_H2JETS {    MC_H2JETS_06_2() : MC_H2JETS("MC_H2JETS_06_LOOSE") { _R = 0.6; _VBFDef=2;}  };
  struct MC_H2JETS_07_2 : public MC_H2JETS {    MC_H2JETS_07_2() : MC_H2JETS("MC_H2JETS_07_LOOSE") { _R = 0.7; _VBFDef=2;}  };
  struct MC_H2JETS_08_2 : public MC_H2JETS {    MC_H2JETS_08_2() : MC_H2JETS("MC_H2JETS_08_LOOSE") { _R = 0.8; _VBFDef=2;}  };
  struct MC_H2JETS_09_2 : public MC_H2JETS {    MC_H2JETS_09_2() : MC_H2JETS("MC_H2JETS_09_LOOSE") { _R = 0.9; _VBFDef=2;}  };
  struct MC_H2JETS_10_2 : public MC_H2JETS {    MC_H2JETS_10_2() : MC_H2JETS("MC_H2JETS_10_LOOSE") { _R = 1.0; _VBFDef=2;}  };
  struct MC_H2JETS_11_2 : public MC_H2JETS {    MC_H2JETS_11_2() : MC_H2JETS("MC_H2JETS_11_LOOSE") { _R = 1.1; _VBFDef=2; }  };
  struct MC_H2JETS_12_2 : public MC_H2JETS {    MC_H2JETS_12_2() : MC_H2JETS("MC_H2JETS_12_LOOSE") { _R = 1.2; _VBFDef=2; }  };
  struct MC_H2JETS_13_2 : public MC_H2JETS {    MC_H2JETS_13_2() : MC_H2JETS("MC_H2JETS_13_LOOSE") { _R = 1.3; _VBFDef=2; }  };
  struct MC_H2JETS_14_2 : public MC_H2JETS {    MC_H2JETS_14_2() : MC_H2JETS("MC_H2JETS_14_LOOSE") { _R = 1.4; _VBFDef=2; }  };
  struct MC_H2JETS_15_2 : public MC_H2JETS {    MC_H2JETS_15_2() : MC_H2JETS("MC_H2JETS_15_LOOSE") { _R = 1.5; _VBFDef=2; }  };
  



  
  DECLARE_RIVET_PLUGIN(MC_H2JETS_01_0);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_02_0);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_03_0);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_04_0);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_05_0);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_06_0);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_07_0);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_08_0);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_09_0);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_10_0);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_11_0);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_12_0);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_13_0);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_14_0);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_15_0);

  DECLARE_RIVET_PLUGIN(MC_H2JETS_01_1);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_02_1);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_03_1);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_04_1);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_05_1);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_06_1);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_07_1);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_08_1);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_09_1);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_10_1);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_11_1);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_12_1);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_13_1);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_14_1);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_15_1);

  DECLARE_RIVET_PLUGIN(MC_H2JETS_01_2);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_02_2);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_03_2);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_04_2);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_05_2);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_06_2);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_07_2);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_08_2);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_09_2);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_10_2);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_11_2);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_12_2);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_13_2);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_14_2);
  DECLARE_RIVET_PLUGIN(MC_H2JETS_15_2);

}
