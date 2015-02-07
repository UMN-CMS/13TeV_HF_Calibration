// -*- C++ -*-
//
// Package:    base/calib
// Class:      calib
// 
/**\class calib calib.cc base/calib/plugins/calib.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Quynh Minh Nguyen
//         Created:  Fri, 15 Aug 2014 18:13:15 GMT
//
//

//my include file

//from Zack
#include "Genconverter.h" 

//#include "z_plots.h"
#include "tower_plotter.h"

// system include files
#include <memory>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// CMSSW
#include "FWCore/ServiceRegistry/interface/Service.h" // edm::Service
#include "CommonTools/UtilAlgos/interface/TFileService.h" // TFileService

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h" // GsfElectron
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h" // reco::RecoEcalCandidate

#include "DataFormats/Math/interface/deltaPhi.h"

#include <iostream>
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h" // reco::RecoEcalCandidateCollection

#include <TH1D.h>

//which one is not needed ?
#include "DataFormats/Common/interface/Handle.h" // edm::Handle
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h" // reco::PhotonCollection
#include "DataFormats/EgammaReco/interface/HFEMClusterShape.h" // reco::HFEMClusterShape
#include "DataFormats/EgammaReco/interface/HFEMClusterShapeAssociation.h" // reco::HFEMClusterShapeAssociationCollection
#include "DataFormats/EgammaReco/interface/HFEMClusterShapeFwd.h" // reco::HFEMClusterShapeRef,
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h" // reco::SuperClusterCollection, reco::SuperClusterRef
#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h" // EgammaCutBasedEleId::PassWP, EgammaCutBasedEleId::*
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h" 
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" // PileupSummaryInfo
#include "DataFormats/HLTReco/interface/TriggerEvent.h" // trigger::TriggerEvent


#include "DataFormats/EgammaCandidates/interface/GsfElectron.h" // reco::GsfElectron
#include "DataFormats/EgammaCandidates/interface/Photon.h" // reco::Photon
#include "DataFormats/HLTReco/interface/TriggerObject.h" // trigger::TriggerObject
#include "DataFormats/HepMCCandidate/interface/GenParticle.h" // reco::GenParticle
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h" // reco::RecoEcalCandidate
#include "FWCore/Framework/interface/Event.h" // edm::Event, edm::EventSetup
#include "FWCore/ParameterSet/interface/ParameterSet.h" // edm::ParameterSet
#include "FWCore/Utilities/interface/InputTag.h" // edm::InputTag
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h" // edm::LumiReWeighting

//
// class declaration
//

class calib : public edm::EDAnalyzer { //how to put this in a separate file ?
   public:
      explicit calib(const edm::ParameterSet&); //constructor
      ~calib(); //destructor

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      double pred_hf_e_energy(const reco::GsfElectron gsf_electron, const reco::RecoEcalCandidate hf_electron) ;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
        
      //make some histo to count events
      TH1D* all_events_histo_;
      TH1D* all_gsf_e_histo_;
      TH1D* all_hf_e_histo_;
        
      // ----------member data ---------------------------
	  TH1D* mass_histo_; //a private variable
	  //GSF
	  TH1D* pt_gsf_e_histo_;
	  TH1D* eta_gsf_e_histo_;
	  TH1D* phi_gsf_e_histo_;
	  //HF
	  TH1D* pt_hf_e_histo_;
	  TH1D* eta_hf_e_histo_;
	  TH1D* phi_hf_e_histo_;
	  //HF predicted energy
	  TH1D* pred_hf_e_energy_histo_;
	  TH1D* reco_over_pred_hf_e_energy_histo_;
	  
	  //GEN (GENERATOR) PARTICLES
	  TH1D* pt_gen_e0_histo_;
	  TH1D* eta_gen_e0_histo_;
	  TH1D* phi_gen_e0_histo_;
	  
	  TH1D* pt_gen_e1_histo_;
	  TH1D* eta_gen_e1_histo_;
	  TH1D* phi_gen_e1_histo_;
	  
	  TH1D* mass_gen_Z_histo_;
	  TH1D* eta_gen_Z_histo_;
	  TH1D* phi_gen_Z_histo_;
	  
	  TowerPlotter* tower_plotter_obj_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
calib::calib(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
    //HISTOGRAM
    edm::Service<TFileService> fs;
    mass_histo_ = fs->make<TH1D>("mass_histo", "Mass_{ee}", 100, 50, 150);// make hist within the constructor
   //some counter
   all_events_histo_ = fs->make<TH1D>("all events statistics", " 1: Events processed - 3:GSF_reco  - 4:GSF_Gen - 6:HF_reco -7:HF_gen - 9:Zreco -10:Z_Gen", 20, 0.5, 10.5);
   all_gsf_e_histo_ = fs->make<TH1D>("all_gsf_e", "Pt_all_{gsf_e}", 100, 0, 100);
   all_hf_e_histo_ = fs->make<TH1D>("all_hf_e", "Pt_all_{hf_e}", 100, 0, 100);
   
   //hf
	pt_hf_e_histo_ = fs->make<TH1D>("pt_hf_e_histo", "Pt_{hf_e}", 100, 0, 90);
	eta_hf_e_histo_ = fs->make<TH1D>("eta_hf_e_histo", "eta_{hf_e}", 100, 1, 5);
	phi_hf_e_histo_ = fs->make<TH1D>("phi_hf_e_histo", "phi_{hf_e}", 100, -3.2, 3.2);//radian
	//hf predicted energy // or ratio
	pred_hf_e_energy_histo_ = fs->make<TH1D>("Predicted_energy_hf_e_histo", "Predicted_energy_{hf_e}", 100, 200, 500);
	
	reco_over_pred_hf_e_energy_histo_ = fs->make<TH1D>("Reco/Predicted_energy_hf_e_histo", "Reco/Predicted_energy_{hf_e}", 100, 0.8, 1.2);
	
   //gsf
  	pt_gsf_e_histo_ = fs->make<TH1D>("pt_gsf_e_histo", "Pt_{gsf_e}", 100, 0, 100);
  	eta_gsf_e_histo_ = fs->make<TH1D>("eta_gsf_e_histo", "eta_{gsf_e}", 100, 1, 5);
  	phi_gsf_e_histo_ = fs->make<TH1D>("phi_gsf_e_histo", "phi_{gsf_e}", 100, -3.2, 3.2);
  	
  	//all Generators particles histo
 	mass_gen_Z_histo_ =fs->make<TH1D>("mass_gen_Z_histo", "mass Z^{0}", 100, 60, 120);
	eta_gen_Z_histo_ =fs->make<TH1D>("eta_gen_Z_histo", "Eta Z^{0}", 100, 1, 5);
	phi_gen_Z_histo_ = fs-> make<TH1D>("phi_gen_Z_histo", "Phi Z^{0}", 100, -3.2, 3.2); 	   
    
    pt_gen_e1_histo_=fs->make<TH1D>("pt_gen_e1_histo", "Pt e1^{0}", 100, 0, 90);
	eta_gen_e1_histo_=fs->make<TH1D>("eta_gen_e1_histo", "Eta e1^{0}", 100, 1, 5);
	phi_gen_e1_histo_= fs-> make<TH1D>("phi_gen_e1_histo", "Phi e1^{0}", 100, -3.2, 3.2); 
	  
	pt_gen_e0_histo_=fs->make<TH1D>("pt_gen_e0_histo", "Pt e0^{0}", 100, 1, 100);
	eta_gen_e0_histo_=fs->make<TH1D>("eta_gen_e0_histo", "Eta e0^{0}", 100, 1, 5);
	phi_gen_e0_histo_= fs-> make<TH1D>("phi_gen_e0_histo", "Phi e0^{0}", 100, -3.2, 3.2); 

    // Set up zplots  // Creat a folder name "Zplots"
    TFileDirectory tdir_for_zplots(fs->mkdir("ZPlots"));// https://github.com/UMN-CMS/ZFinder/blob/master/ZFinder/Event/src/ZFinder.cc:150   // this main folder is created in main() 

	// pass "Zplot" to tower plotter obj
	//can't use z_mass_histo_(tdir_for_zplots)
//	TowerPlotter* tower_plotter_obj_= new TowerPlotter(tdir_for_zplots);
	tower_plotter_obj_= new TowerPlotter(tdir_for_zplots);
}
 
//destructor
calib::~calib()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   delete tower_plotter_obj_;
}


//
// member functions
//

//calculated predicted energy of HF electrons
double calib::pred_hf_e_energy(const reco::GsfElectron gsf_electron, const reco::RecoEcalCandidate hf_electron)//const ?
{ 
    using namespace edm;//neccessary ?
    const double ELECTRON_MASS = 5.109989e-4;
    math::PtEtaPhiMLorentzVector hf_e_lv(hf_electron.pt(), hf_electron.eta(), hf_electron.phi(), ELECTRON_MASS);//not sure. why cant I use pointer here?
    math::PtEtaPhiMLorentzVector gsf_e_lv(gsf_electron.pt(),gsf_electron.eta(), gsf_electron.phi(), ELECTRON_MASS);//particle physics: review LorentzVector
    math::PtEtaPhiMLorentzVector zlv;
    zlv = hf_e_lv + gsf_e_lv;
    double m = zlv.mass();//Z mass
    // formula from Perrie Cole's thesis, page 16
    double numerator = m*  m*             cosh(  gsf_electron.eta()  )*     cosh(  hf_electron.eta()  );
    double coshTerm = cosh(     gsf_electron.eta()  -      hf_electron.eta() );
    double cosTerm = cos(  reco::deltaPhi( gsf_electron.phi()  ,  hf_electron.phi())    );
    double E_Ecal = gsf_electron.p();
    double denominator = 2.0 * E_Ecal  * (  coshTerm - cosTerm  );
    double pred_hf_e_energy = numerator / denominator;
    return pred_hf_e_energy;
}
   
// ------------ method called for each event  ------------
void
calib::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)


{
    using namespace edm;
    //https://github.com/UMN-CMS/ZFinder/blob/master/ZFinder/ZSkimmer/src/ZSkimmer.cc
    
   bool GSF_LOOSE_MEDIUM_TIGHT=0;//SCOPE VARIABLE 
    
    //THE HF BLOCK---------------------------------------  
    reco::RecoEcalCandidate high_pt_hf_electron;//VARIABLE OUTSIDE OF FOR LOOP
    double highest_hf_pt = 0.;

    edm::Handle<reco::RecoEcalCandidateCollection> hf_electrons;//does <argument> means using a template ?
    iEvent.getByLabel("hfRecoEcalCandidate", hf_electrons);

    //VARIABLES USED FOR SORTING HF ELECTRONS - 
    // HF Superclusters
    edm::Handle<reco::SuperClusterCollection> scs_h;//what's super cluster --> different parts of the Detector
    //iEvent.getByLabel(inputtags_.hf_clusters, scs_h);//change name
    iEvent.getByLabel("hfEMClusters", scs_h);
    edm::Handle<reco::HFEMClusterShapeAssociationCollection> scas_h;
    //iEvent.getByLabel(inputtags_.hf_clusters, scas_h);
    iEvent.getByLabel("hfEMClusters", scas_h);//change name to label


    for(unsigned int i = 0; i < hf_electrons->size(); ++i) {
        reco::RecoEcalCandidate hf_electron = hf_electrons->at(i);
        //  math::PtEtaPhiMLorentzVector hf_e_lv(hf_electron.pt(), hf_electron.eta(), hf_electron.phi(), ELECTRON_MASS);//temp
        //COUNTER
        all_hf_e_histo_ ->Fill(hf_electron.pt());
       // all_events_histo_->Fill(4); 
      
        if ((hf_electron.pt() >= highest_hf_pt)&&(hf_electron.pt()>20)) {//carried from Alex code. //may be lower is the future
            highest_hf_pt = hf_electron.pt();
            high_pt_hf_electron = hf_electron;

            reco::SuperClusterRef cluster_ref = hf_electron.superCluster();
            const reco::HFEMClusterShapeRef CLUSTER_SHAPE_REF = scas_h->find(cluster_ref)->val;
            const reco::HFEMClusterShape& CLUSTER_SHAPE = *CLUSTER_SHAPE_REF;
            const double ECE9 = CLUSTER_SHAPE.eCOREe9();
            const double ESEL = CLUSTER_SHAPE.eSeL();
            //unused now. may be later//const double E9E25 = (CLUSTER_SHAPE.eLong3x3() * 1.0 / CLUSTER_SHAPE.eLong5x5());
            // e9e25 cut
            //unused now--//const bool PASS_E9E25 = (E9E25 > 0.94); //WHAT ARE THESE NUMBER ?
            // HF Tight (as defined in hfRecoEcalCandidate_cfi.py in ZShape)
            //const double TIGHT2D = (ECE9 - (ESEL * 0.20));
            //std::cout << TIGHT2D <<std::endl;
            //const bool HFTIGHT = (TIGHT2D > 0.92);
            // HF Medium
            //unused now//const double MEDIUM2D = (ECE9 - (ESEL * 0.275));
            //unused now//const bool HFMEDIUM = (MEDIUM2D > 0.875);
            // HF Loose
            const double LOOSE2D = (ECE9 - (ESEL * 0.475));//slope
            std::cout <<"Loose2D: " <<LOOSE2D <<std::endl;
            const bool HFLOOSE = (LOOSE2D > 0.815);//these numbers need to retune for 13T

            //adding an action to test output here  
            if (HFLOOSE) {std::cout << "HFLOOSE is true  ";}
        }
    }
    
    //THE GSF BLOCK---------------------------------------------------------------------------------------------------
    // conversions
    edm::Handle<reco::ConversionCollection> conversions_h;
    //iEvent.getByLabel(inputtags_.conversion, conversions_h);
    iEvent.getByLabel("allConversions", conversions_h);

    //   turned these on for 13 TeV 
    typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;
    IsoDepositVals isoVals(3); 
    iEvent.getByLabel("elPFIsoValueCharged03PFIdPFIso", isoVals[0]);//OK ??? charge
    iEvent.getByLabel("elPFIsoValueGamma03PFIdPFIso", isoVals[1]);//OK ??? photon
    iEvent.getByLabel("elPFIsoValueNeutral03PFIdPFIso", isoVals[2]);//OK ??? neutral hadron

    // beam spot
    edm::Handle<reco::BeamSpot> beamspot_h;
    //    iEvent.getByLabel(inputtags_.beamspot, beamspot_h);
    iEvent.getByLabel("offlineBeamSpot", beamspot_h);
    const reco::BeamSpot &beamSpot = *(beamspot_h.product());

    // vertices
    edm::Handle<reco::VertexCollection> vtx_h;
    //    iEvent.getByLabel(inputtags_.vertex, vtx_h);
    iEvent.getByLabel("offlinePrimaryVertices", vtx_h);

    // rho for isolation
    // The python uses:
    // cms.InputTag("kt6PFJetsForIsolation", "rho")
    edm::Handle<double> rho_iso_h;
    iEvent.getByLabel("kt6PFJetsForIsolation", "rho", rho_iso_h);
    const double RHO_ISO = *(rho_iso_h.product());

    reco::GsfElectron high_pt_gsf_electron;//VARIABLE FOR USE OUTSIDE OF FOR LOOP
    double highest_gsf_pt = 0.;

    //SELECTING HIGHEST PT GSF ELECTRON    
    edm::Handle<reco::GsfElectronCollection> gsf_electrons;
    //iEvent.getByLabel("gsfElectrons", gsf_electrons);
    iEvent.getByLabel("gedGsfElectrons", gsf_electrons);//the MC of 13Tev doesn't have label gsfElectrons
    for(unsigned int i = 0; i < gsf_electrons->size(); ++i) {
        reco::GsfElectron gsf_electron = gsf_electrons->at(i); 
        //math::PtEtaPhiMLorentzVector gsf_e_lv(gsf_electron.pt(), gsf_electron.eta(), gsf_electron.phi(), ELECTRON_MASS);//temp
        
        //COUNTER
        all_gsf_e_histo_ ->Fill(gsf_electron.pt());
        //all_events_histo_->Fill(3); 
        
        if ((gsf_electron.pt() >= highest_gsf_pt)&&(gsf_electron.pt()>20)) { ///excercise 2
            highest_gsf_pt = gsf_electron.pt();    
            high_pt_gsf_electron = gsf_electron;//THIS MIGHT HAVE TO MOVE DOWN AFTER THE ADDITIONAL CUT ?

            // ADDITIONAL CUTS (ISO, SHAPE, R9 )GONNA BE HERE; code stolen from Alex:

            

            // get reference to electron and the electron, from the collection
            reco::GsfElectronRef ele_ref(gsf_electrons, i);
            // get particle flow isolation
            const double ISO_CH = (*(isoVals)[0])[ele_ref];
            const double ISO_EM = (*(isoVals)[1])[ele_ref];
            const double ISO_NH = (*(isoVals)[2])[ele_ref];
            //const double ISO_CH = 0;
           //const double ISO_EM = 0;
            //const double ISO_NH = 0;

            // test ID
            // working points // veto, loose, medium are unused now
            const bool VETO = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::VETO, ele_ref, conversions_h, beamSpot, vtx_h, ISO_CH, ISO_EM, ISO_NH, RHO_ISO, ElectronEffectiveArea::kEleEANoCorr);
            const bool LOOSE = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::LOOSE, ele_ref, conversions_h, beamSpot, vtx_h, ISO_CH, ISO_EM, ISO_NH, RHO_ISO, ElectronEffectiveArea::kEleEANoCorr);
            const bool MEDIUM = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::MEDIUM, ele_ref, conversions_h, beamSpot, vtx_h, ISO_CH, ISO_EM, ISO_NH, RHO_ISO, ElectronEffectiveArea::kEleEANoCorr);
            const bool TIGHT = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::TIGHT, ele_ref, conversions_h, beamSpot, vtx_h, ISO_CH, ISO_EM, ISO_NH, RHO_ISO, ElectronEffectiveArea::kEleEANoCorr);
            if (MEDIUM) {std::cout << "MEDIUM" << std::endl;}
            if (LOOSE) {std::cout << "LOOSE " << std::endl;}//use this line to avoid unused variable
            if (TIGHT) {std::cout << "TIGHT" << std::endl;}
            if (VETO) {std::cout << "VETO" << std::endl;}
            
            
           GSF_LOOSE_MEDIUM_TIGHT = 1;
         
        }
    }

    //GENERATOR PARTICLE ---------------------------------------------------------------------------------
    //: from Alex's initTruth https://github.com/UMN-CMS/ZFinder/blob/master/ZFinder/Event/src/ZFinderEvent.cc
    /* Alex:
     * We don't need to select electrons with cuts, because in Monte Carlo we
     * can just ask for the Z.
     */
    edm::Handle<reco::GenParticleCollection> mc_particles;
    iEvent.getByLabel("genParticles", mc_particles);

    /* Alex: Finding the Z and daughter electrons
     *
     * We loop over all gen particles. If it is a Z, we check its daughters
     * until we find an electron, then we know that it is a Z->ee decay. If
     * this is the first Z we save it. If the particle is an electron, we
     * make sure it came from a Z. This might have problems in ZZ->eeee
     * decays, but we expect those to be impossibly rare.
     */
    const reco::GenParticle* electron_0 = NULL;//used for outside the loop
    const reco::GenParticle* electron_1 = NULL;//are these reset for each event?
    const reco::GenParticle* z_boson = NULL;  
    for(unsigned int i = 0; i < mc_particles->size(); ++i) {
        const reco::GenParticle* gen_particle = &mc_particles->at(i);//loop over generator particles. syntax ? don't understand
        // Is a Z
        if (gen_particle->pdgId() == 23 && z_boson == NULL) { //23: Z zero. Why z_boson ==NULL ? not sure.
            for (size_t j = 0; j < gen_particle->numberOfDaughters(); ++j) {//loop over daughters
                if (fabs(gen_particle->daughter(j)->pdgId() == 11)) {
                    z_boson = gen_particle;//save the first Z why? --> i won't have ZZ events 1%
                    
                    //COUNTING EVENTS - MESSY TEMPORARY CODE
                    //-----------------------------------
                    
                     // count generator Z
                    if ((z_boson->mass() >= 60.0) && (z_boson->mass() <= 120.0 ) )
                    {
                        all_events_histo_->Fill(10);
                       /* for(unsigned int i = 0; i < gsf_electrons->size(); ++i) {  //count GSF if Generator Z exist.
                            if (highest_gsf_pt >0)
                            {
                                    all_events_histo_->Fill(3); 
                                    for(unsigned int i = 0; i < hf_electrons->size(); ++i) { //count HF only if GSF was found
                                        if (highest_hf_pt > 0
                                            && (fabs(high_pt_hf_electron.eta()) >=3.1) && (fabs(high_pt_hf_electron.eta())<=4.6)
                                            )
                                            {all_events_histo_->Fill(6); }//if(highest_hf_pt > 0)
             							
             							}//  for(unsigned int i = 0; i < hf_electrons->size(); ++i) 
                          
                            }//if(highest_gsf_pt >0)
                          
                          }//for(unsigned int i = 0; i < gsf_electrons->size(); ++i) 
                       */
                     
                     } //if ((z_boson->mass() >= 60.0) && (z_boson->mass() <= 120.0 ))                   
                     //----------------------------------
                     
                    break;
                }//f (fabs(gen_particle->daughter(j)->pdgId() == 11)) 
            }//for (size_t j = 0; j < gen_particle->numberOfDaughters(); ++j)
            
            // Is an e+ / e-. //
        }//if (gen_particle->pdgId() == 23 && z_boson == NULL) 
        else if (fabs(gen_particle->pdgId()) == 11// In pdgId, fabs(POSITRON) == ELECTRON
                && (electron_0 == NULL || electron_1 == NULL)// verify this line
                ) {
            	for (size_t j = 0; j < gen_particle->numberOfMothers(); ++j) {//loop over mothers
                	if (gen_particle->mother(j)->pdgId() == 23) {
                    	if (electron_0 == NULL) {
                        electron_0 = gen_particle;}
                    	else { electron_1 = gen_particle;}
                	}
                	
                	
                	
            }
        }//end of else if (electron)
    }//the end of generator loop to find z
    
    // Continue only if all particles have been found
    if (z_boson != NULL && electron_0 != NULL && electron_1 != NULL) {
        // We set electron_0 to the higher pt electron
        if (electron_0->pt() < electron_1->pt()) {std::swap(electron_0, electron_1);} 
    } 
    else {
        return;
    }

    //COMMON ACTIONS ON HF AND GSF ELECTRONS------------------------------------------------------------------------
    //Mass calculation https://github.com/UMN-CMS/ZFinder/blob/master/ZFinder/Event/src/ZFinderEvent.cc: 467
    const double ELECTRON_MASS = 5.109989e-4;//
    math::PtEtaPhiMLorentzVector hf_e_lv(high_pt_hf_electron.pt(), high_pt_hf_electron.eta(), high_pt_hf_electron.phi(), ELECTRON_MASS);//not sure. why cant I use pointer here?
    math::PtEtaPhiMLorentzVector gsf_e_lv(high_pt_gsf_electron.pt(),high_pt_gsf_electron.eta(), high_pt_gsf_electron.phi(), ELECTRON_MASS);//not sure. review LorentzVector
    math::PtEtaPhiMLorentzVector zlv;
    zlv = hf_e_lv + gsf_e_lv;
    double m = zlv.mass();//Z mass
   //the above block of code is also tower_plotter. should clean up.


    double PredHfElecEnergy = pred_hf_e_energy(high_pt_gsf_electron, high_pt_hf_electron);
    double RecoOverPredHfElecEnergy = high_pt_hf_electron.p()/PredHfElecEnergy;

    
    //testing the GenConverter. 
    Genconverter* hf_e_ieta_gen_convert=NULL;
    Genconverter* hf_e_iphi_gen_convert=NULL;
    
    //GETTING IETA AND IPHI TO PASS ON TO TOWERPLOTTER OBJ
    int hf_e_ieta = hf_e_ieta_gen_convert->Eta2IEta(high_pt_hf_electron.eta());
    int hf_e_iphi = hf_e_iphi_gen_convert->Phi2Iphi(high_pt_hf_electron.phi(),high_pt_hf_electron.eta());
    
    //HISTOGRAM OUTPUT AND MASS PRINT OUT
    
    //COUNTER BLOCK
    
   all_events_histo_->Fill(1); //count all event process
   
   
 //Counting maximum generator electron
    if  (
            (
            ( ( fabs((electron_0->eta()))  >= 3.1 ) && ( fabs(electron_0->eta()) <=4.6) && (fabs(electron_0->pt()) >=20) )
            || ( ( fabs(electron_1->eta())  >= 3.1 ) && ( fabs(electron_1->eta()) <=4.6) && fabs(electron_1->pt()>=20))
            )
            && (  ( z_boson->mass() >= 60.0 ) && ( z_boson->mass() <= 120.0 ) )   
        )
        {
            if (
                (fabs(high_pt_hf_electron.eta()) >=3.1) 
                && (fabs(high_pt_hf_electron.eta())<=4.6)
                && (highest_hf_pt > 0)
                ) {all_events_histo_->Fill(6);}//counting reco hf
                
            all_events_histo_->Fill(7);//counting gen hf
        }
        
    if  (
        ( ( fabs(electron_0->eta())  >= -1.0 ) && ( fabs(electron_0->eta()) <=2.5) && fabs(electron_0->pt()>=20)) //no lower bound
        || ( ( fabs(electron_1->eta())  >= -1.0 ) && ( fabs(electron_1->eta()) <=2.5) && fabs(electron_1->pt()>=20))
        
        )
        { 
            if (highest_gsf_pt >0) {all_events_histo_->Fill(3);}//counting reco gsf
            all_events_histo_->Fill(4);//counting gen gsf
        }
        
   
 //GEN PARTICLE HISTO
        mass_gen_Z_histo_->Fill(z_boson->mass());
        eta_gen_Z_histo_->Fill(z_boson->eta());			
        phi_gen_Z_histo_->Fill(z_boson->phi());	
        
        pt_gen_e1_histo_->Fill(electron_1->pt());
	    eta_gen_e1_histo_->Fill(electron_1->eta());
	    phi_gen_e1_histo_->Fill(electron_1->phi());
	  
	    pt_gen_e0_histo_->Fill(electron_0->pt());
	    eta_gen_e0_histo_->Fill(electron_0->eta());
	    phi_gen_e0_histo_->Fill(electron_0->phi());
    // A BUNCH OF CUTS
    if (    (highest_hf_pt > 20.0)&&(highest_gsf_pt >20.0) //not sure. greater or equal or just greater ? is it possible to have pt=0 
    		
    		&& GSF_LOOSE_MEDIUM_TIGHT  // LOOSE OR MEDIUM OR TIGHT, SET ABOVE
    		
    		&& ( //Jeremy It would be very useful to have a requirement that there is a _generator_ electron with |eta| in HF.

    		        ((fabs(electron_0->eta()) >= 3.1) && (fabs(electron_0->eta()) <=4.6 )) 
    		        || ((fabs(electron_1->eta()) >= 3.1) && (fabs(electron_1->eta()) <=4.6 ))
    		    )  //HAVING A GENERATOR e in HF with the right eta
    		    
    		&& (z_boson->mass() >= 60.0) && (z_boson->mass() <= 120.0 ) // GENERATOR Z MASS 
    		
    		&& (fabs(high_pt_hf_electron.eta()) >=3.1) && (fabs(high_pt_hf_electron.eta())<=4.6)  //RECO HF ELECTRON ETA RANGE
    		
    		//&& (fabs(high_pt_gsf_electron.eta())>=1.56) && (fabs(high_pt_gsf_electron.eta())<=2.85) // no eta cut for gsf in any case.(Jeremiah)
    		
    	)   
    {   
       
       std::cout << "found highest transverse momentum HF electron: "<< highest_hf_pt << " and GSF electron: " << highest_gsf_pt<<"  " << std::endl;  
       std::cout <<"test cosh of Pi/3:" <<cosh(3.14/3.0) <<std::endl;
				
       std::cout << "predicted energy of HF electron: " << PredHfElecEnergy <<std::endl;
       std::cout << "RECO p of HF electron: " << high_pt_hf_electron.p() << std:: endl;
       std::cout << "iEta: " << hf_e_ieta << std::endl;
       std::cout << "iPhi: " << hf_e_iphi << std ::endl;
       std::cout << "------------------------------------" << std::endl;
        //z mass calculated from highest pt e

        mass_histo_->Fill(m);
       all_events_histo_->Fill(9); // should this be moved inside the 

        //HF HISTO
        pt_hf_e_histo_->Fill(high_pt_hf_electron.pt());
        eta_hf_e_histo_->Fill(high_pt_hf_electron.eta());
        phi_hf_e_histo_->Fill(high_pt_hf_electron.phi());
        
        //predicted energy of hf e
        pred_hf_e_energy_histo_->Fill(PredHfElecEnergy);//(Perrie thesis formula)
        reco_over_pred_hf_e_energy_histo_->Fill(RecoOverPredHfElecEnergy);
       
          
        //GSF HISTO			
        pt_gsf_e_histo_->Fill(high_pt_gsf_electron.pt());
        eta_gsf_e_histo_->Fill(high_pt_gsf_electron.eta());
        phi_gsf_e_histo_->Fill(high_pt_gsf_electron.phi());
        
        
	    
	    //Testing towerplot
	    //towerplot make histos groupd index by ieta and iphi
	    tower_plotter_obj_->Fill(hf_e_ieta,hf_e_iphi,high_pt_hf_electron,high_pt_gsf_electron);
	    
	    
                    
    }//closing the if

}//closing the member function: calib
//}//missing one } ?

// ------------ method called once each job just before starting event loop  ------------
void calib::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void calib::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
calib::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
calib::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
calib::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
calib::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
calib::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(calib);
