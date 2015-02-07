#ifndef Z_PLOTS_H_
#define Z_PLOTS_H_

// Root
#include <TH1D.h>  // TH1D

// CMSSW
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Include the necessary files for hf and gsf electrons
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h" // GsfElectron
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h" // reco::RecoEcalCandidate


class ZPlots {
    public:
        // Constructor
        ZPlots(TFileDirectory& tdir);

        // Add events
        void Fill(
            const reco::RecoEcalCandidate& hf_electron,
            const reco::GsfElectron& gsf_electron
        );

    protected:
        // Histograms
        TH1D* z_mass_histo_;
        TH1D* z_pt_histo_;//? and add predicted energy ratio
        TH1D* reco_over_pred_hf_e_energy_histo_;
     
        
        TH1D* gsf_e_eta_histo_;
        TH1D* gsf_e_phi_histo_;
        TH1D* gsf_e_pt_histo_;
        TH1D* hf_e_eta_histo_;
        TH1D* hf_e_phi_histo_;
        TH1D* hf_e_pt_histo_;
 /*        TH1D* z_rapidity_;//not very important
*/        
};
#endif  // Z_PLOTS_H_
