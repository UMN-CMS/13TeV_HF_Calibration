#ifndef Z_PLOTS_H_
#define Z_PLOTS_H_

// Root
#include <TH1D.h>  // TH1D

// CMSSW
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Include the necessary files for hf and gsf electrons

class ZPlots {
    public:
        // Constructor
        ZPlots(TFileDirectory& tdir);

        // Add events
        void Fill(
            const type& hf_electron,
            const type& gsf_electron
        );

    protected:
        // Histograms
        TH1D* z_mass_;
        TH1D* z_rapidity_;
        TH1D* z_pt_;
        TH1D* gsf_e_eta_;
        TH1D* gsf_e_phi_;
        TH1D* gsf_e_pt_;
        TH1D* hf_e_eta_;
        TH1D* hf_e_phi_;
        TH1D* hf_e_pt_;
};
#endif  // Z_PLOTS_H_
