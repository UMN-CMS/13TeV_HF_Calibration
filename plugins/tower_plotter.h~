#ifndef TOWER_PLOTTER_H_
#define TOWER_PLOTTER_H_

// Standard Library
#include <map>  // std::map 
#include <utility>  // std::pair

// CMSSW
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// ZPlots
#include "z_plots.h"  // ZPlots


class TowerPlotter {
    public:
        // Constructor
        TowerPlotter(TFileDirectory& tdir);//* or & ?

        // Add events
        void Fill(
            const int ieta,
            const int iphi,
            const reco::RecoEcalCandidate& hf_electron,
            const reco::GsfElectron& gsf_electron
        );

    protected:
        // The map
        std::map<std::pair<int, int>, ZPlots> zplot_map_;// mapping between ieta, iphi as key and value is zplots

        // Internal TFileDirectory
        TFileDirectory tdir_;
};
#endif  // TOWER_PLOTTER_H_
