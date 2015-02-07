#include "tower_plotter.h"

// Standard Library
#include <string>  // std::string
#include <sstream>
//cmssw
#include "FWCore/ServiceRegistry/interface/Service.h" // edm::Service
#include "CommonTools/UtilAlgos/interface/TFileService.h" // TFileService

// ZPlots
#include "z_plots.h"  // ZPlots

TowerPlotter::TowerPlotter(TFileDirectory& tdir) {//take dir as argument
    // Set up internal tdir_
    tdir_ = tdir.mkdir("TowerPlotter");//create a subfolder under tdir, which will survive with the towerplooter obj (one way to avoid segfault)
}

void TowerPlotter::Fill(//take in ieta, iphi, 2 electrons to do plot and set the folder name correspondingly
    const int ieta,
    const int iphi,
    const reco::RecoEcalCandidate& hf_electron,
    const reco::GsfElectron& gsf_electron
) {
    //edm::Service<TFileService> fs;//trying to fix seg fault
    
    // The pair to use in the map
    const std::pair<int, int> ipair = std::make_pair(ieta, iphi);

    // We try to get the entry. If we fail then we know we have to make the
    // ZPlots object first, otherwise we just fill it.
    auto iter = zplot_map_.find(ipair);//Quynh: auto: the type is deduced from variable initializer

    // In this case the ZPlots must be created first
    if (iter == zplot_map_.end()) {//Quynh: map::end Returns an iterator referring to the past-the-end element in the map container.
    
        // cook up The Name of the directory
       //const std::string name = "ZPlots: iEta: " + ieta + ", iPhi: " + iphi; //Alex's: this gave error
       std::string ieta_string;//Convert ieta int to string , lengthy but works
       std::ostringstream convert_ieta;
       convert_ieta << ieta;
       ieta_string = convert_ieta.str();
       
       std::string iphi_string;//Convert iphi int to string lengthy but works
       std::ostringstream convert_iphi;
       convert_iphi << iphi;
       iphi_string = convert_iphi.str(); 
       
       const std::string name = "ZPlots iEta " + ieta_string + " iPhi " + iphi_string; 
        
        TFileDirectory zp_dir = tdir_.mkdir(name.c_str());//Quynh: create zp_dir (with name indicating ieta, iphi) as a subdir of tdir_
        
        // Make a ZPlots and add them to the map
        ZPlots zplots(zp_dir);//Zplots is a class, zplots is the object initiated
        
        //zplot_map_.insert(ipair, zplots);//Alex's- error: no matching function
        zplot_map_.insert(std::make_pair(ipair, zplots));//Quynh: fixed it
        
        // Move the iterator to the newly created entry
        iter = zplot_map_.find(ipair);
      
    }

    // Fill the plots
    //iter.second.Fill(hf_electron, gsf_electron);//struct std::_Rb_tree_iterator<std::pair<const std::pair<int, int>, ZPlots> >' has no member named 'second'
        iter->second.Fill(hf_electron, gsf_electron);//use pointer to second. "second" refers to the value of the map (vs. the key)

}
