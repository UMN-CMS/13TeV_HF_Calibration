#include "z_plots.h"
#include "DataFormats/Math/interface/deltaPhi.h" 

// You'll need to include the 4-vector files

ZPlots::ZPlots(TFileDirectory& tdir) {//when to use passing by ref, when to use pointer ? // take directory as argument, doesn't make one
    // Make plots like this:
    z_mass_histo_ = tdir.make<TH1D>("z_mass", "Z Mass", 100., 50., 150.);
    z_mass_histo_->GetXaxis()->SetTitle("m_{ee} [GeV]");
    z_mass_histo_->GetYaxis()->SetTitle("Counts / GeV");//why counts /GeV ?

   //Z BOSON
    z_pt_histo_ = tdir.make<TH1D>("z_pt","Z P_t", 100., 0, 150.);//the range need to be change
    z_pt_histo_->GetXaxis()->SetTitle("Z P_{t} [GeV]");
    z_pt_histo_->GetYaxis()->SetTitle("Counts / GeV");
    
    reco_over_pred_hf_e_energy_histo_ = tdir.make<TH1D>("Reco/Predicted_energy_hf_e_histo", "Reco/Predicted_energy_{hf_e}", 100, 0, 2);
    
    //HF
    hf_e_pt_histo_ = tdir.make<TH1D>("pt_hf_e_histo", "Pt_{hf_e}", 100, 0, 100);
	hf_e_eta_histo_ = tdir.make<TH1D>("eta_hf_e_histo", "eta_{hf_e}", 100, 2.5, 5);
	hf_e_phi_histo_ = tdir.make<TH1D>("phi_hf_e_histo", "phi_{hf_e}", 100, -3.2, 3.2);//radian
	//GSF
	gsf_e_pt_histo_ = tdir.make<TH1D>("pt_gsf_e_histo", "Pt_{gsf_e}", 100, 0, 100);
  	gsf_e_eta_histo_ = tdir.make<TH1D>("eta_gsf_e_histo", "eta_{gsf_e}", 100, 0, 3);
  	gsf_e_phi_histo_ = tdir.make<TH1D>("phi_gsf_e_histo", "phi_{gsf_e}", 100, -3.2, 3.2);
}

void ZPlots::Fill(

    const reco::RecoEcalCandidate& hf_electron,
    const reco::GsfElectron& gsf_electron
) {
    // Make the Z using 4-vectors
    
    const double ELECTRON_MASS = 5.109989e-4;
    math::PtEtaPhiMLorentzVector hf_e_lv(hf_electron.pt(), hf_electron.eta(), hf_electron.phi(), ELECTRON_MASS);
    math::PtEtaPhiMLorentzVector gsf_e_lv(gsf_electron.pt(),gsf_electron.eta(), gsf_electron.phi(), ELECTRON_MASS);
    math::PtEtaPhiMLorentzVector zlv;
    zlv = hf_e_lv + gsf_e_lv;
    double z_mass = zlv.mass();//Z mass
    double z_pt = zlv.pt();
    
    //calculate predicted energy
    // formula from Perrie Cole's thesis, page 16
    double numerator = z_mass*  z_mass*             cosh(  gsf_electron.eta()  )*     cosh(  hf_electron.eta()  );
    double coshTerm = cosh(     gsf_electron.eta()  -      hf_electron.eta() );
    double cosTerm = cos(  reco::deltaPhi( gsf_electron.phi()  ,  hf_electron.phi())    );
    double E_Ecal = gsf_electron.p();
    double denominator = 2.0 * E_Ecal  * (  coshTerm - cosTerm  );
    double pred_hf_e_energy = numerator / denominator;
    double RecoOverPredHfElecEnergy = hf_electron.p()/pred_hf_e_energy;
    
    // Fill the plots
    //Z boson
    z_mass_histo_->Fill(z_mass);
    z_pt_histo_->Fill(z_pt);
    reco_over_pred_hf_e_energy_histo_->Fill(RecoOverPredHfElecEnergy);
    
     //HF HISTO
    hf_e_pt_histo_->Fill(hf_electron.pt());
    hf_e_eta_histo_->Fill(hf_electron.eta());
    hf_e_phi_histo_->Fill(hf_electron.phi());
        
    //GSF HISTO			
    gsf_e_pt_histo_->Fill(gsf_electron.pt());
    gsf_e_eta_histo_->Fill(gsf_electron.eta());
    gsf_e_phi_histo_->Fill(gsf_electron.phi());
}
