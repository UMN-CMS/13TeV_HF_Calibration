#include "z_plots.h"

// You'll need to include the 4-vector files

ZPlots::ZPlots(TFileDirectory& tdir) {//take a dir as argument
    // Make plots like this:
    z_mass_ = tdir.make<TH1D>("z_mass", "Z Mass", 100, 50., 150.);
    z_mass_->GetXaxis()->SetTitle("m_{ee} [GeV]");
    z_mass_->GetYaxis()->SetTitle("Counts / GeV");

    // Fill in the rest of the plots

}

void ZPlots::Fill(
    const type& hf_electron,
    const type& gsf_electron
) {
    // Make the Z using 4-vectors

    // Fill the plots
    z_mass_->Fill(mass);
}
