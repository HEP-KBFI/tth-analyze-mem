#include "tthAnalysis/tthMEM/interface/MEMOutput_3l1tau.h"

#include <iomanip> // std::setprecision()

MEMOutput_3l1tau::MEMOutput_3l1tau()
  : prob_tth(-1.)
  , prob_tth_err(-1.)
  , prob_ttz(-1.)
  , prob_ttz_err(-1.)
  , prob_tth_h2ww(-1.)
  , prob_tth_h2ww_err(-1.)
  , sig(-1.)
  , sig_err(-1.)
  , sig_up(-1.)
  , sig_down(-1.)
  , bkg(-1.)
  , bkg_err(-1.)
  , bkg_up(-1.)
  , bkg_down(-1.)
  , lr(-1.)
  , lr_err(-1.)
  , lr_up(-1.)
  , lr_down(-1.)
{}

std::ostream &
operator<<(std::ostream & os,
           const MEMOutput_3l1tau & output)
{
  std::ios_base::fmtflags prev_flags(os.flags());
  os << std::scientific << std::setprecision(4)
     << "P(tth,h->tautau)              = " << output.prob_tth      << " (" << output.prob_tth_err      << ")\n"
     << "P(tth,h->WW)                  = " << output.prob_tth_h2ww << " (" << output.prob_tth_h2ww_err << ")\n"
     << "P(ttz,z->tautau)              = " << output.prob_ttz      << " (" << output.prob_ttz_err      << ")\n"
     << "P(signal hypothesis) weighted = " << output.sig           << " (" << output.sig_err           << ")\n"
     << "P(bkg hypothesis) weighted    = " << output.bkg           << " (" << output.bkg_err           << ")\n"
     << std::fixed << std::setprecision(6)
     << "Likelihood ratio              = " << output.lr            << " (" << output.lr_err            << ")\n";
  os.flags(prev_flags);
  return os;
}
