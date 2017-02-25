#ifndef MVAVARIABLES_H
#define MVAVARIABLES_H

#include <Rtypes.h> // Double_t

class TTree;
class TBranch;

namespace tthMEM
{
  /**
   * @brief Simple class for storing MVA input-output variables
   *
   * The class is implemented only to copy the variables from
   * an input tree to a new tree. The variables stored here are
   * not used in MEM in any way. Therefore no getters or setters
   * are defined here as the user is expected to not use the class
   * for anything other than copying the data from one tree to another.
   */
  class MVAVariables
  {
  public:
    MVAVariables() = default;

    void
    setBranches(TTree * t);

    void
    initNewBranches(TTree * t);

  private:
    Double_t max2LeptonEta;
    Double_t MT_met_lep1;
    Double_t nJet25_Recl;
    Double_t mindr_lep1_jet;
    Double_t mindr_lep2_jet;
    Double_t lep1_cone_pt;
    Double_t lep3_cone_pt;
    Double_t avg_dr_jet;
    Double_t mhtJet25_Recl;
    Double_t mvaOutput_3l_ttV;
    Double_t mvaOutput_3l_ttbar;

    TBranch * branch_max2LeptonEta      = nullptr;
    TBranch * branch_MT_met_lep1        = nullptr;
    TBranch * branch_nJet25_Recl        = nullptr;
    TBranch * branch_mindr_lep1_jet     = nullptr;
    TBranch * branch_mindr_lep2_jet     = nullptr;
    TBranch * branch_lep1_cone_pt       = nullptr;
    TBranch * branch_lep3_cone_pt       = nullptr;
    TBranch * branch_avg_dr_jet         = nullptr;
    TBranch * branch_mhtJet25_Recl      = nullptr;
    TBranch * branch_mvaOutput_3l_ttV   = nullptr;
    TBranch * branch_mvaOutput_3l_ttbar = nullptr;
  };
}

#endif // MVAVARIABLES_H
