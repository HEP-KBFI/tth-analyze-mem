#include "tthAnalysis/tthMEM/interface/objects/MVAVariables.h" // tthMEM::MVAVariables

using namespace tthMEM;

void
MVAVariables::setBranches(TTree * t)
{
  t -> SetBranchAddress("max2LeptonEta",      &max2LeptonEta);
  t -> SetBranchAddress("MT_met_lep1",        &MT_met_lep1);
  t -> SetBranchAddress("nJet25_Recl",        &nJet25_Recl);
  t -> SetBranchAddress("mindr_lep1_jet",     &mindr_lep1_jet);
  t -> SetBranchAddress("mindr_lep2_jet",     &mindr_lep2_jet);
  t -> SetBranchAddress("lep1_cone_pt",       &lep1_cone_pt);
  t -> SetBranchAddress("lep3_cone_pt",       &lep3_cone_pt);
  t -> SetBranchAddress("avg_dr_jet",         &avg_dr_jet);
  t -> SetBranchAddress("mhtJet25_Recl",      &mhtJet25_Recl);
  t -> SetBranchAddress("mvaOutput_3l_ttV",   &mvaOutput_3l_ttV);
  t -> SetBranchAddress("mvaOutput_3l_ttbar", &mvaOutput_3l_ttbar);
}

void
MVAVariables::initNewBranches(TTree * t)
{
  branch_max2LeptonEta      = t -> Branch("max2LeptonEta",
                                          &max2LeptonEta,      "max2LeptonEta/D");
  branch_MT_met_lep1        = t -> Branch("MT_met_lep1",
                                          &MT_met_lep1,        "MT_met_lep1/D");
  branch_nJet25_Recl        = t -> Branch("nJet25_Recl",
                                          &nJet25_Recl,        "nJet25_Recl/D");
  branch_mindr_lep1_jet     = t -> Branch("mindr_lep1_jet",
                                          &mindr_lep1_jet,     "mindr_lep1_jet/D");
  branch_mindr_lep2_jet     = t -> Branch("mindr_lep2_jet",
                                          &mindr_lep2_jet,     "mindr_lep2_jet/D");
  branch_lep1_cone_pt       = t -> Branch("lep1_cone_pt",
                                          &lep1_cone_pt,       "lep1_cone_pt/D");
  branch_lep3_cone_pt       = t -> Branch("lep3_cone_pt",
                                          &lep3_cone_pt,       "lep3_cone_pt/D");
  branch_avg_dr_jet         = t -> Branch("avg_dr_jet",
                                          &avg_dr_jet,         "avg_dr_jet/D");
  branch_mhtJet25_Recl      = t -> Branch("mhtJet25_Recl",
                                          &mhtJet25_Recl,      "mhtJet25_Recl/D");
  branch_mvaOutput_3l_ttV   = t -> Branch("mvaOutput_3l_ttV",
                                          &mvaOutput_3l_ttV,   "mvaOutput_3l_ttV/D");
  branch_mvaOutput_3l_ttbar = t -> Branch("mvaOutput_3l_ttbar",
                                          &mvaOutput_3l_ttbar, "mvaOutput_3l_ttbar/D");
}

