#include "tthAnalysis/tthMEM/interface/MEMInterface_3l1tau.h" // MEMInterface_3l1tau
#include "tthAnalysis/tthMEM/interface/MeasuredEvent_3l1tau.h" // MeasuredEvent_3l1tau

MEMInterface_3l1tau::MEMInterface_3l1tau()
  : pdfName                ("MSTW2008lo68cl")
  , madgraphFileName       ("tthAnalysis/tthMEM/data/param_card.dat")
  , integrationMode        ("markovchain")
  , maxObjFunctionCalls    (100000)
  , higgsWidth             (-1.)
  , useBJetTransferFunction(true)
  , useAvgBjetCombo        (true)
  , mxMode                 ("uniform")
  , nofBatches             (100)
  , nofChains              (1)
  , maxCallsStartingPos    (1000)
  , epsilon0               (1.e-2)
  , T0                     (15.)
  , nu                     (0.71)
  , err                    (0)
  , mem_tt_HandZ           (pdfName, findFile(madgraphFileName), vm)
{}

void
MEMInterface_3l1tau::initialize()
{
  try
  {
    mem_tt_HandZ.setIntegrationMode(integrationMode);
    mem_tt_HandZ.setMaxObjFunctionCalls(maxObjFunctionCalls);
    mem_tt_HandZ.setBJetTransferFunction(useBJetTransferFunction);
    mem_tt_HandZ.useAvgBjetCombo(useAvgBjetCombo);
    if(mem_tt_HandZ.isMarkovChainIntegrator())
       mem_tt_HandZ.setMarkovChainParams(mxMode, nofBatches, nofChains,
                                         maxCallsStartingPos, epsilon0, T0, nu);
    if(higgsWidth > 0.)
      mem_tt_HandZ.setHiggsWidth(higgsWidth);
    lr_computation = LikelihoodRatio_3l1tau(
      { LikelihoodRatio_3l1tau::Hypothesis::tth }, // signal hypotheses
      { LikelihoodRatio_3l1tau::Hypothesis::ttz }  // background hypotheses
    );
  } catch(const tthMEMexception & exception)
  {
    err = exception.getErrCode();
  } catch(...)
  {
    err = TTHEXCEPTION_ERR_CODE_DEFAULT;
  }
}

MEMOutput_3l1tau
MEMInterface_3l1tau::operator()(const std::vector<MeasuredJet> & selectedJets,
                                const MeasuredLepton           & leadingLepton,
                                const MeasuredLepton           & subLeadingLepton,
                                const MeasuredLepton           & thirdLepton,
                                const MeasuredHadronicTau      & selectedHadronicTau,
                                const MeasuredMET              & measuredMET)
{
  MeasuredEvent_3l1tau event;

  if(! err)
  {
    try
    {
      if(! (selectedJets.size() >= MIN_NOF_RECO_JETS &&
            selectedJets.size() <= MAX_NOF_RECO_JETS))
        throw_line_ext("MEMInterface_3l1tau", TTHEXCEPTION_ERR_CODE_INVALID_NOF_JETS)
          << "Invalid number of selected jets passed: " << selectedJets.size()
          << " (should be between " << MIN_NOF_RECO_JETS << " and " << MAX_NOF_RECO_JETS << ')';

      event.leptons[0] = leadingLepton;
      event.leptons[1] = subLeadingLepton;
      event.leptons[2] = thirdLepton;
      for(unsigned i = 0; i < selectedJets.size(); ++i)
        event.allJets[i] = selectedJets[i];
      event.njets = selectedJets.size();
      event.htau = selectedHadronicTau;
      event.met  = measuredMET;
      event.rle.reset();
      event.initialize();
    } catch(const tthMEMexception & exception)
    {
      err = exception.getErrCode();
    } catch(...)
    {
      err = TTHEXCEPTION_ERR_CODE_DEFAULT;
    }
  }

  const std::array<double, 2> probSignalResult = [&event, this]()
    -> std::array<double, 2>
  {
    if(! err)
    {
      try
      {
        return mem_tt_HandZ.integrate(event, ME_mg5_3l1tau::kTTH);
      } catch(const tthMEMexception & exception)
      {
        err = exception.getErrCode();
      } catch(...)
      {
        err = TTHEXCEPTION_ERR_CODE_DEFAULT;
      }
    }
    return {{0., 0.}}; // never used (see below)
  }();
  const std::array<double, 2> probBackgroundResult_ttz = [&event, this]()
    -> std::array<double, 2>
  {
    if(! err)
    {
      try
      {
        return mem_tt_HandZ.integrate(event, ME_mg5_3l1tau::kTTZ);
      } catch(const tthMEMexception & exception)
      {
        err = exception.getErrCode();
      } catch(...)
      {
        err = TTHEXCEPTION_ERR_CODE_DEFAULT;
      }
    }
    return {{0., 0.}}; // never used (see below)
  }();

  const MEMOutput_3l1tau result =
    [&probSignalResult, &probBackgroundResult_ttz, this]()
  {
    if(! err)
    {
      try
      {
        return lr_computation.compute(
          { { LikelihoodRatio_3l1tau::Hypothesis::tth, probSignalResult         } },
          { { LikelihoodRatio_3l1tau::Hypothesis::ttz, probBackgroundResult_ttz } }
        );
      } catch(const tthMEMexception & exception)
      {
        err = exception.getErrCode();
      } catch(...)
      {
        err = TTHEXCEPTION_ERR_CODE_DEFAULT;
      }
   }
   return MEMOutput_3l1tau(err);
  }();

  return result;
}
