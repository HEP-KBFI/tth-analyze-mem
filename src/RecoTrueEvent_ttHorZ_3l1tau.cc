#include "tthAnalysis/tthMEM/interface/RecoTrueEvent_ttHorZ_3l1tau.h"

#include <Math/VectorUtil.h> // ROOT::Math::VectorUtil::boost()

using namespace tthMEM;
namespace VectorUtil = ROOT::Math::VectorUtil;

RecoTrueEvent_ttHorZ_3l1tau
RecoTrueEvent_ttHorZ_3l1tau::boost(Vector boostVector) const
{
  RecoTrueEvent_ttHorZ_3l1tau boosted;

  boosted.hTauLepton = VectorUtil::boost(hTauLepton, boostVector);
  boosted.lTauLepton = VectorUtil::boost(lTauLepton, boostVector);
  boosted.nuHtau     = VectorUtil::boost(nuHtau,     boostVector);
  boosted.nuLtau     = VectorUtil::boost(nuLtau,     boostVector);
  boosted.hTau       = VectorUtil::boost(hTau,       boostVector);
  boosted.lTau       = VectorUtil::boost(lTau,       boostVector);
  boosted.higgsOrZ   = VectorUtil::boost(higgsOrZ,   boostVector);

  for(unsigned i = 0; i < 2; ++i)
  {
    boosted.lW[i]  = VectorUtil::boost(lW[i],  boostVector);
    boosted.nuW[i] = VectorUtil::boost(nuW[i], boostVector);
    boosted.W[i]   = VectorUtil::boost(W[i],   boostVector);
    boosted.b[i]   = VectorUtil::boost(b[i],   boostVector);
    boosted.t[i]   = VectorUtil::boost(t[i],   boostVector);

    boosted.g[i] = g[i];
  }

  return boosted;
}

LorentzVector
RecoTrueEvent_ttHorZ_3l1tau::getNuSum() const
{
  return nuHtau + nuLtau + nuW[0] + nuW[1];
}

LorentzVector
RecoTrueEvent_ttHorZ_3l1tau::getTTHorZ() const
{
  return higgsOrZ + t[0] + t[1];
}

std::vector<LorentzVector>
RecoTrueEvent_ttHorZ_3l1tau::getForMG() const
{
  return { g[0], g[1], b[0], lW[0], nuW[0], b[1], lW[1], nuW[1], lTau, hTau };
}
