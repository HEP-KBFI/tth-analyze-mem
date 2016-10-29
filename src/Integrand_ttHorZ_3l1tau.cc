#include "tthAnalysis/tthMEM/interface/Integrand_ttHorZ_3l1tau.h"
#include "tthAnalysis/tthMEM/interface/me_tth_3l1tau_mg5.h"
#include "tthAnalysis/tthMEM/interface/me_ttz_3l1tau_mg5.h"
#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h"
#include "tthAnalysis/tthMEM/interface/tthMEMrecFunctions.h"
#include "tthAnalysis/tthMEM/interface/tthMEMconstants.h"
#include "tthAnalysis/tthMEM/interface/BJetTransferFunction.h"
#include "tthAnalysis/tthMEM/interface/Logger.h"
#include "tthAnalysis/tthMEM/interface/Exception.h" // throw_line()

#include <cmath> // std::sqrt()
#include <cstring> // std::memset()
#include <algorithm> // std::for_each(), std::copy()
#include <sstream> // std::ostringstream
#include <iomanip> // std::setprecision()

#include <TMath.h> // TMath::IsNaN() ...
 // ... (why not use std here: http://stackoverflow.com/a/570694)

using namespace tthMEM;

const Integrand_ttHorZ_3l1tau * Integrand_ttHorZ_3l1tau::gIntegrand = 0;

Integrand_ttHorZ_3l1tau::Integrand_ttHorZ_3l1tau(const std::string & pdfName,
                                                 const std::string & madgraphFilename,
                                                 const VariableManager_3l1tau & vm)
  : beamAxis_(0., 0., 1.)
  , pdf_(0)
  , currentME_(ME_mg5_3l1tau::kTTH) // default to tth
  , me_madgraph_{0, 0}
  , measuredEvent_(0)
  , vm_(vm)
  , bJetTF_(functions::deltaFunction)
{
  LOGTRC;

  if(! pdf_ && pdfName != "") pdf_ = LHAPDF::mkPDF(pdfName.c_str(), 0);
  else
  {
    throw_line("invalid argument") << "PDF file name empty!";
  }

  me_madgraph_[ME_mg5_3l1tau::kTTH] = new me_tth_3l1tau_mg5();
  me_madgraph_[ME_mg5_3l1tau::kTTZ] = new me_ttz_3l1tau_mg5();

  if(madgraphFilename != "")
  {
    for(unsigned i = 0; i < 2; ++i)
      me_madgraph_[i] -> initProc(madgraphFilename);
  }
  else
  {
    throw_line("invalid argument") << "Madgraph file name empty!";
  }

  for(unsigned i = 0; i < 10; ++i)
    mgMomenta_.push_back(new double[4]);
  std::for_each(mgMomenta_.begin(), mgMomenta_.end(),
    [](double * & d) { std::memset(d, 0., 4 * sizeof(double)); }
  );

  gIntegrand = this;
}

Integrand_ttHorZ_3l1tau::~Integrand_ttHorZ_3l1tau()
{
  LOGTRC;

  if(pdf_) delete pdf_;

  for(unsigned i = 0; i < 2; ++i)
    if(me_madgraph_[i]) delete me_madgraph_[i];

  measuredEvent_ = 0; // no allocation, just the address

  std::for_each(mgMomenta_.begin(), mgMomenta_.end(),
    [](double * & d) { delete d; d = 0; }
  );
}

Integrand_ttHorZ_3l1tau &
Integrand_ttHorZ_3l1tau::setCurrentME(ME_mg5_3l1tau currentME)
{
  LOGTRC;
  currentME_ = currentME;
  Q_ = currentME_ == ME_mg5_3l1tau::kTTH ?
                     constants::resolutionScaleTTH :
                     constants::resolutionScaleTTZ;
  return *this;
}

Integrand_ttHorZ_3l1tau &
Integrand_ttHorZ_3l1tau::setBJetTransferFunction(bool setTF)
{
  LOGTRC;
  bJetTF_ = setTF ? functions::bJetTF : functions::deltaFunction;
  return *this;
}

void
Integrand_ttHorZ_3l1tau::setHiggsWidth(double higgsWidth) const
{
  me_madgraph_[ME_mg5_3l1tau::kTTH] -> setHiggsWidth(higgsWidth);
  LOGVRB_S << "Set Higgs width to = " << higgsWidth << " GeV";
}

Integrand_ttHorZ_3l1tau &
Integrand_ttHorZ_3l1tau::setEvent(const MeasuredEvent_3l1tau & measuredEvent)
{
  LOGTRC;
  measuredEvent_ = &measuredEvent;
  mgMomentaIdxs_.clear();

//--- set the variables related to the hadronic tau
  const MeasuredHadronicTau & htau = measuredEvent_ -> htau;
  recoEvent.hTauLepton = htau.p4();
  LOGVRB << lextvrap("htau p4", recoEvent.hTauLepton);
  const TMatrixD nuHtauLocalSystem = functions::nuLocalSystem(
    beamAxis_, htau.p3().unit(), "htau"
  );

//--- bind the functional arguments, so that no explicit storage is needed
//--- for the variables which remain constant during the integration
  const double hTauEnergy = htau.energy();
  const double hTauMassSquared = pow2(htau.mass());
  const double hTauMomentum = htau.p();
  const double hTauInvBeta = 1. / htau.p4().Beta();
  nuHtauCosTheta_ = [hTauEnergy, hTauMassSquared, hTauMomentum]
                    (double nuHtau_en) -> double
  {
    return functions::nuHtauCosTheta(
      nuHtau_en, hTauEnergy, hTauMassSquared, hTauMomentum
    );
  };
  hadTauPSJacobiFactor_ = [hTauMassSquared, hTauInvBeta]
                          (double z1) -> double
  {
    return functions::hadTauPSJacobiFactor(z1, hTauMassSquared, hTauInvBeta);
  };
  nuHtauEnergy_ = [hTauEnergy] (double z1) -> double
  {
    return functions::nuTauEnergy(z1, hTauEnergy);
  };
  nuHTau_ = [nuHtauLocalSystem]
            (double nuHtau_theta, double nuHtau_phi, double nuHtau_en) -> LorentzVector
  {
    return functions::nuP4( //      <- notice ->
      nuHtau_theta, nuHtau_phi, nuHtau_en, nuHtau_en, nuHtauLocalSystem
    );
  };

//--- set the variables related to the MET/hadronic recoil TF
  const MeasuredMET & met = measuredEvent_ -> met;
  const double MET_x = met.px();
  const double MET_y = met.py();
  TMatrixDSym invCovMET = met.covMET();
  const double covDet = invCovMET.Determinant();
  if(covDet != 0.)
  {
    invCovMET.Invert();
    const double MET_TF_denom = 1. / (2. * pi() * std::sqrt(invCovMET.Determinant()));
    LOGVRB << "MET TF denominator = " << MET_TF_denom;
    MET_TF_ = [MET_x, MET_y, MET_TF_denom, invCovMET]
              (double MET_x_, double MET_y_) -> double
    {
      return functions::MET_TF(
        MET_x_, MET_y_, MET_x, MET_y, MET_TF_denom, invCovMET
      );
    };
  }
  else
  {
    MET_TF_ = 0;
    LOGERR << "Cannot invert MET covariance matrix b/c det = 0";
  }

  return *this;
}

void
Integrand_ttHorZ_3l1tau::renewInputs()
{
  LOGTRC;
//--- set the variables related to the leptonic tau
  const unsigned & complLeptIdx = measuredEvent_ -> complLeptonIdx;
  const MeasuredLepton & complLepton = measuredEvent_-> leptons[complLeptIdx];

  recoEvent.lTauLepton = complLepton.p4();
  LOGVRB << lextvrap("lept p4", recoEvent.lTauLepton);
  const int complLeptCharge = complLepton.charge();
  const TMatrixD nuLtauLocalSystem = functions::nuLocalSystem(
    beamAxis_, complLepton.p3().unit(), "ltau"
  );

//--- find the measured mass of both visible tau decay products
  const double measuredVisMassSquared =
    (recoEvent.hTauLepton + recoEvent.lTauLepton).mass2();
  const double massHiggsOrZsquared = currentME_ == ME_mg5_3l1tau::kTTH ?
    constants::massHiggsSquared : constants::massZSquared;
  LOGVRB << "measured mass of visible tau decay products: "
         << std::sqrt(measuredVisMassSquared);

//--- find the kinematic variables of the rest of the leptons (from t decay)
//--- and bind the functional arguments, so that no explicit storage is needed
//--- for the variables which remain constant during the integration
  const double complLeptMassSquared = pow2(complLepton.mass());
  const double complLeptEnergy = complLepton.energy();
  const double complLeptP = complLepton.p();
  nuLeptTauCosTheta_ = [complLeptEnergy, complLeptMassSquared, complLeptP]
                       (double nuLTau_en, double mInvSquared, double nuLTau_p) -> double
  {
    return functions::nuLeptTauCosTheta(
      nuLTau_en, mInvSquared, nuLTau_p, complLeptEnergy, complLeptMassSquared, complLeptP
    );
  };
  leptTauPSJacobiFactor_ = [complLeptMassSquared, complLeptP]
                           (double mInvSquared, double z2) -> double
  {
    return functions::leptTauPSJacobiFactor(
      mInvSquared, z2, complLeptMassSquared, complLeptP
    );
  };
  nuLTauEnergy_ = [complLeptEnergy] (double z2) -> double
  {
    return functions::nuTauEnergy(z2, complLeptEnergy);
  };
  nuLTau_ = [nuLtauLocalSystem]
            (double nuLTau_theta, double nuLTau_phi, double nuLTau_en, double nuLTau_p)
              -> LorentzVector
  {
    return functions::nuP4(
      nuLTau_theta, nuLTau_phi, nuLTau_en, nuLTau_p, nuLtauLocalSystem
    );
  };
  z2_ = [measuredVisMassSquared, massHiggsOrZsquared] (double z1) -> double
  {
    return functions::z2(z1, measuredVisMassSquared, massHiggsOrZsquared);
  };
  int lept1Charge = 0;
  for(unsigned i = 0; i < 2; ++i)
  {
    const MeasuredJet & bJet_i = measuredEvent_ -> jets[i];
    const double bJetEnergy_i = bJet_i.energy();
    const double bJetEta_i = bJet_i.eta();
    const Vector bJetP3unit_i = bJet_i.p3().unit();

    const unsigned & leptIdx_i = measuredEvent_ -> bjetLeptonIdxs[i];
    const MeasuredLepton & lept_i = measuredEvent_ -> leptons[leptIdx_i];
//--- set the lepton energy equal to its momentum, thus setting it massless
    const double leptEnergy_i = lept_i.p();
    const Vector leptP3_i = lept_i.p3();
    const Vector leptP3unit_i = leptP3_i.unit();
    recoEvent.lW[i] = getLorentzVector(leptP3_i, leptEnergy_i);
    LOGVRB << lextvrap("W lept " + std::to_string(i + 1) + " p4", recoEvent.lW[i]);
    if(i == 0) lept1Charge = lept_i.charge();

    bQuarkEnergy_[i] = [bJetP3unit_i, bJetEnergy_i]
                       (const LorentzVector & W_i) -> double
    {
      return functions::bQuarkEnergy(W_i, bJetP3unit_i, bJetEnergy_i);
    };
    bJetTFBound_[i] = [bJetEnergy_i, bJetEta_i](double bEnergy_i) -> double
    {
      return functions::bJetTF(bEnergy_i, bJetEnergy_i, bJetEta_i);
    };
    tDecayJacobiFactor_[i] = [leptEnergy_i, bJetP3unit_i]
                             (const LorentzVector & W_i, double bQuarkEnergy_i,
                              double bQuarkP_i, double nuWEnergy_i) -> double
    {
      return functions::tDecayJacobiFactor(
        W_i, bQuarkEnergy_i, bQuarkP_i, nuWEnergy_i, leptEnergy_i, bJetP3unit_i
      );
    };
    nuWEnergy_[i] = [leptP3unit_i, leptEnergy_i]
                      (const VectorSpherical & nuTopPunit_i) -> double
    {
      return functions::nuWEnergy(nuTopPunit_i, leptP3unit_i, leptEnergy_i);
    };
  }

//--- MadGraph momenta legend:
// 0, g
// 1, g
// 2, b, associated W+ => associated lepton +
// 3, e+, associated W+ => associated lepton +
// 4, ve, associated W+ => associated lepton +
// 5, b~, associated W- => associated lepton -
// 6, e-, associated W- => associated lepton -
// 7, ve~, associated W- => associated lepton -
// 8, ta+ => associated lepton +
// 9, ta- => associated lepton -
  mgMomentaIdxs_.clear();
  if(lept1Charge == +1 && complLeptCharge == +1)
    mgMomentaIdxs_ = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  else
  if(lept1Charge == +1 && complLeptCharge == -1)
    mgMomentaIdxs_ = { 0, 1, 2, 3, 4, 5, 6, 7, 9, 8 };
  else
    throw_line("integrand") << "This permutation shouldn't happen";

//--- for debugging purposes plot some variables
  if(DebugPlotter_ttHorZ_3l1tau * dPlotter = measuredEvent_ -> debugPlotter)
  {
    dPlotter -> write();
    std::string measuredEventStr = measuredEvent_ -> str();
    measuredEventStr += std::string("_") +
      (currentME_ == ME_mg5_3l1tau::kTTH ? "tth" : "ttz");
    dPlotter -> initialize(measuredEventStr, vm_);
  }
}

void
Integrand_ttHorZ_3l1tau::setMGmomenta(const std::vector<LorentzVector> & memVector_p4) const
{
  for(unsigned i = 0; i < mgMomentaIdxs_.size(); ++i)
    setMGmomentum(memVector_p4[i], mgMomenta_[mgMomentaIdxs_[i]]);
}

double
Integrand_ttHorZ_3l1tau::eval(const double * x) const
{
  LOGTRC;
  if(! pdf_)                      LOGERR << "PDF not initialized!";
  if(! me_madgraph_[currentME_])  LOGERR << "Madgraph's ME not initialized!";
  if(! measuredEvent_)            LOGERR << "Measured event not specified!";
  if(mgMomentaIdxs_.size() != 10) LOGERR << "Number of MG momenta indexes not equal to 10";
  if(! pdf_ || ! me_madgraph_[currentME_] || ! measuredEvent_ || mgMomentaIdxs_.size() != 10)
    throw_line("invalid argument") << "Insufficient data for eval()";

  LOGVRB << "Current MG5 ME: " << me_madgraph_[currentME_] -> name();
  if(! MET_TF_) return 0.;
  const bool isTTH = currentME_ == ME_mg5_3l1tau::kTTH;

//--- read the sampled values
  LOGVRB << "x = { " << vm_.getArrayString(x) << " }";

  const double cosTheta[2] = { vm_.get(Var_3l1tau::kBcosTheta1, x),
                               vm_.get(Var_3l1tau::kBcosTheta2, x) };
  const double varphi[2]   = { vm_.get(Var_3l1tau::kBphi1, x),
                               vm_.get(Var_3l1tau::kBphi2, x) };
  const double z1          = vm_.get(Var_3l1tau::kZ1, x);
  const double nuHtau_phi  = vm_.get(Var_3l1tau::kTauPhi, x);
  const double nuLTau_phi  = vm_.get(Var_3l1tau::kTauPhiInv, x);
  const double mInvSquared = vm_.get(Var_3l1tau::kTauMinvSquared, x);

//--- confirm that the energy fraction carried by the tau is indeed in (0,1)
  const double z2 = roundToNdigits(z2_(z1));
  if(! (z2 >= 1.e-5 && z2 <= 1.))
  {
    LOGVRB << "z2 = " << z2 << " not in (0, 1) => p = 0";
    return 0.;
  }
  else LOGTRC << "z2 = " << z2;

//--- compute the neutrino and tau lepton 4-vector from hadronic tau
  const double nuHtau_en = nuHtauEnergy_(z1);
  const double nuHtau_cosTheta = roundToNdigits(nuHtauCosTheta_(nuHtau_en));
  if(! (nuHtau_cosTheta >= -1. && nuHtau_cosTheta <= +1.))
  {
    LOGVRB << "nuHtau_cosTheta = " << nuHtau_cosTheta << " not in (-1, 1) => p = 0";
    return 0.;
  }
  else LOGTRC << "nuHtau_cosTheta = " << nuHtau_cosTheta;
  recoEvent.nuHtau = nuHTau_(std::acos(nuHtau_cosTheta), nuHtau_phi, nuHtau_en);
  recoEvent.hTau = recoEvent.hTauLepton + recoEvent.nuHtau;
  LOGTRC << lextvrap("htau nu",      recoEvent.nuHtau);
  if(measuredEvent_ -> generatorLevel)
  {
    LOGALL << lextvrap("gen htau nu", measuredEvent_ -> generatorLevel -> genNuFromHtau.p4());
    LOGALL << std::string(20, '-');
    LOGALL << lextvrap("htau lept", recoEvent.hTauLepton);
    LOGALL << lextvrap("gen htau lept", measuredEvent_ -> generatorLevel -> genHtau.p4());
    LOGALL << std::string(20, '-');
  }
  LOGTRC << lextvrap("hadronic tau", recoEvent.hTau);
  if(measuredEvent_ -> generatorLevel)
  {
    const std::size_t & hadronicTauDecayIdx = measuredEvent_ -> generatorLevel -> hadronicTauDecayIdx;
    LOGALL << lextvrap("gen hadr tau", measuredEvent_ -> generatorLevel -> genTau[hadronicTauDecayIdx].p4());
    LOGALL << std::string(20, '-');
  }

//--- compute the neutrino and tau lepton 4-vector from leptonic tau
  const double nuLTau_en = nuLTauEnergy_(z2);
  const double nuLTau_p = std::sqrt(std::max(0., pow2(nuLTau_en) - mInvSquared));
  const double nuLTau_cosTheta = roundToNdigits(nuLeptTauCosTheta_(nuLTau_en, mInvSquared, nuLTau_p));
  if(! (nuLTau_cosTheta >= -1. && nuLTau_cosTheta <= +1.))
  {
    LOGVRB << "nuLTau_en = "  << nuLTau_en   << ", "
           << "nuLTau_p = "   << nuLTau_p    << ", "
           << "nuLTau_phi = " << nuLTau_phi  << " and "
           << "nuLTau_m = "   << std::sqrt(mInvSquared) << "; but";
    LOGVRB << "nuLTau_cosTheta = " << nuLTau_cosTheta << " not in (-1, 1) => p = 0";
    return 0.;
  }
  recoEvent.nuLtau = nuLTau_(std::acos(nuLTau_cosTheta), nuLTau_phi, nuLTau_en, nuLTau_p);
  recoEvent.lTau = recoEvent.lTauLepton + recoEvent.nuLtau;
  LOGTRC << lextvrap("lept tau di nu",  recoEvent.nuLtau);
  if(measuredEvent_ -> generatorLevel)
  {
    LOGALL << lextvrap("gen di nu", measuredEvent_ -> generatorLevel -> genDiNuFromLtau.p4());
    LOGALL << std::string(20, '-');
    LOGALL << lextvrap("lept tau lept", recoEvent.lTauLepton);
    LOGALL << lextvrap("gen ltau lept", measuredEvent_ -> generatorLevel -> genLepFromTau.p4());
    LOGALL << std::string(20, '-');
  }
  LOGTRC << lextvrap("leptonic tau", recoEvent.lTau);
  if(measuredEvent_ -> generatorLevel)
  {
    const std::size_t & leptonicTauDecayIdx = measuredEvent_ -> generatorLevel -> leptonicTauDecayIdx;
    LOGALL << lextvrap("gen lept tau", measuredEvent_ -> generatorLevel -> genTau[leptonicTauDecayIdx].p4());
    LOGALL << std::string(20, '-');
  }

  recoEvent.higgsOrZ = recoEvent.hTau + recoEvent.lTau;
  LOGTRC << lextvrap(isTTH ? "higgs" : "Z", recoEvent.higgsOrZ);
  if(measuredEvent_ -> generatorLevel)
  {
    LOGALL << lextvrap(Form("gen %s", (isTTH ? "higgs" : "Z")),
                       measuredEvent_ -> generatorLevel -> genHorZ.p4());
    LOGALL << std::string(20, '-');
  }

//--- get 4-momenta of top decay products (and top itself)
  double bEnergyTF[2];  // value of b-quark energy transfer function
  for(unsigned i = 0; i < 2; ++i)
  {
    const std::string i_str = std::to_string(i + 1);
    const VectorSpherical nuWP3unit_i(1., std::acos(cosTheta[i]), varphi[i]);
    const double nuWenergy_i = nuWEnergy_[i](nuWP3unit_i);
    recoEvent.nuW[i] = getLorentzVector(Vector(nuWenergy_i * nuWP3unit_i), nuWenergy_i);
    recoEvent.W[i] = recoEvent.nuW[i] + recoEvent.lW[i];

    LOGTRC << lextvrap("W nu " + i_str, recoEvent.nuW[i]);
    if(measuredEvent_ -> generatorLevel)
    {
      LOGALL << lextvrap("gen W nu " + i_str,
                         measuredEvent_ -> generatorLevel -> genNuFromTop[i].p4());
      LOGALL << std::string(20, '-');
    }
    LOGTRC << lextvrap("W " + i_str, recoEvent.W[i]);
    if(measuredEvent_ -> generatorLevel)
    {
      LOGALL << lextvrap("gen W " + i_str,
                         measuredEvent_ -> generatorLevel -> genWBoson[i].p4());
      LOGALL << std::string(20, '-');
    }

    const double bEnergy_i = roundToNdigits(bQuarkEnergy_[i](recoEvent.W[i]));
    if(bEnergy_i == 0.) return 0.;
    const MeasuredJet & bJet_i = measuredEvent_ -> jets[i];
    const Vector bP3_i = std::sqrt(pow2(bEnergy_i) - constants::massBSquared) *
                         bJet_i.p3().unit();
    recoEvent.b[i] = getLorentzVector(bP3_i, bEnergy_i);
    recoEvent.t[i] = recoEvent.b[i] + recoEvent.W[i];

    LOGTRC << lextvrap("b " + i_str,           recoEvent.b[i]);
    LOGTRC << lextvrap("b " + i_str + " reco", bJet_i.p4());
    if(measuredEvent_ -> generatorLevel) LOGALL << std::string(20, '-');
    LOGTRC << lextvrap("t " + i_str,           recoEvent.t[i]);
    if(measuredEvent_ -> generatorLevel)
    {
      LOGALL << lextvrap("gen t " + i_str,
                         measuredEvent_ -> generatorLevel -> genTop[i].p4());
      LOGALL << std::string(20, '-');
    }

//--- b-jet energy transfer function
    bEnergyTF[i] = bJetTFBound_[i](bEnergy_i);
  }

  for(unsigned i = 0; i < 2; ++i)
    LOGTRC_S << "TF for #" << i << " b-quark energy = " << bEnergyTF[i];

//--- hadronic recoil transfer function; simplification: use only neutrinos
//--- in the difference of ,,measured'' and ,,true'' hadronic recoil, because
//--- leptons are assumed to be measured perfectly whereas the pT of jets/quarks
//--- introduces uncertainties difficult to handle here
  const LorentzVector nuSum = recoEvent.getNuSum();
  const LorentzVector bJetDifference = [this]() -> LorentzVector
  {
    LorentzVector result(0., 0., 0., 0.);
    for(unsigned i = 0; i < 2; ++i)
      result += recoEvent.b[i] - (measuredEvent_ -> jets[i]).p4();
    return result;
  }();
  const double MET_TF = MET_TF_(nuSum.x() + bJetDifference.x(),
                                nuSum.y() + bJetDifference.y());

//--- compute Bjorken x variables
//--- assume that hadronic recoil has only transverse component
  const double hadRecE = 0.;
  const double hadRecPz = 0.;
  const LorentzVector tthOrZ = recoEvent.getTTHorZ();
  const double xa = roundToNdigits((hadRecE + tthOrZ.e() + hadRecPz + tthOrZ.pz()) *
                    constants::invSqrtS);
  const double xb = roundToNdigits((hadRecE + tthOrZ.e() - hadRecPz - tthOrZ.pz()) *
                    constants::invSqrtS);
  LOGTRC << "xa = " << xa << "; xb = " << xb;
  if(xa <= 0. || xa >= 1. || xb <= 0. || xb >= 1.)
  {
    LOGVRB << "xa or xb have unphysical values";
    return 0.;
  }

  const double fa = pdf_ -> xfxQ(21, xa, Q_) / xa;
  const double fb = pdf_ -> xfxQ(21, xb, Q_) / xb;
  const double probPDF = fa * fb;
  const double flux = constants::invS / (xa * xb);

//--- boost all MG momenta into frame where pT(tth) = 0
//--- note: the boost vector must be a 3-vector of velocities, hence
//---       the division by energy
  const Vector boost(-tthOrZ.px() / tthOrZ.e(), -tthOrZ.py() / tthOrZ.e(), 0.);
  LOGTRC << cvrap("boost vector", boost);
  recoEvent.g[0] = LorentzVector(0., 0., +0.5 * xa * constants::sqrtS,
                                          0.5 * xa * constants::sqrtS);
  recoEvent.g[1] = LorentzVector(0., 0., -0.5 * xb * constants::sqrtS,
                                          0.5 * xb * constants::sqrtS);
  const RecoTrueEvent_ttHorZ_3l1tau recoEventmem = recoEvent.boost(boost);
  const LorentzVector tthOrZ_mem = recoEventmem.getTTHorZ();
  if(isTTH) LOGTRC << lmvrap("tth mem", tthOrZ_mem);
  else      LOGTRC << lmvrap("ttz mem", tthOrZ_mem);

  const double z1_mem = roundToNdigits(recoEventmem.hTauLepton.e() / recoEventmem.hTau.e());
  const double z2_mem = roundToNdigits(recoEventmem.lTauLepton.e() / recoEventmem.lTau.e());
  LOGTRC << "z1_mem = " << z1_mem << "; z2_mem = " << z2_mem;
  if(! (z1_mem >= 1.e-5 && z1_mem <= 1.) ||
     ! (z2_mem >= 1.e-5 && z2_mem <= 1.))
  {
    LOGVRB << "The tau energy fraction z1 and z2 have unphysical values "
           << "when the lab frame is boosted such that pT(tth) = 0";
    return 0.;
  }

//--- set MG momenta & calculate the matrix element
  setMGmomenta(recoEventmem.getForMG());
  me_madgraph_[currentME_] -> setMomenta(mgMomenta_);
  me_madgraph_[currentME_] -> sigmaKin();
  const double prob_ME_mg = me_madgraph_[currentME_] -> getMatrixElements()[0];
  if(TMath::IsNaN(prob_ME_mg) || prob_ME_mg < 0.)
  {
    LOGWARN_S << "MadGraph5 returned NaN or is zero: "
              << "|M|^2 = " << prob_ME_mg << " => skipping event";
    return 0.;
  }
  LOGTRC_S << "prob(PDF) = " << probPDF << "; flux factor = " << flux;
  LOGVRB_S << "|M|^2 = " << prob_ME_mg;

//--- jacobi factors
  double tDecayJacobiFactors[2];
  for(unsigned i = 0; i < 2; ++i)
  {
    tDecayJacobiFactors[i] = tDecayJacobiFactor_[i](
      recoEvent.W[i], recoEvent.b[i].e(), recoEvent.b[i].P(), recoEvent.nuW[i].e()
    );
    LOGTRC_S << "Jacobi factors arising from the top decay: #" << (i + 1) << " = "
             << tDecayJacobiFactors[i];
    if(tDecayJacobiFactors[i] == 0.) return 0.;
  }

  const double hTauPSJacobiFactor = hadTauPSJacobiFactor_(z1);
  LOGTRC_S << "PS x Jacobi factor for hadronic tau decay = "
           << hTauPSJacobiFactor;
  const double lTauPSJacobiFactor = leptTauPSJacobiFactor_(mInvSquared, z2);
  if(lTauPSJacobiFactor == 0.) return 0.;
  LOGTRC_S << "PS x Jacobi factor for leptonic tau decay = "
           << lTauPSJacobiFactor;

  const double jacobiFactor = tDecayJacobiFactors[0] * tDecayJacobiFactors[1] *
    hTauPSJacobiFactor * lTauPSJacobiFactor * z2 *
    (isTTH ? constants::ttHfactor : constants::ttZfactor);
  LOGTRC_S << "Product of all Jacobi factors = " << jacobiFactor;

//--- assemble the integrand
  const double p = prob_ME_mg * probPDF * flux * MET_TF *
                   bEnergyTF[0] * bEnergyTF[1] * jacobiFactor;
  LOGVRB_S << "p = " << p;

//--- for debugging purposes plot some variables
  if(DebugPlotter_ttHorZ_3l1tau * dPlotter = measuredEvent_ -> debugPlotter)
    (*dPlotter).fill(hVar_3l1tau::kZ2,         z2)
               .fill(hVar_3l1tau::kMassHorZ,   recoEvent.higgsOrZ.mass())
               .fill(hVar_3l1tau::kMassHtau,   recoEvent.hTau.mass())
               .fill(hVar_3l1tau::kMassLtau,   recoEvent.lTau.mass())
               .fill(hVar_3l1tau::kB1en,       recoEvent.b[0].e())
               .fill(hVar_3l1tau::kB2en,       recoEvent.b[1].e())
               .fill(hVar_3l1tau::kB1RecoEn,   measuredEvent_ -> jets[0].energy())
               .fill(hVar_3l1tau::kB2RecoEn,   measuredEvent_ -> jets[1].energy())
               .fill(hVar_3l1tau::kMETtf,      MET_TF)
               .fill(hVar_3l1tau::kB1energyTF, bEnergyTF[0])
               .fill(hVar_3l1tau::kB2energyTF, bEnergyTF[1])
               .fill(hVar_3l1tau::kTdecayJF1,  tDecayJacobiFactors[0])
               .fill(hVar_3l1tau::kTdecayJF2,  tDecayJacobiFactors[1])
               .fill(hVar_3l1tau::kHtauPSF,    hTauPSJacobiFactor)
               .fill(hVar_3l1tau::kLtauPSF,    lTauPSJacobiFactor)
               .fill(hVar_3l1tau::kJacobiF,    jacobiFactor)
               .fill(hVar_3l1tau::kXa,         xa)
               .fill(hVar_3l1tau::kXb,         xb)
               .fill(hVar_3l1tau::kFlux,       flux)
               .fill(hVar_3l1tau::kProbPDF,    probPDF)
               .fill(hVar_3l1tau::kMsquared,   prob_ME_mg)
               .fill(hVar_3l1tau::kProb,       p)
               .fill(vm_, x)
               .fill();

  return p;
}
