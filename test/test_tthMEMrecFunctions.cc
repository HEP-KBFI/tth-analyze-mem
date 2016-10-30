#include <cppunit/extensions/HelperMacros.h> // CppUnit::TestFixture, CPPUNIT_ASSERT_*

#include <cmath> // std::acos()

#include "tthAnalysis/tthMEM/interface/tthMEMconstants.h" // tthMEM::constants::
#include "tthAnalysis/tthMEM/interface/tthMEMrecFunctions.h" // tthMEM::functions::

using namespace tthMEM;
typedef math::PtEtaPhiMLorentzVector LV;

class Test_tthMEMrecFunctions
  : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(Test_tthMEMrecFunctions);
//-------------------------------------------------------------------------------------------------
  CPPUNIT_TEST(testHtauMass);
  CPPUNIT_TEST(testLtauMass);
  CPPUNIT_TEST(testHmass);
  CPPUNIT_TEST(testW1mass);
  CPPUNIT_TEST(testW2mass);
  CPPUNIT_TEST(testTop1mass);
  CPPUNIT_TEST(testTop2mass);
//-------------------------------------------------------------------------------------------------
  CPPUNIT_TEST(testDiNeutrinoSystem_e);
  CPPUNIT_TEST(testDiNeutrinoSystem_m);
  CPPUNIT_TEST(testDiNeutrinoSystem_px);
  CPPUNIT_TEST(testDiNeutrinoSystem_py);
  CPPUNIT_TEST(testDiNeutrinoSystem_pz);
  CPPUNIT_TEST(testLeptonicTauBranch_e);
  CPPUNIT_TEST(testLeptonicTauBranch_m);
  CPPUNIT_TEST(testLeptonicTauBranch_px);
  CPPUNIT_TEST(testLeptonicTauBranch_py);
  CPPUNIT_TEST(testLeptonicTauBranch_pz);
  CPPUNIT_TEST(testHadronicTauBranch_e);
  CPPUNIT_TEST(testHadronicTauBranch_m);
  CPPUNIT_TEST(testHadronicTauBranch_px);
  CPPUNIT_TEST(testHadronicTauBranch_py);
  CPPUNIT_TEST(testHadronicTauBranch_pz);
  CPPUNIT_TEST(testHiggsReconstruction_e);
  CPPUNIT_TEST(testHiggsReconstruction_m);
  CPPUNIT_TEST(testHiggsReconstruction_px);
  CPPUNIT_TEST(testHiggsReconstruction_py);
  CPPUNIT_TEST(testHiggsReconstruction_pz);
  CPPUNIT_TEST(testW1reconstruction_e);
  CPPUNIT_TEST(testW1reconstruction_m);
  CPPUNIT_TEST(testW1reconstruction_px);
  CPPUNIT_TEST(testW1reconstruction_py);
  CPPUNIT_TEST(testW1reconstruction_pz);
  CPPUNIT_TEST(testW2reconstruction_e);
  CPPUNIT_TEST(testW2reconstruction_m);
  CPPUNIT_TEST(testW2reconstruction_px);
  CPPUNIT_TEST(testW2reconstruction_py);
  CPPUNIT_TEST(testW2reconstruction_pz);
  CPPUNIT_TEST(testTop1Reconstruction_e);
  CPPUNIT_TEST(testTop1Reconstruction_m);
  CPPUNIT_TEST(testTop1Reconstruction_px);
  CPPUNIT_TEST(testTop1Reconstruction_py);
  CPPUNIT_TEST(testTop1Reconstruction_pz);
  CPPUNIT_TEST(testTop2Reconstruction_e);
  CPPUNIT_TEST(testTop2Reconstruction_m);
  CPPUNIT_TEST(testTop2Reconstruction_px);
  CPPUNIT_TEST(testTop2Reconstruction_py);
  CPPUNIT_TEST(testTop2Reconstruction_pz);
//-------------------------------------------------------------------------------------------------
  CPPUNIT_TEST(testZ2);
  CPPUNIT_TEST(testNuHTauEnergy);
  CPPUNIT_TEST(testNuLTauEnergy);
  CPPUNIT_TEST(testNuHtauCosTheta);
  CPPUNIT_TEST(testNuLeptTauCosTheta);
  CPPUNIT_TEST(testHtauPhi_e);
  CPPUNIT_TEST(testHtauPhi_px);
  CPPUNIT_TEST(testHtauPhi_py);
  CPPUNIT_TEST(testHtauPhi_pz);
  CPPUNIT_TEST(testLtauPhi_e);
  CPPUNIT_TEST(testLtauPhi_px);
  CPPUNIT_TEST(testLtauPhi_py);
  CPPUNIT_TEST(testLtauPhi_pz);
  CPPUNIT_TEST(testNuWEnergy);
  CPPUNIT_TEST_SUITE_END();

public:

  /* Test event: 1:3347:665066 from the sample
   * /ttHToNonbb_M125_13TeV_powheg_pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM
   */
  const Vector beamAxis{0., 0., 1.};
  const LorentzVector hTauLepton = getLorentzVector(LV{ 77.240, -0.818, -2.835,   0.623}); // lepton          from tau decaying hadronically (aka hadronic tau)
  const LorentzVector hTauNu     = getLorentzVector(LV{352.600, -0.812, -2.834,   0.000}); // neutrino        from tau decaying hadronically
  const LorentzVector hTau       = getLorentzVector(LV{429.900, -0.813, -2.834,   1.777}); // tau decaying hadronically
  const LorentzVector lTauLepton = getLorentzVector(LV{ 70.060, -1.228, -2.661,   0.106}); // lepton             from tau decaying leptonically
  const LorentzVector lTauLeptNu = getLorentzVector(LV{ 77.080, -1.216, -2.655,   0.000}); // lepton neutrino    from tau decaying leptonically
  const LorentzVector lTauTauNu  = getLorentzVector(LV{ 34.030, -1.214, -2.641,   0.000}); // tau neutrino       from tau decaying leptonically
  const LorentzVector lTauDiNu   = getLorentzVector(LV{111.100, -1.215, -2.651,   0.724}); // di-neutrino system from tau decaying leptonically
  const LorentzVector lTau       = getLorentzVector(LV{181.200, -1.220, -2.655,   1.777}); // tau decaying leptonically
  const LorentzVector H          = getLorentzVector(LV{609.100, -0.950, -2.781, 124.900}); // Higgs
  const LorentzVector lW0        = getLorentzVector(LV{ 67.670, -0.107, -3.006,   0.001}); // lepton   from W boson decay (1)
  const LorentzVector nuW0       = getLorentzVector(LV{ 94.340, -0.699, -2.198,   0.000}); // neutrino from W boson decay (1)
  const LorentzVector b0         = getLorentzVector(LV{124.900, -0.708,  2.777,   4.750}); // b quark  from t quark decay (1)
  const LorentzVector W0         = getLorentzVector(LV{149.300, -0.505, -2.532,  79.070}); // W boson  from t quark decay (1)
  const LorentzVector t0         = getLorentzVector(LV{242.600, -0.669, -2.972, 170.200}); // t quark                     (1)
  const LorentzVector lW1        = getLorentzVector(LV{137.000, -0.204,  0.520,   0.106}); // lepton   from W boson decay (2)
  const LorentzVector nuW1       = getLorentzVector(LV{366.000,  0.087,  0.302,   0.000}); // neutrino from W boson decay (2)
  const LorentzVector b1         = getLorentzVector(LV{305.200, -0.249,  0.112,   4.750}); // b quark  from t quark decay (2)
  const LorentzVector W1         = getLorentzVector(LV{500.600,  0.007,  0.361,  81.530}); // W boson  from t quark decay (2)
  const LorentzVector t1         = getLorentzVector(LV{799.900, -0.092,  0.267, 174.400}); // t quark                     (2)
  const double measuredVisMassSquared = (hTauLepton + lTauLepton).mass2();
//-----------------------------------------------------------------------------------------------------------------
  const LorentzVector lTauDiNuReco = lTauLeptNu + lTauTauNu;
  const LorentzVector lTauReco     = lTauDiNuReco + lTauLepton;
  const LorentzVector hTauReco     = hTauNu + hTauLepton;
  const LorentzVector hReco        = lTauReco + hTauReco;
  const LorentzVector W0reco       = lW0 + nuW0;
  const LorentzVector W1reco       = lW1 + nuW1;
  const LorentzVector T0reco       = W0reco + b0;
  const LorentzVector T1reco       = W1reco + b1;
//-----------------------------------------------------------------------------------------------------------------
  double z1, z2;
  double cosThetaH, cosThetaL;
  double thetaH, thetaL;
  double nuHTauPhi, nuLTauPhi;
  Vector eX_htau, eY_htau, eZ_htau;
  Vector eX_ltau, eY_ltau, eZ_ltau;
  LorentzVector hTauNu_2, lTauNu_2;

  VectorSpherical nuW0unit;
  Vector lW0unit;
  double lW0energy;

  void
  setUp() override
  {
    z1 = functions::z(hTau, hTauLepton);
    z2 = functions::z(lTau, lTauLepton);
    cosThetaH = functions::cosTheta(hTauLepton, hTauNu);
    cosThetaL = functions::cosTheta(lTauLepton, lTauDiNu);
    thetaH = std::acos(cosThetaH);
    thetaL = std::acos(cosThetaL);

    const TMatrixD nuHlocalSystem = functions::nuLocalSystem(
      beamAxis, getVector(hTauLepton).unit()
    );
    const TMatrixD nuLlocalSystem = functions::nuLocalSystem(
      beamAxis, getVector(lTauLepton).unit()
    );

    nuHTauPhi = functions::phiFromLabMomenta(hTau, hTauLepton, beamAxis);
    hTauNu_2 = functions::nuP4(
      thetaH, nuHTauPhi, hTauNu.e(), hTauNu.P(), nuHlocalSystem
    );

    nuLTauPhi = functions::phiFromLabMomenta(lTau, lTauLepton, beamAxis);
    lTauNu_2 = functions::nuP4(
      thetaL, nuLTauPhi, lTauDiNu.e(), lTauDiNu.P(), nuLlocalSystem
    );

    const Vector nuW0unitCart = getVector(nuW0).unit();
    nuW0unit = VectorSpherical(1., nuW0unitCart.theta(), nuW0unitCart.phi());
    lW0unit = getVector(lW0).unit();
    lW0energy = lW0.e();
  }

  void
  tearDown() override
  {
    /* empty */
  }

  /**
   * @brief Test whether reconstructed tau mass decaying hadronically is good enough for the unit tests
   */
  void
  testHtauMass()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(hTauReco.mass(), constants::massTau, constants::gammaTau);
  }

  /**
   * @brief Test whether reconstructed tau mass decaying leptonically is good enough for the unit tests
   */
  void
  testLtauMass()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(lTauReco.mass(), constants::massTau, constants::gammaTau);
  }

  /**
   * @brief Test whether reconstructed H mass is good enough for the unit tests
   */
  void
  testHmass()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(hReco.mass(), constants::massHiggs, constants::gammaHiggs);
  }

  /**
   * @brief Test whether reconstructed W mass in the 1st leg is good enough for the unit tests
   */
  void
  testW1mass()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(W0reco.mass(), constants::massW, constants::gammaW);
  }

  /**
   * @brief Test whether reconstructed W mass in the 2nd leg is good enough for the unit tests
   */
  void
  testW2mass()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(W1reco.mass(), constants::massW, constants::gammaW);
  }

  /**
   * @brief Test whether reconstructed top mass in the 1st leg is good enough for the unit tests
   */
  void
  testTop1mass()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(T0reco.mass(), constants::massT, constants::gammaT);
  }

  /**
   * @brief Test whether reconstructed top mass in the 2nd leg is good enough for the unit tests
   */
  void
  testTop2mass()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(T1reco.mass(), constants::massT, constants::gammaT);
  }

//-------------------------------------------------------------------------------------------------

  /**
   * @brief Test whether hardcoded and reconstructed di-neutrino system have the same energy
   */
  void
  testDiNeutrinoSystem_e()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(lTauDiNu.e(), lTauDiNuReco.e(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed di-neutrino system have the same mass
   */
  void
  testDiNeutrinoSystem_m()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(lTauDiNu.mass(), lTauDiNuReco.mass(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed di-neutrino system have the same px
   */
  void
  testDiNeutrinoSystem_px()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(lTauDiNu.px(), lTauDiNuReco.px(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed di-neutrino system have the same py
   */
  void
  testDiNeutrinoSystem_py()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(lTauDiNu.py(), lTauDiNuReco.py(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed di-neutrino system have the same pz
   */
  void
  testDiNeutrinoSystem_pz()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(lTauDiNu.pz(), lTauDiNuReco.pz(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed tau decaying leptonically have the same energy
   */
  void
  testLeptonicTauBranch_e()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(lTau.e(), lTauReco.e(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed tau decaying leptonically have the same mass
   */
  void
  testLeptonicTauBranch_m()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(lTau.mass(), lTauReco.mass(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed tau decaying leptonically have the same px
   */
  void
  testLeptonicTauBranch_px()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(lTau.px(), lTauReco.px(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed tau decaying leptonically have the same py
   */
  void
  testLeptonicTauBranch_py()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(lTau.py(), lTauReco.py(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed tau decaying leptonically have the same pz
   */
  void
  testLeptonicTauBranch_pz()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(lTau.pz(), lTauReco.pz(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed tau decaying hadronically have the same energy
   */
  void
  testHadronicTauBranch_e()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(hTau.e(), hTauReco.e(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed tau decaying hadronically have the same mass
   */
  void
  testHadronicTauBranch_m()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(hTau.mass(), hTauReco.mass(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed tau decaying hadronically have the same px
   */
  void
  testHadronicTauBranch_px()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(hTau.px(), hTauReco.px(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed tau decaying hadronically have the same py
   */
  void
  testHadronicTauBranch_py()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(hTau.py(), hTauReco.py(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed tau decaying hadronically have the same pz
   */
  void
  testHadronicTauBranch_pz()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(hTau.pz(), hTauReco.pz(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed Higgs have the same energy
   */
  void
  testHiggsReconstruction_e()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(H.e(), hReco.e(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed Higgs have the same mass
   */
  void
  testHiggsReconstruction_m()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(H.mass(), hReco.mass(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed Higgs have the same px
   */
  void
  testHiggsReconstruction_px()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(H.px(), hReco.px(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed Higgs have the same py
   */
  void
  testHiggsReconstruction_py()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(H.py(), hReco.py(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed Higgs have the same pz
   */
  void
  testHiggsReconstruction_pz()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(H.pz(), hReco.pz(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed W in the 1st leg have the same energy
   */
  void
  testW1reconstruction_e()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(W0.e(), W0reco.e(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed W in the 1st leg have the same mass
   */
  void
  testW1reconstruction_m()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(W0.mass(), W0reco.mass(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed W in the 1st leg have the same px
   */
  void
  testW1reconstruction_px()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(W0.px(), W0reco.px(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed W in the 1st leg have the same py
   */
  void
  testW1reconstruction_py()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(W0.py(), W0reco.py(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed W in the 1st leg have the same pz
   */
  void
  testW1reconstruction_pz()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(W0.pz(), W0reco.pz(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed W in the 2nd leg have the same energy
   */
  void
  testW2reconstruction_e()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(W1.e(), W1reco.e(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed W in the 2nd leg have the same mass
   */
  void
  testW2reconstruction_m()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(W1.mass(), W1reco.mass(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed W in the 2nd leg have the same px
   */
  void
  testW2reconstruction_px()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(W1.px(), W1reco.px(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed W in the 2nd leg have the same py
   */
  void
  testW2reconstruction_py()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(W1.py(), W1reco.py(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed W in the 2nd leg have the same pz
   */
  void
  testW2reconstruction_pz()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(W1.pz(), W1reco.pz(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed top in the 1st leg have the same energy
   */
  void
  testTop1Reconstruction_e()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(t0.e(), T0reco.e(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed top in the 1st leg have the same mass
   */
  void
  testTop1Reconstruction_m()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(t0.mass(), T0reco.mass(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed top in the 1st leg have the same px
   */
  void
  testTop1Reconstruction_px()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(t0.px(), T0reco.px(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed top in the 1st leg have the same py
   */
  void
  testTop1Reconstruction_py()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(t0.py(), T0reco.py(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed top in the 1st leg have the same pz
   */
  void
  testTop1Reconstruction_pz()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(t0.pz(), T0reco.pz(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed top in the 2nd leg have the same energy
   */
  void
  testTop2Reconstruction_e()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(t1.e(), T1reco.e(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed top in the 2nd leg have the same mass
   */
  void
  testTop2Reconstruction_m()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(t1.mass(), T1reco.mass(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed top in the 2nd leg have the same px
   */
  void
  testTop2Reconstruction_px()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(t1.px(), T1reco.px(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed top in the 2nd leg have the same py
   */
  void
  testTop2Reconstruction_py()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(t1.py(), T1reco.py(), +1e-3);
  }

  /**
   * @brief Test whether hardcoded and reconstructed top in the 2nd leg have the same pz
   */
  void
  testTop2Reconstruction_pz()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(t1.pz(), T1reco.pz(), +1e-3);
  }

//-------------------------------------------------------------------------------------------------

  /**
   * @brief Test energy fraction for letponic branch in tau decay system
   */
  void
  testZ2()
  {
    const double z2_2 = functions::z2(z1, measuredVisMassSquared, constants::massZSquared);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(z2, z2_2, +1.e-3);
  }

  /**
   * @brief Test the energy fraction carried by invisible hadronic tau decay products
   */
  void
  testNuHTauEnergy()
  {
    const double nuHtauEnergy = functions::nuTauEnergy(z1, hTauLepton.e());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(hTauNu.e(), nuHtauEnergy, +1.e-3);
  }

  /**
   * @brief Test the energy fraction carried by invisible leptonic tau decay products
   */
  void
  testNuLTauEnergy()
  {
    const double nuLtauEnergy = functions::nuTauEnergy(z2, lTauLepton.e());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(lTauDiNu.e(), nuLtauEnergy, +1.e-3);
  }

  /**
   * @brief Test cosine of the opening angle in hadronic tau decay system
   */
  void
  testNuHtauCosTheta()
  {
    const double cosThetaH_2 = functions::nuHtauCosTheta(
      hTauNu.e(), hTauLepton.e(), hTauLepton.M2(), hTauLepton.P()
    );
    CPPUNIT_ASSERT_DOUBLES_EQUAL(cosThetaH, cosThetaH_2, +1.e-3);
  }

  /**
   * @brief Test cosinge of the opening angle in leptonic tau decay system
   */
  void
  testNuLeptTauCosTheta()
  {
    const double cosThetaL_2 = functions::nuLeptTauCosTheta(
      lTauDiNu.e(), lTauDiNu.M2(), lTauDiNu.P(), lTauLepton.e(), lTauLepton.M2(), lTauLepton.P()
    );
    CPPUNIT_ASSERT_DOUBLES_EQUAL(cosThetaL, cosThetaL_2, +1.e-3);
  }

  /**
   * @brief Test rotation angle correctness in hadronic tau decay system
   *        (neutrino energy equality)
   */
  void
  testHtauPhi_e()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(hTauNu_2.e(),  hTauNu.e(),  +1.e-3);
  }

  /**
   * @brief Test rotation angle correctness in hadronic tau decay system
   *        (neutrino px equality)
   */
  void
  testHtauPhi_px()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(hTauNu_2.px(), hTauNu.px(), +1.e-3);
  }

  /**
   * @brief Test rotation angle correctness in hadronic tau decay system
   *        (neutrino py equality)
   */
  void
  testHtauPhi_py()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(hTauNu_2.py(), hTauNu.py(), +1.e-3);
  }

  /**
   * @brief Test rotation angle correctness in hadronic tau decay system
   *        (neutrino pz equality)
   */
  void
  testHtauPhi_pz()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(hTauNu_2.pz(), hTauNu.pz(), +1.e-3);
  }

  /**
   * @brief Test rotation angle correctness in leptonic tau decay system
   *        (neutrino energy equality)
   */
  void
  testLtauPhi_e()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(lTauNu_2.e(),  lTauDiNu.e(),  +1.e-3);
  }

  /**
   * @brief Test rotation angle correctness in leptonic tau decay system
   *        (neutrino px equality)
   */
  void
  testLtauPhi_px()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(lTauNu_2.px(), lTauDiNu.px(), +1.e-3);
  }

  /**
   * @brief Test rotation angle correctness in leptonic tau decay system
   *        (neutrino py equality)
   */
  void
  testLtauPhi_py()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(lTauNu_2.py(), lTauDiNu.py(), +1.e-3);
  }

  /**
   * @brief Test rotation angle correctness in leptonic tau decay system
   *        (neutrino pz equality)
   */
  void
  testLtauPhi_pz()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(lTauNu_2.pz(), lTauDiNu.pz(), +1.e-3);
  }

  /**
   * @brief Test neutrino energy coming from top/associated W decay
   */
  void
  testNuWEnergy()
  {
    const double nuWEnergy = functions::nuWEnergy(nuW0unit, lW0unit, lW0energy);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(nuWEnergy, nuW0.e(), +1.e-3);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(Test_tthMEMrecFunctions);
