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
  CPPUNIT_TEST(testHmass);
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

  /* Test event: 1:13167:2616317 from the sample
   * /ttHToNonbb_M125_13TeV_powheg_pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM
   */
  const Vector beamAxis{0., 0., 1.};
  const LorentzVector hTauLepton = getLorentzVector(LV{107.000, -0.828, -1.700,   0.739}); // lepton          from tau decaying hadronically (aka hadronic tau)
  const LorentzVector hTauNu     = getLorentzVector(LV{ 39.100, -0.813, -1.690,   0.000}); // neutrino        from tau decaying hadronically
  const LorentzVector hTau       = getLorentzVector(LV{146.000, -0.824, -1.700,   1.780}); // tau decaying hadronically
  const LorentzVector lTauLepton = getLorentzVector(LV{ 40.000, -0.092, -2.290,   0.106}); // lepton             from tau decaying leptonically
  const LorentzVector lTauLeptNu = getLorentzVector(LV{ 45.900, -0.097, -2.280,   0.000}); // lepton neutrino    from tau decaying leptonically
  const LorentzVector lTauTauNu  = getLorentzVector(LV{ 36.700, -0.116, -2.260,   0.000}); // tau neutrino       from tau decaying leptonically
  const LorentzVector lTauDiNu   = getLorentzVector(LV{ 82.600, -0.106, -2.270,   1.130}); // di-neutrino system from tau decaying leptonically
  const LorentzVector lTau       = getLorentzVector(LV{123.000, -0.101, -2.280,   1.780}); // tau decaying leptonically
  const LorentzVector H          = getLorentzVector(LV{258.000, -0.542, -1.960, 125.000}); // Higgs
  const LorentzVector lW0        = getLorentzVector(LV{100.000, -0.416,  0.971,   0.106}); // lepton   from W boson decay (1)
  const LorentzVector nuW0       = getLorentzVector(LV{ 18.200, -1.590, -0.605,   0.000}); // neutrino from W boson decay (1)
  const LorentzVector b0         = getLorentzVector(LV{207.000, -0.380,  1.500,   4.750}); // b quark  from t quark decay (1)
  const LorentzVector W0         = getLorentzVector(LV{102.000, -0.766,  0.791,  80.400}); // W boson  from t quark decay (1)
  const LorentzVector t0         = getLorentzVector(LV{292.000, -0.543,  1.270, 173.000}); // t quark                     (1)
  const LorentzVector lW1        = getLorentzVector(LV{ 13.900, -0.198,  2.240,   0.106}); // lepton   from W boson decay (2)
  const LorentzVector nuW1       = getLorentzVector(LV{125.000, -0.494, -1.560,   0.000}); // neutrino from W boson decay (2)
  const LorentzVector b1         = getLorentzVector(LV{ 43.900, -1.250,  0.753,   4.750}); // b quark  from t quark decay (2)
  const LorentzVector W1         = getLorentzVector(LV{114.000, -0.557, -1.630,  79.900}); // W boson  from t quark decay (2)
  const LorentzVector t1         = getLorentzVector(LV{ 87.500, -1.230, -1.280, 173.000}); // t quark                     (2)
  const double measuredVisMassSquared = (hTauLepton + lTauLepton).mass2();

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
   * @brief Test whether reconstructed H mass is good enough for the unit tests
   */
  void
  testHmass()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(H.mass(), constants::massHiggs, constants::gammaHiggs);
  }

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
