#include <cppunit/extensions/HelperMacros.h> // CppUnit::TestFixture, CPPUNIT_ASSERT_*

#include <cmath> // std::acos()

#include "tthAnalysis/tthMEM/interface/tthMEMconstants.h" // tthMEM::constants::
#include "tthAnalysis/tthMEM/interface/tthMEMrecFunctions.h" // tthMEM::functions::

using namespace tthMEM;

class Test_tthMEMrecFunctions
  : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(Test_tthMEMrecFunctions);
  CPPUNIT_TEST(testZmass);
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

  /* testing TTZ system for now */
  /* @todo find better Z system to test with (testZmass fails) */
  const Vector beamAxis{0., 0., 1.};
  const LorentzVector hTauLepton{-19.639, 14.979, -0.766, 24.716};
  const LorentzVector hTauNu    {-11.378,  9.964, -1.294, 15.179};
  const LorentzVector hTau = hTauLepton + hTauNu;
  const LorentzVector lTauLepton{-18.860,  77.334,  -39.071,  88.672};
  const LorentzVector lTauNu    {-59.853, 237.870, -117.872, 272.137};
  const LorentzVector lTau = lTauLepton + lTauNu;
  const LorentzVector Z = hTau + lTau;
  const double measuredVisMassSquared = (hTauLepton + lTauLepton).mass2();
  /* testing only one top decay leg */
  const LorentzVector lW0 {90.038,  19.851, 171.645, 194.841};
  const LorentzVector nuW0{29.901, -13.992,   5.014,  33.392};

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
    cosThetaL = functions::cosTheta(lTauLepton, lTauNu);
    thetaH = std::acos(cosThetaH);
    thetaL = std::acos(cosThetaL);

    const Vector eZ_htau_ = getVector(hTauLepton).unit();
    const Vector eY_htau_ = eZ_htau_.Cross(beamAxis).unit();
    const Vector eX_htau_ = eY_htau_.Cross(eZ_htau_).unit();
    eX_htau = Vector(eX_htau_.x(), eY_htau_.x(), eZ_htau_.x());
    eY_htau = Vector(eX_htau_.y(), eY_htau_.y(), eZ_htau_.y());
    eZ_htau = Vector(eX_htau_.z(), eY_htau_.z(), eZ_htau_.z());

    const Vector eZ_ltau_ = getVector(lTauLepton).unit();
    const Vector eY_ltau_ = eZ_ltau_.Cross(beamAxis).unit();
    const Vector eX_ltau_ = eY_ltau_.Cross(eZ_ltau_).unit();
    eX_ltau = Vector(eX_ltau_.x(), eY_ltau_.x(), eZ_ltau_.x());
    eY_ltau = Vector(eX_ltau_.y(), eY_ltau_.y(), eZ_ltau_.y());
    eZ_ltau = Vector(eX_ltau_.z(), eY_ltau_.z(), eZ_ltau_.z());

    nuHTauPhi = functions::phiFromLabMomenta(hTau, hTauLepton, beamAxis);
    hTauNu_2 = functions::nuP4(
      thetaH, nuHTauPhi, hTauNu.e(), hTauNu.P(), eX_htau, eY_htau, eZ_htau
    );

    nuLTauPhi = functions::phiFromLabMomenta(lTau, lTauLepton, beamAxis);
    lTauNu_2 = functions::nuP4(
      thetaL, nuLTauPhi, lTauNu.e(), lTauNu.P(), eX_ltau, eY_ltau, eZ_ltau
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
   * @brief Test whether reconstructed Z mass is good enough for the unit tests
   */
  void
  testZmass()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Z.mass(), constants::massZ, +1.);
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
    CPPUNIT_ASSERT_DOUBLES_EQUAL(lTauNu.e(), nuLtauEnergy, +1.e-3);
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
      lTauNu.e(), lTauNu.M2(), lTauNu.P(), lTauLepton.e(), lTauLepton.M2(), lTauLepton.P()
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
    CPPUNIT_ASSERT_DOUBLES_EQUAL(lTauNu_2.e(),  lTauNu.e(),  +1.e-3);
  }

  /**
   * @brief Test rotation angle correctness in leptonic tau decay system
   *        (neutrino px equality)
   */
  void
  testLtauPhi_px()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(lTauNu_2.px(), lTauNu.px(), +1.e-3);
  }

  /**
   * @brief Test rotation angle correctness in leptonic tau decay system
   *        (neutrino py equality)
   */
  void
  testLtauPhi_py()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(lTauNu_2.py(), lTauNu.py(), +1.e-3);
  }

  /**
   * @brief Test rotation angle correctness in leptonic tau decay system
   *        (neutrino pz equality)
   */
  void
  testLtauPhi_pz()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(lTauNu_2.pz(), lTauNu.pz(), +1.e-3);
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
