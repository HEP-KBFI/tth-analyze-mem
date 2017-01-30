#include "tthAnalysis/tthMEM/interface/MEMInterface_3l1tau.h" // MEMInterface_3l1tau, MEMOutput_3l1tau
#include "tthAnalysis/tthMEM/interface/Logger.h" // Logger::

#include <cppunit/extensions/HelperMacros.h> // CppUnit::TestFixture, CPPUNIT_ASSERT_*

class Test_MEMInterface_3l1tau
  : public CppUnit::TestFixture
{

  CPPUNIT_TEST_SUITE(Test_MEMInterface_3l1tau);
  CPPUNIT_TEST(test_very_quick);
  CPPUNIT_TEST(test_quick);
  CPPUNIT_TEST(test_medium);
  CPPUNIT_TEST(test_error_too_many_jets);
  CPPUNIT_TEST(test_error_too_few_jets);
  CPPUNIT_TEST(test_error_lepton_charge_sum);
  CPPUNIT_TEST_SUITE_END();

public:

  // declare the MEM environment
  MEMInterface_3l1tau mem;

  // declare the objects used to evaluate the MEM score
  std::vector<MeasuredJet> jets;
  MeasuredJet firstJet, secondJet, thirdJet;
  MeasuredLepton leadingLepton, subLeadingLepton, thirdLepton;
  MeasuredHadronicTau tau;
  MeasuredMET met;

  void
  setUp() override
  {
    Logger::enableLogging(false); // disable logging

//                               pT     eta     phi    mass
    firstJet  = MeasuredJet(184.900, -0.022,  0.460, 14.140);
    secondJet = MeasuredJet(286.300,  1.050, -2.045, 23.510);
    thirdJet  = MeasuredJet( 42.340, -1.028,  0.086,  7.765);
    jets = {firstJet, secondJet, thirdJet};

//                                         pT     eta     phi   mass charge
    leadingLepton    = MeasuredLepton(117.000,  1.140,  2.136, 0.106,    1);
    subLeadingLepton = MeasuredLepton(107.800, -0.225, -0.157, 0.106,   -1);
    thirdLepton      = MeasuredLepton( 89.290,  1.639, -1.706, 0.106,    1);

//                                pT    eta    phi   mass charge  decayMode
    tau = MeasuredHadronicTau(25.600, 0.410, 2.469, 0.855,    -1,        1);

//                         pT    phi   cov_XX   cov_XY   cov_YY
    met = MeasuredMET(164.900, 2.406, 671.900, 356.200, 775.100);
  }

  void
  tearDown() override
  {
    /* empty */
  }

  void
  test_very_quick()
  {
    // set MEM environment variables
    mem.maxCallsStartingPos = 100;
    mem.maxObjFunctionCalls = 1000;

    // initialize the MEM environment
    mem.initialize();

    // compute the MEM score
    const MEMOutput_3l1tau result = mem(jets, leadingLepton, subLeadingLepton, thirdLepton, tau, met);

    // assert
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.prob_tth,           5.8570e-55, 1.e-59);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.prob_tth_err,       1.7339e-55, 1.e-59);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.prob_tth_h2ww,     -1.0000e+00, 0.    );
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.prob_tth_h2ww_err, -1.0000e+00, 0.    );
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.prob_ttz,           1.7738e-51, 1.e-55);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.prob_ttz_err,       2.2831e-52, 1.e-56);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.sig,                5.8570e-55, 1.e-59);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.sig_err,            1.7339e-55, 1.e-59);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.bkg,                1.7738e-51, 1.e-55);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.bkg_err,            2.2831e-52, 1.e-56);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.lr,                 0.000330,   1.e-6 );
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.lr_err,             0.128623,   1.e-6 );
    CPPUNIT_ASSERT_EQUAL(result.err, 0);
  }

  void
  test_quick()
  {
    // set MEM environment variables
    mem.maxCallsStartingPos = 5000;
    mem.maxObjFunctionCalls = 2000;

    // initialize the MEM environment
    mem.initialize();

    // compute the MEM score
    const MEMOutput_3l1tau result = mem(jets, leadingLepton, subLeadingLepton, thirdLepton, tau, met);

    // assert
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.prob_tth,           8.4016e-53, 1.e-57);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.prob_tth_err,       1.6246e-53, 1.e-57);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.prob_tth_h2ww,     -1.0000e+00, 0.    );
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.prob_tth_h2ww_err, -1.0000e+00, 0.    );
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.prob_ttz,           6.6082e-51, 1.e-55);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.prob_ttz_err,       1.0005e-51, 1.e-55);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.sig,                8.4016e-53, 1.e-57);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.sig_err,            1.6246e-53, 1.e-57);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.bkg,                6.6082e-51, 1.e-55);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.bkg_err,            1.0005e-51, 1.e-55);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.lr,                 0.012554,   1.e-6 );
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.lr_err,             0.147621,   1.e-6 );
    CPPUNIT_ASSERT_EQUAL(result.err, 0);
  }

  void
  test_medium()
  {
    // set MEM environment variables
    mem.maxCallsStartingPos = 5000;
    mem.maxObjFunctionCalls = 10000;

    // initialize the MEM environment
    mem.initialize();

    // compute the MEM score
    const MEMOutput_3l1tau result = mem(jets, leadingLepton, subLeadingLepton, thirdLepton, tau, met);

    // assert
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.prob_tth,           1.6883e-51, 1.e-55);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.prob_tth_err,       2.6524e-52, 1.e-56);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.prob_tth_h2ww,     -1.0000e+00, 0.    );
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.prob_tth_h2ww_err, -1.0000e+00, 0.    );
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.prob_ttz,           6.6437e-51, 1.e-55);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.prob_ttz_err,       1.3602e-51, 1.e-55);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.sig,                1.6883e-51, 1.e-55);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.sig_err,            2.6524e-52, 1.e-56);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.bkg,                6.6437e-51, 1.e-55);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.bkg_err,            1.3602e-51, 1.e-55);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.lr,                 0.202628,   1.e-6 );
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result.lr_err,             0.130335,   1.e-6 );
    CPPUNIT_ASSERT_EQUAL(result.err, 0);
  }

  void
  test_error_lepton_charge_sum()
  {
    // flip lepton charge
    const MeasuredLepton subLeadingLepton_flip(subLeadingLepton.p4(), -subLeadingLepton.charge());

    // set MEM environment variables
    mem.maxCallsStartingPos = 100;
    mem.maxObjFunctionCalls = 1000;
    // initialize the MEM environment
    mem.initialize();

    // compute the MEM score
    const MEMOutput_3l1tau result = mem(jets, leadingLepton, subLeadingLepton_flip, thirdLepton, tau, met);

    // assert
    CPPUNIT_ASSERT_EQUAL(result.err, TTHEXCEPTION_ERR_CODE_NONZERO_CHARGE_SUM);
  }

  void
  test_error_too_many_jets()
  {
    // add a fourth jet
    MeasuredJet fourthJet(28.940, -0.037, 1.802, 6.749);
    const std::vector<MeasuredJet> jets_4 = { firstJet, secondJet, thirdJet, fourthJet };

    // set MEM environment variables
    mem.maxCallsStartingPos = 100;
    mem.maxObjFunctionCalls = 1000;

    // initialize the MEM environment
    mem.initialize();

    // compute the MEM score
    const MEMOutput_3l1tau result = mem(jets_4, leadingLepton, subLeadingLepton, thirdLepton, tau, met);

    // assert
    CPPUNIT_ASSERT_EQUAL(result.err, TTHEXCEPTION_ERR_CODE_INVALID_NOF_JETS);
  }

  void
  test_error_too_few_jets()
  {
    // use only one jet
    const std::vector<MeasuredJet> jets_1 = { firstJet };

    // set MEM environment variables
    mem.maxCallsStartingPos = 100;
    mem.maxObjFunctionCalls = 1000;

    // initialize the MEM environment
    mem.initialize();

    // compute the MEM score
    const MEMOutput_3l1tau result = mem(jets_1, leadingLepton, subLeadingLepton, thirdLepton, tau, met);

    // assert
    CPPUNIT_ASSERT_EQUAL(result.err, TTHEXCEPTION_ERR_CODE_INVALID_NOF_JETS);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(Test_MEMInterface_3l1tau);



