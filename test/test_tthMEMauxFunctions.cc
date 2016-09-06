#include <string> // std::string

#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // tthMEM::roundToNdigits()
#include <cppunit/extensions/HelperMacros.h> // CppUnit::TestFixture, CPPUNIT_ASSERT_DOUBLES_EQUAL,
                                             // CPPUNIT_ASSERT_ASSERTION_FAIL

class Test_tthMEMauxFunctions
  : public CppUnit::TestFixture
{

  CPPUNIT_TEST_SUITE(Test_tthMEMauxFunctions);
  CPPUNIT_TEST(testRoundToNdigits);
  CPPUNIT_TEST_SUITE_END();

public:

  void
  setUp() override
  {
    /* empty */
  }

  void
  tearDown() override
  {
    /* empty */
  }

  void
  testRoundToNdigits()
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tthMEM::roundToNdigits(123.456, 1), 100.000, 1e-8);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tthMEM::roundToNdigits(123.456, 2), 120.000, 1e-8);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tthMEM::roundToNdigits(123.456, 3), 123.000, 1e-8);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tthMEM::roundToNdigits(123.456, 4), 123.500, 1e-8);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tthMEM::roundToNdigits(123.456, 5), 123.460, 1e-8);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tthMEM::roundToNdigits(123.456, 6), 123.456, 1e-8);

    CPPUNIT_ASSERT_ASSERTION_FAIL(
      CPPUNIT_ASSERT_DOUBLES_EQUAL(tthMEM::roundToNdigits(123.456, 1), 123.456, 1e-8)
    );
    CPPUNIT_ASSERT_ASSERTION_FAIL(
      CPPUNIT_ASSERT_DOUBLES_EQUAL(tthMEM::roundToNdigits(123.456, 2), 123.456, 1e-8)
    );
    CPPUNIT_ASSERT_ASSERTION_FAIL(
      CPPUNIT_ASSERT_DOUBLES_EQUAL(tthMEM::roundToNdigits(123.456, 3), 123.456, 1e-8)
    );
    CPPUNIT_ASSERT_ASSERTION_FAIL(
      CPPUNIT_ASSERT_DOUBLES_EQUAL(tthMEM::roundToNdigits(123.456, 4), 123.456, 1e-8)
    );
    CPPUNIT_ASSERT_ASSERTION_FAIL(
      CPPUNIT_ASSERT_DOUBLES_EQUAL(tthMEM::roundToNdigits(123.456, 5), 123.456, 1e-8)
    );
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(Test_tthMEMauxFunctions);

