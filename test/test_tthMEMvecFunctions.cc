#include <string> // std::string
#include <vector> // std::vector<>
#include <algorithm> // std::equal()
#include <cmath> // std::sqrt()

#include <cppunit/extensions/HelperMacros.h> // CppUnit::TestFixture, CPPUNIT_ASSERT_*

#include "tthAnalysis/tthMEM/interface/tthMEMvecFunctions.h" // tthMEM::vec::, tthMEMexception

// enable operator overloading for std::vector<>'s that we unit-test here
using namespace tthMEM;

class Test_tthMEMvecFunctions
  : public CppUnit::TestFixture
{

  CPPUNIT_TEST_SUITE(Test_tthMEMvecFunctions);
  CPPUNIT_TEST(testSubv);
  CPPUNIT_TEST(testAvg);
  CPPUNIT_TEST(testL2);
  CPPUNIT_TEST(testStdev);
  CPPUNIT_TEST(testVectorHadamardMultiplication);
  CPPUNIT_TEST(testVectorAddition);
  CPPUNIT_TEST(testVectorSubtraction);
  CPPUNIT_TEST(testVectorScalarAddition);
  CPPUNIT_TEST(testVectorScalarSubtraction);
  CPPUNIT_TEST(testVectorScalarMultiplication);
  CPPUNIT_TEST(testVectorScalarDivision);
  CPPUNIT_TEST(testVectorScalarAdditionAssignment);
  CPPUNIT_TEST(testVectorScalarSubtractionAssignment);
  CPPUNIT_TEST(testVectorScalarMultiplicationAssignment);
  CPPUNIT_TEST(testVectorScalarDivisionAssignment);
  CPPUNIT_TEST(testVectorAdditionAssignment);
  CPPUNIT_TEST(testVectorSubtractionAssignment);
  CPPUNIT_TEST_SUITE_END();

public:

  /* for function operating on (sub)vectors */
  const std::vector<double> simpleVector   { +1., +2., +3., +4., +5. };
  const std::vector<double> subVector0to3  { +1., +2., +3.           };
  const std::vector<double> subVector2toEnd{           +3., +4., +5. };
  const std::vector<double> subVector1to4  {      +2., +3., +4.      };
  const double subVector0to3Avg   = +2.0;
  const double subVector2toEndAvg = +4.0;
  const double subVector1to4Avg   = +3.0;
  const double subVector0to3L2        = std::sqrt(14.);
  const double subVector0to3L2ShiftM1 = std::sqrt(5.);
  const double subVector2toEndL2        = std::sqrt(50.);
  const double subVector2toEndL2ShiftM1 = std::sqrt(29.);
  const double subVector1to4L2        = std::sqrt(29.);
  const double subVector1to4L2ShiftM1 = std::sqrt(14.);
  const double subVector0to3Stdev   = std::sqrt(1. / 3);
  const double subVector2toEndStdev = std::sqrt(1. / 3);
  const double subVector1to4Stdev   = std::sqrt(1. / 3);
  const unsigned simpleVectorSize = simpleVector.size();

  /* general purpose vectors */
  const std::vector<double> emptyVector  {     };
  const std::vector<double> vectorWElem0 {  0. };
  const std::vector<double> vectorWElem1 { +1. };
  const std::vector<double> vectorWElemM1{ -1. };

  /* general purpose scalars (prevent implicit casts */
  const int Int0  =  0;
  const int Int1  = +1;
  const int IntM1 = -1;
  const int Int2  = +2;
  const int IntM2 = -2;
  const double Double0  =  0.;
  const double Double1  = +1.;
  const double DoubleM1 = -1.;
  const double Double2  = +2.;
  const double DoubleM2 = -2.;

  /* vectors for testing Hadamard product operation */
  const std::vector<double> vectorProdOperand0{ -2., -1.,  0., +1., +2., +3. };
  const std::vector<double> vectorProdOperand1{ +1.,  0., -1., -2., +2., +3. };
  const std::vector<double> vectorProdResult  { -2.,  0.,  0., -2., +4., +9. };

  /* vectors for testing addition operation */
  const std::vector<double> vectorAddOperand0{ -2., -1., 0., +1., +2., +3. };
  const std::vector<double> vectorAddOperand1{ +2., +1., 0., -1., -2., -3. };
  const std::vector<double> vectorAddResult  {  0.,  0., 0.,  0.,  0.,  0. };

  /* vectors for testing subtraction operation */
  const std::vector<double> vectorSubOperand0{ +1., +2., +3., +4., +5., +6. };
  const std::vector<double> vectorSubOperand1{ +1., -1., +1., -1., +1., -1. };
  const std::vector<double> vectorSubResult  {  0., +3., +2., +5., +4., +7. };
  const std::vector<double> vectorSubResultR {  0., -3., -2., -5., -4., -7. };
  const std::vector<double> vectorSubEmpty   {  0.,  0.,  0.,  0.,  0.,  0. };

  /* vectors for testing scalar-vector addition */
  const std::vector<double> vectorScalarAdd      { -2., -1.,  0., +1., +2., +3. };
  const std::vector<double> vectorScalarAddResult{  0., +1., +2., +3., +4., +5. }; // add 2

  /* vectors for testing scalar-vector subtraction */
  const std::vector<double> vectorScalarSub       { -2., -1.,  0., +1., +2., +3. };
  const std::vector<double> vectorScalarSubResult { -4., -3., -2., -1.,  0., +1. }; // subtract 2
  const std::vector<double> vectorScalarSubResultR{ +4., +3., +2., +1.,  0., -1. }; // subtract from 2

  /* vectors for testing scalar-vector multiplication */
  const std::vector<double> vectorScalarProd      { -2., -1., 0., +1., +2., +3. };
  const std::vector<double> vectorScalarProdResult{ +4., +2., 0., -2., -4., -6. }; // multiply by -2

  /* vectors for testing scalar-vector division */
  const std::vector<double> vectorScalarDiv        { -2., -1.0, 0., +1.0, +2., +4. };
  const std::vector<double> vectorScalarDivResult  { +1., +0.5, 0., -0.5, -1., -2. }; // divide by -2
  const std::vector<double> vectorScalarDiv2       { -2., -1.0, +1.0, +2., +4.0 };
  const std::vector<double> vectorScalarDiv2Result { +1., +0.5, -0.5, -1., -2.0 }; // divide by -2
  const std::vector<double> vectorScalarDiv2ResultR{ +1., +2.0, -2.0, -1., -0.5 }; // inverse of previous operation

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

  /**
   * @brief tests:
   *
   * std::vector<VectorElementType>
   * subv(const std::vector<VectorElementType> & vec,
   *      unsigned shiftFromBegin_begin,
   *      unsigned shiftFromBegin_end);
   */
  void
  testSubv()
  {
    /* test with conflicting shifts */
    CPPUNIT_ASSERT_THROW(vec::subv(simpleVector, 1,                       0                      ), tthMEMexception);
    CPPUNIT_ASSERT_THROW(vec::subv(simpleVector, 0,                       simpleVector.size() + 1), tthMEMexception);
    CPPUNIT_ASSERT_THROW(vec::subv(simpleVector, simpleVector.size() + 1, simpleVector.size() + 2), tthMEMexception);
    /* test actual functionality */
    CPPUNIT_ASSERT_NO_THROW(vec::subv(simpleVector, 0, 3));
    const auto simpleVector0to3 = vec::subv(simpleVector, 0, 3);
    CPPUNIT_ASSERT(std::equal(simpleVector0to3.begin(), simpleVector0to3.end(), subVector0to3.begin()));
    CPPUNIT_ASSERT_NO_THROW(vec::subv(simpleVector, 2, simpleVector.size()));
    const auto simpleVector2toEnd = vec::subv(simpleVector, 2, simpleVector.size());
    CPPUNIT_ASSERT(std::equal(simpleVector2toEnd.begin(), simpleVector2toEnd.end(), subVector2toEnd.begin()));
    CPPUNIT_ASSERT_NO_THROW(vec::subv(simpleVector, 1, 4));
    const auto simpleVector1to4 = vec::subv(simpleVector, 1, 4);
    CPPUNIT_ASSERT(std::equal(simpleVector1to4.begin(), simpleVector1to4.end(), subVector1to4.begin()));
  }

  /**
   * @brief tests:
   *
   * double
   * avg(const std::vector<double> & vec,
   *     unsigned shiftFromBegin_begin,
   *     unsigned shiftFromBegin_end);
   */
  void
  testAvg()
  {
    /* test with conflicting shifts */
    CPPUNIT_ASSERT_THROW(vec::avg(simpleVector, 1,                    0                   ), tthMEMexception);
    CPPUNIT_ASSERT_THROW(vec::avg(simpleVector, 0,                    simpleVectorSize + 1), tthMEMexception);
    CPPUNIT_ASSERT_THROW(vec::avg(simpleVector, simpleVectorSize + 1, simpleVectorSize + 2), tthMEMexception);
    /* test actual functionality */
    CPPUNIT_ASSERT_NO_THROW(     vec::avg(simpleVector, 0, 3));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(vec::avg(simpleVector, 0, 3),                subVector0to3Avg,   +1.e-8);
    CPPUNIT_ASSERT_NO_THROW(     vec::avg(simpleVector, 2, simpleVectorSize));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(vec::avg(simpleVector, 2, simpleVectorSize), subVector2toEndAvg, +1.e-8);
    CPPUNIT_ASSERT_NO_THROW(     vec::avg(simpleVector, 1, 4));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(vec::avg(simpleVector, 1, 4),                subVector1to4Avg,   +1.e-8);
  }

  /**
   * @brief tests:
   *
   * double
   * l2(const std::vector<double> & vec,
   *    unsigned shiftFromBegin_begin,
   *    unsigned shiftFromBegin_end,
   *    double shiftValue = 0.);
   */
  void
  testL2()
  {
    /* test with conflicting shifts */
    CPPUNIT_ASSERT_THROW(vec::l2(simpleVector, 1u,                    0u                   ), tthMEMexception);
    CPPUNIT_ASSERT_THROW(vec::l2(simpleVector, 0u,                    simpleVectorSize + 1u), tthMEMexception);
    CPPUNIT_ASSERT_THROW(vec::l2(simpleVector, simpleVectorSize + 1u, simpleVectorSize + 2u), tthMEMexception);
    /* test actual functionality */
    CPPUNIT_ASSERT_NO_THROW(     vec::l2(simpleVector, 0u, 3u));
    CPPUNIT_ASSERT_NO_THROW(     vec::l2(simpleVector, 0u, 3u, -1.));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(vec::l2(simpleVector, 0u, 3u),      subVector0to3L2,        +1.e-8);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(vec::l2(simpleVector, 0u, 3u, -1.), subVector0to3L2ShiftM1, +1.e-8);
    CPPUNIT_ASSERT_NO_THROW(     vec::l2(simpleVector, 2u, simpleVectorSize));
    CPPUNIT_ASSERT_NO_THROW(     vec::l2(simpleVector, 2u, simpleVectorSize, -1.));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(vec::l2(simpleVector, 2u, simpleVectorSize),      subVector2toEndL2,        +1.e-8);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(vec::l2(simpleVector, 2u, simpleVectorSize, -1.), subVector2toEndL2ShiftM1, +1.e-8);
    CPPUNIT_ASSERT_NO_THROW(     vec::l2(simpleVector, 1u, 4u));
    CPPUNIT_ASSERT_NO_THROW(     vec::l2(simpleVector, 1u, 4u, -1.));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(vec::l2(simpleVector, 1u, 4u),      subVector1to4L2,        +1.e-8);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(vec::l2(simpleVector, 1u, 4u, -1.), subVector1to4L2ShiftM1, +1.e-8);
  }

  /**
   * @brief tests:
   *
   * double
   * stdev(const std::vector<double> & vec,
   *       unsigned shiftFromBegin_begin,
   *       unsigned shiftFromBegin_end,
   *       double average);
   */
  void
  testStdev()
  {
    /* test with conflicting shifts */
    CPPUNIT_ASSERT_THROW(vec::stdev(simpleVector, 1,                    0                   , subVector0to3Avg  ), tthMEMexception);
    CPPUNIT_ASSERT_THROW(vec::stdev(simpleVector, 0,                    simpleVectorSize + 1, subVector2toEndAvg), tthMEMexception);
    CPPUNIT_ASSERT_THROW(vec::stdev(simpleVector, simpleVectorSize + 1, simpleVectorSize + 2, subVector1to4Avg  ), tthMEMexception);
    /* test actual functionality */
    CPPUNIT_ASSERT_NO_THROW(     vec::stdev(simpleVector, 0, 3,                subVector0to3Avg));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(vec::stdev(simpleVector, 0, 3,                subVector0to3Avg),   subVector0to3Stdev,    +1.e-8);
    CPPUNIT_ASSERT_NO_THROW(     vec::stdev(simpleVector, 2, simpleVectorSize, subVector2toEndAvg));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(vec::stdev(simpleVector, 2, simpleVectorSize, subVector2toEndAvg), subVector2toEndStdev, +1.e-8);
    CPPUNIT_ASSERT_NO_THROW(     vec::stdev(simpleVector, 1, 4,                subVector1to4Avg));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(vec::stdev(simpleVector, 1, 4,                subVector1to4Avg),   subVector1to4Stdev,   +1.e-8);
  }

  /**
   * @brief tests:
   *
   * std::vector<VectorElementType>
   * operator*(const std::vector<VectorElementType> & lhs,
   *           const std::vector<VectorElementType> & rhs);
   */
  void
  testVectorHadamardMultiplication()
  {
    /* multiply empty vectors */
    CPPUNIT_ASSERT_NO_THROW(emptyVector * emptyVector);
    const auto emptyVectorProduct = emptyVector * emptyVector;
    CPPUNIT_ASSERT(std::equal(emptyVectorProduct.begin(), emptyVectorProduct.end(), emptyVector.begin()));
    /* multiply empty vector with a non-empty one */
    CPPUNIT_ASSERT_THROW(emptyVector * vectorWElem0, tthMEMexception);
    /* simple 1-element product: { 0 } x { 1 } = { 0 } */
    CPPUNIT_ASSERT_NO_THROW(vectorWElem0 * vectorWElem1);
    const auto vectorProduct0Times1 = vectorWElem0 * vectorWElem1;
    CPPUNIT_ASSERT_NO_THROW(std::equal(vectorProduct0Times1.begin(), vectorProduct0Times1.end(), vectorWElem0.begin()));
    /* simple 1-element product: { 1 } x { 1 } = { 1 } */
    CPPUNIT_ASSERT_NO_THROW(vectorWElem1 * vectorWElem1);
    const auto vectorProduct1Times1 = vectorWElem1 * vectorWElem1;
    CPPUNIT_ASSERT_NO_THROW(std::equal(vectorProduct1Times1.begin(), vectorProduct1Times1.end(), vectorWElem1.begin()));
    /* simple 1-element product: { 1 } x { -1 } = { -1 } */
    CPPUNIT_ASSERT_NO_THROW(vectorWElem1 * vectorWElemM1);
    const auto vectorProduct1TimesM1 = vectorWElem1 * vectorWElemM1;
    CPPUNIT_ASSERT_NO_THROW(std::equal(vectorProduct1TimesM1.begin(), vectorProduct1TimesM1.end(), vectorWElemM1.begin()));
    /* simple 1-element product: { -1 } x { 1 } = { -1 } */
    CPPUNIT_ASSERT_NO_THROW(vectorWElemM1 * vectorWElem1);
    const auto vectorProductM1Times1 = vectorWElemM1 * vectorWElem1;
    CPPUNIT_ASSERT_NO_THROW(std::equal(vectorProductM1Times1.begin(), vectorProductM1Times1.end(), vectorWElemM1.begin()));
    /* complex multiplication */
    CPPUNIT_ASSERT_NO_THROW(vectorProdOperand0 * vectorProdOperand1);
    const auto vectorProdActualResult = vectorProdOperand0 * vectorProdOperand1;
    CPPUNIT_ASSERT(std::equal(vectorProdActualResult.begin(), vectorProdActualResult.end(), vectorProdResult.begin()));
    /* complex multiplication: commutativity test */
    CPPUNIT_ASSERT_NO_THROW(vectorProdOperand1 * vectorProdOperand0);
    const auto vectorProdActualResultR = vectorProdOperand1 * vectorProdOperand0;
    CPPUNIT_ASSERT(std::equal(vectorProdActualResultR.begin(), vectorProdActualResultR.end(), vectorProdResult.begin()));
  }

  /**
   * @brief tests:
   *
   * std::vector<VectorElementType>
   * operator+(const std::vector<VectorElementType> & lhs,
   *           const std::vector<VectorElementType> & rhs);
   */
  void
  testVectorAddition()
  {
    /* add empty vectors: { } + { } = { } */
    CPPUNIT_ASSERT_NO_THROW(emptyVector + emptyVector);
    const auto emptyVectorSum = emptyVector + emptyVector;
    CPPUNIT_ASSERT(std::equal(emptyVectorSum.begin(), emptyVectorSum.end(), emptyVector.begin()));
    /* add empty vector w/ non-empty one: { } + { 0 } */
    CPPUNIT_ASSERT_THROW(emptyVector + vectorWElem0, tthMEMexception);
    /* simple 1-element addition: { 0 } + { 1 } = { 1 } */
    CPPUNIT_ASSERT_NO_THROW(vectorWElem0 + vectorWElem1);
    const auto vectorSum0plus1 = vectorWElem0 + vectorWElem1;
    CPPUNIT_ASSERT(std::equal(vectorSum0plus1.begin(), vectorSum0plus1.end(), vectorWElem1.begin()));
    /* complex multi-element addition */
    CPPUNIT_ASSERT_NO_THROW(vectorAddOperand0 + vectorAddOperand1);
    const auto vectorAddActualResult = vectorAddOperand0 + vectorAddOperand1;
    CPPUNIT_ASSERT(std::equal(vectorAddActualResult.begin(), vectorAddActualResult.end(), vectorAddResult.begin()));
    /* the same as the previous assertion, but in reverse order (tests commutation) */
    CPPUNIT_ASSERT_NO_THROW(vectorAddOperand1 + vectorAddOperand0);
    const auto vectorAddActualResultR = vectorAddOperand1 + vectorAddOperand0;
    CPPUNIT_ASSERT(std::equal(vectorAddActualResultR.begin(), vectorAddActualResultR.end(), vectorAddResult.begin()));
  }

  /**
   * @brief tests:
   *
   * std::vector<VectorElementType>
   * operator-(const std::vector<VectorElementType> & lhs,
   *           const std::vector<VectorElementType> & rhs);
   */
  void
  testVectorSubtraction()
  {
    /* subtract empty vectors: { } - { } = { } */
    CPPUNIT_ASSERT_NO_THROW(emptyVector - emptyVector);
    const auto emptyVectorDiff = emptyVector - emptyVector;
    CPPUNIT_ASSERT(std::equal(emptyVectorDiff.begin(), emptyVectorDiff.end(), emptyVector.begin()));
    /* subtract non-empty vector from an empty one: { } - { 0 } */
    CPPUNIT_ASSERT_THROW(emptyVector - vectorWElem0, tthMEMexception);
    /* subtract 1-element vector from itself: { 0 } - { 0 } = { } */
    CPPUNIT_ASSERT_NO_THROW(vectorWElem0 - vectorWElem0);
    const auto vectorSub0minus0 = vectorWElem0 - vectorWElem0;
    CPPUNIT_ASSERT(std::equal(vectorSub0minus0.begin(), vectorSub0minus0.end(), vectorWElem0.begin()));
    /* subtract two different 1-element vectors: { 0 } - { 1 } = { -1 } */
    CPPUNIT_ASSERT_NO_THROW(vectorWElem0 - vectorWElem1);
    const auto vectorSub0minus1 = vectorWElem0 - vectorWElem1;
    CPPUNIT_ASSERT(std::equal(vectorSub0minus1.begin(), vectorSub0minus1.end(), vectorWElemM1.begin()));
    /* subtract two different 1-element vectors, but in reverse order: { 1 } - { 0 } = { 1 } */
    CPPUNIT_ASSERT_NO_THROW(vectorWElem1 - vectorWElem0);
    const auto vectorSub1minus0 = vectorWElem1 - vectorWElem0;
    CPPUNIT_ASSERT(std::equal(vectorSub1minus0.begin(), vectorSub1minus0.end(), vectorWElem1.begin()));
    /* more complex vector subtrcation: large empty check */
    CPPUNIT_ASSERT_NO_THROW(vectorSubOperand0 - vectorSubOperand0);
    const auto vectorSubOperand0minus0 = vectorSubOperand0 - vectorSubOperand0;
    CPPUNIT_ASSERT(std::equal(vectorSubOperand0minus0.begin(), vectorSubOperand0minus0.end(), vectorSubEmpty.begin()));
    /* more complex vector subtrcation: different operands */
    CPPUNIT_ASSERT_NO_THROW(vectorSubOperand0 - vectorSubOperand1);
    const auto vectorSubOperand0minus1 = vectorSubOperand0 - vectorSubOperand1;
    CPPUNIT_ASSERT(std::equal(vectorSubOperand0minus1.begin(), vectorSubOperand0minus1.end(), vectorSubResult.begin()));
    /* more complex vector subtrcation: different operands in reverse order */
    CPPUNIT_ASSERT_NO_THROW(vectorSubOperand1 - vectorSubOperand0);
    const auto vectorSubOperand1minus0 = vectorSubOperand1 - vectorSubOperand0;
    CPPUNIT_ASSERT(std::equal(vectorSubOperand1minus0.begin(), vectorSubOperand1minus0.end(), vectorSubResultR.begin()));
  }

  /**
   * @brief tests:
   *
   * std::vector<VectorElementType>
   * operator+(ScalarType lhs,
   *           const std::vector<VectorElementType> & rhs);
   *
   * std::vector<VectorElementType>
   * operator+(const std::vector<VectorElementType> & lhs,
   *           ScalarType rhs);
   */
  void
  testVectorScalarAddition()
  {
    /* adding a scalar to an empty vector gives us an empty vector: { } + 1 = { } */
    /* scalar is an int */
    CPPUNIT_ASSERT_NO_THROW(emptyVector + Int1);
    const auto emptyVectorPlus1Int = emptyVector + Int1;
    CPPUNIT_ASSERT(std::equal(emptyVectorPlus1Int.begin(), emptyVectorPlus1Int.end(), emptyVector.begin()));
    /* scalar is a double */
    CPPUNIT_ASSERT_NO_THROW(emptyVector + Double1);
    const auto emptyVectorPlus1Double = emptyVector + Double1;
    CPPUNIT_ASSERT(std::equal(emptyVectorPlus1Double.begin(), emptyVectorPlus1Double.end(), emptyVector.begin()));
    /* scalar is an int, reverse order */
    CPPUNIT_ASSERT_NO_THROW(Int1 + emptyVector);
    const auto Int1PlusEmptyVector = Int1 + emptyVector;
    CPPUNIT_ASSERT(std::equal(Int1PlusEmptyVector.begin(), Int1PlusEmptyVector.end(), emptyVector.begin()));
    /* scalar is a double, reverse order */
    CPPUNIT_ASSERT_NO_THROW(Double1 + emptyVector);
    const auto Double1PlusEmptyVector = Double1 + emptyVector;
    CPPUNIT_ASSERT(std::equal(Double1PlusEmptyVector.begin(), Double1PlusEmptyVector.end(), emptyVector.begin()));
    /* add scalar to a simple 1-element vector (int): { 0 } + 1 = { 1 } */
    CPPUNIT_ASSERT_NO_THROW(vectorWElem0 + Int1);
    const auto vectorAdd0Vplus1Sint = vectorWElem0 + Int1;
    CPPUNIT_ASSERT(std::equal(vectorAdd0Vplus1Sint.begin(), vectorAdd0Vplus1Sint.end(), vectorWElem1.begin()));
    /* add scalar to a simple 1-element vector (int, reverse): 1 + { 0 } = { 1 } */
    CPPUNIT_ASSERT_NO_THROW(Int1 + vectorWElem0);
    const auto vectorAdd0Vplus1SintR = Int1 + vectorWElem0;
    CPPUNIT_ASSERT(std::equal(vectorAdd0Vplus1SintR.begin(), vectorAdd0Vplus1SintR.end(), vectorWElem1.begin()));
    /* add scalar to a simple 1-element vector (double): { 0 } + 1 = { 1 } */
    CPPUNIT_ASSERT_NO_THROW(vectorWElem0 + Double1);
    const auto vectorAdd0Vplus1Sdouble = vectorWElem0 + Double1;
    CPPUNIT_ASSERT(std::equal(vectorAdd0Vplus1Sdouble.begin(), vectorAdd0Vplus1Sdouble.end(), vectorWElem1.begin()));
    /* add scalar to a simple 1-element vector (double, reverse): 1 + { 0 } = { 1 } */
    CPPUNIT_ASSERT_NO_THROW(Double1 + vectorWElem0);
    const auto vectorAdd0Vplus1SdoubleR = Double1 + vectorWElem0;
    CPPUNIT_ASSERT(std::equal(vectorAdd0Vplus1SdoubleR.begin(), vectorAdd0Vplus1SdoubleR.end(), vectorWElem1.begin()));
    /* complex addition (int) */
    CPPUNIT_ASSERT_NO_THROW(vectorScalarAdd + Int2);
    const auto vectorScalarAddResultInt = vectorScalarAdd + Int2;
    CPPUNIT_ASSERT(std::equal(vectorScalarAddResultInt.begin(), vectorScalarAddResultInt.end(), vectorScalarAddResult.begin()));
    /* complex addition (int, reverse) */
    CPPUNIT_ASSERT_NO_THROW(Int2 + vectorScalarAdd);
    const auto vectorScalarAddResultIntR = Int2 + vectorScalarAdd;
    CPPUNIT_ASSERT(std::equal(vectorScalarAddResultIntR.begin(), vectorScalarAddResultIntR.end(), vectorScalarAddResult.begin()));
    /* complex addition (double) */
    CPPUNIT_ASSERT_NO_THROW(vectorScalarAdd + Double2);
    const auto vectorScalarAddResultDouble = vectorScalarAdd + Double2;
    CPPUNIT_ASSERT(std::equal(vectorScalarAddResultDouble.begin(), vectorScalarAddResultDouble.end(), vectorScalarAddResult.begin()));
    /* complex addition (double, reverse) */
    CPPUNIT_ASSERT_NO_THROW(Double2 + vectorScalarAdd);
    const auto vectorScalarAddResultDoubleR = Double2 + vectorScalarAdd;
    CPPUNIT_ASSERT(std::equal(vectorScalarAddResultDoubleR.begin(), vectorScalarAddResultDoubleR.end(), vectorScalarAddResult.begin()));
  }

  /**
   * @brief tests:
   *
   * std::vector<VectorElementType>
   * operator-(ScalarType lhs,
   *           const std::vector<VectorElementType> & rhs);
   *
   * std::vector<VectorElementType>
   * operator-(const std::vector<VectorElementType> & lhs,
   *           ScalarType rhs);
   */
  void
  testVectorScalarSubtraction()
  {
    /* subtracting a scalar from an empty vector gives us an empty vector: { } - 1 = { } */
    /* scalar is an int */
    CPPUNIT_ASSERT_NO_THROW(emptyVector - Int1);
    const auto emptyVectorMinus1Int = emptyVector - Int1;
    CPPUNIT_ASSERT(std::equal(emptyVectorMinus1Int.begin(), emptyVectorMinus1Int.end(), emptyVector.begin()));
    /* scalar is a double */
    CPPUNIT_ASSERT_NO_THROW(emptyVector - Double1);
    const auto emptyVectorMinus1Double = emptyVector - Double1;
    CPPUNIT_ASSERT(std::equal(emptyVectorMinus1Double.begin(), emptyVectorMinus1Double.end(), emptyVector.begin()));
    /* subtracting an empty vector from a scalar gives us an empty vector: 1 - { } = { } */
    CPPUNIT_ASSERT_NO_THROW(Int1 - emptyVector);
    const auto emptyVectorMinus1IntR = Int1 - emptyVector;
    CPPUNIT_ASSERT(std::equal(emptyVectorMinus1IntR.begin(), emptyVectorMinus1IntR.end(), emptyVector.begin()));
    /* scalar is a double */
    CPPUNIT_ASSERT_NO_THROW(Double1 - emptyVector);
    const auto emptyVectorMinus1DoubleR = Double1 - emptyVector;
    CPPUNIT_ASSERT(std::equal(emptyVectorMinus1DoubleR.begin(), emptyVectorMinus1DoubleR.end(), emptyVector.begin()));
    /* subtract a scalar from a simple 1-element vector (int): { 0 } - 1 = { -1 } */
    CPPUNIT_ASSERT_NO_THROW(vectorWElem0 - Int1);
    const auto vectorSubVminus1Sint = vectorWElem0 - Int1;
    CPPUNIT_ASSERT(std::equal(vectorSubVminus1Sint.begin(), vectorSubVminus1Sint.end(), vectorWElemM1.begin()));
    /* subtract a simple 1-element vector from a scalar (int): 1 - { 0 } = { 1 } */
    CPPUNIT_ASSERT_NO_THROW(Int1 - vectorWElem0);
    const auto vectorSubVminus1SintR = Int1 - vectorWElem0;
    CPPUNIT_ASSERT(std::equal(vectorSubVminus1SintR.begin(), vectorSubVminus1SintR.end(), vectorWElem1.begin()));
    /* subtract a scalar from a simple 1-element vector (double): { 0 } - 1 = { -1 } */
    CPPUNIT_ASSERT_NO_THROW(vectorWElem0 - Double1);
    const auto vectorSubVminus1Sdouble = vectorWElem0 - Double1;
    CPPUNIT_ASSERT(std::equal(vectorSubVminus1Sdouble.begin(), vectorSubVminus1Sdouble.end(), vectorWElemM1.begin()));
    /* subtract a simple 1-element vector from a scalar (double): 1 - { 0 } = { 1 } */
    CPPUNIT_ASSERT_NO_THROW(Double1 - vectorWElem0);
    const auto vectorSubVminus1SdoubleR = Double1 - vectorWElem0;
    CPPUNIT_ASSERT(std::equal(vectorSubVminus1SdoubleR.begin(), vectorSubVminus1SdoubleR.end(), vectorWElem1.begin()));
    /* complex subtraction (int) */
    CPPUNIT_ASSERT_NO_THROW(vectorScalarSub - Int2);
    const auto vectorSubVminus2Int = vectorScalarSub - Int2;
    CPPUNIT_ASSERT(std::equal(vectorSubVminus2Int.begin(), vectorSubVminus2Int.end(), vectorScalarSubResult.begin()));
    /* complex subtraction (int, reverse) */
    CPPUNIT_ASSERT_NO_THROW(Int2 - vectorScalarSub);
    const auto vectorSubVminus2IntR = Int2 - vectorScalarSub;
    CPPUNIT_ASSERT(std::equal(vectorSubVminus2IntR.begin(), vectorSubVminus2IntR.end(), vectorScalarSubResultR.begin()));
    /* complex subtraction (double) */
    CPPUNIT_ASSERT_NO_THROW(vectorScalarSub - Double2);
    const auto vectorSubVminus2Double = vectorScalarSub - Double2;
    CPPUNIT_ASSERT(std::equal(vectorSubVminus2Double.begin(), vectorSubVminus2Double.end(), vectorScalarSubResult.begin()));
    /* complex subtraction (double, reverse) */
    CPPUNIT_ASSERT_NO_THROW(Double2 - vectorScalarSub);
    const auto vectorSubVminus2DoubleR = Double2 - vectorScalarSub;
    CPPUNIT_ASSERT(std::equal(vectorSubVminus2DoubleR.begin(), vectorSubVminus2DoubleR.end(), vectorScalarSubResultR.begin()));
  }

  /**
   * @brief tests:
   *
   * std::vector<VectorElementType>
   * operator*(ScalarType lhs,
   *           const std::vector<VectorElementType> & rhs);
   *
   * std::vector<VectorElementType>
   * operator*(const std::vector<VectorElementType> & lhs,
   *           ScalarType rhs);
   */
  void
  testVectorScalarMultiplication()
  {
    /* multiplying a scalar with an empty vector gives us an empty vector: { } x 1 = { } */
    /* scalar is an int */
    CPPUNIT_ASSERT_NO_THROW(emptyVector * Int1);
    const auto emptyVectorTimes1Int = emptyVector * Int1;
    CPPUNIT_ASSERT(std::equal(emptyVectorTimes1Int.begin(), emptyVectorTimes1Int.end(), emptyVector.begin()));
    /* scalar is a double */
    CPPUNIT_ASSERT_NO_THROW(emptyVector * Double1);
    const auto emptyVectorTimes1Double = emptyVector * Double1;
    CPPUNIT_ASSERT(std::equal(emptyVectorTimes1Double.begin(), emptyVectorTimes1Double.end(), emptyVector.begin()));
    /* scalar is an int, reverse order */
    CPPUNIT_ASSERT_NO_THROW(Int1 * emptyVector);
    const auto emptyVectorTimes1IntR = Int1 * emptyVector;
    CPPUNIT_ASSERT(std::equal(emptyVectorTimes1IntR.begin(), emptyVectorTimes1IntR.end(), emptyVector.begin()));
    /* scalar is a double, reverse order */
    CPPUNIT_ASSERT_NO_THROW(Double1 * emptyVector);
    const auto emptyVectorTimes1DoubleR = Double1 * emptyVector;
    CPPUNIT_ASSERT(std::equal(emptyVectorTimes1DoubleR.begin(), emptyVectorTimes1DoubleR.end(), emptyVector.begin()));
    /* multiply scalar with a simple 1-element vector (int): { 1 } x (-1) = { -1 } */
    CPPUNIT_ASSERT_NO_THROW(vectorWElem1 * IntM1);
    const auto vectorProdV1timesM1Sint = vectorWElem1 * IntM1;
    CPPUNIT_ASSERT(std::equal(vectorProdV1timesM1Sint.begin(), vectorProdV1timesM1Sint.end(), vectorWElemM1.begin()));
    /* multiply scalar with a simple 1-element vector (double): { 1 } x (-1) = { -1 } */
    CPPUNIT_ASSERT_NO_THROW(vectorWElem1 * DoubleM1);
    const auto vectorProdV1timesM1Sdouble = vectorWElem1 * DoubleM1;
    CPPUNIT_ASSERT(std::equal(vectorProdV1timesM1Sdouble.begin(), vectorProdV1timesM1Sdouble.end(), vectorWElemM1.begin()));
    /* multiply scalar with a simple 1-element vector (int, reverse): (-1) x { 1 } = { -1 } */
    CPPUNIT_ASSERT_NO_THROW(IntM1 * vectorWElem1);
    const auto vectorProdV1timesM1SintR = IntM1 * vectorWElem1;
    CPPUNIT_ASSERT(std::equal(vectorProdV1timesM1SintR.begin(), vectorProdV1timesM1SintR.end(), vectorWElemM1.begin()));
    /* multiply scalar with a simple 1-element vector (double, reverse): (-1) x { 1 } = { -1 } */
    CPPUNIT_ASSERT_NO_THROW(DoubleM1 * vectorWElem1);
    const auto vectorProdV1timesM1SdoubleR = DoubleM1 * vectorWElem1;
    CPPUNIT_ASSERT(std::equal(vectorProdV1timesM1SdoubleR.begin(), vectorProdV1timesM1SdoubleR.end(), vectorWElemM1.begin()));
    /* complex multiplication (int) */
    CPPUNIT_ASSERT_NO_THROW(vectorScalarProd * IntM2);
    const auto vectorScalarProdResultInt = vectorScalarProd * IntM2;
    CPPUNIT_ASSERT(std::equal(vectorScalarProdResultInt.begin(), vectorScalarProdResultInt.end(), vectorScalarProdResult.begin()));
    /* complex multiplication (double) */
    CPPUNIT_ASSERT_NO_THROW(vectorScalarProd * DoubleM2);
    const auto vectorScalarProdResultDouble = vectorScalarProd * DoubleM2;
    CPPUNIT_ASSERT(std::equal(vectorScalarProdResultDouble.begin(), vectorScalarProdResultDouble.end(), vectorScalarProdResult.begin()));
    /* complex multiplication (int, reverse) */
    CPPUNIT_ASSERT_NO_THROW(IntM2 * vectorScalarProd);
    const auto vectorScalarProdResultIntR = IntM2 * vectorScalarProd;
    CPPUNIT_ASSERT(std::equal(vectorScalarProdResultIntR.begin(), vectorScalarProdResultIntR.end(), vectorScalarProdResult.begin()));
    /* complex multiplication (double, reverse) */
    CPPUNIT_ASSERT_NO_THROW(DoubleM2 * vectorScalarProd);
    const auto vectorScalarProdResultDoubleR = DoubleM2 * vectorScalarProd;
    CPPUNIT_ASSERT(std::equal(vectorScalarProdResultDoubleR.begin(), vectorScalarProdResultDoubleR.end(), vectorScalarProdResult.begin()));
  }

  /**
   * @brief tests:
   *
   * std::vector<VectorElementType>
   * operator/(ScalarType lhs,
   *           const std::vector<VectorElementType> & rhs);
   *
   * std::vector<VectorElementType>
   * operator/(const std::vector<VectorElementType> & lhs,
   *           ScalarType rhs);
   */
  void
  testVectorScalarDivision()
  {
    /* dividing an empty vector with a non-zero scalar gives us an empty vector: { } / 1 = { } */
    /* scalar is an int */
    CPPUNIT_ASSERT_NO_THROW(emptyVector / Int1);
    const auto emptyVectorDivInt = emptyVector / Int1;
    CPPUNIT_ASSERT(std::equal(emptyVectorDivInt.begin(), emptyVectorDivInt.end(), emptyVector.begin()));
    /* scalar is a double */
    CPPUNIT_ASSERT_NO_THROW(emptyVector / Double1);
    const auto emptyVectorDivDouble = emptyVector / Double1;
    CPPUNIT_ASSERT(std::equal(emptyVectorDivDouble.begin(), emptyVectorDivDouble.end(), emptyVector.begin()));
    /* dividing any scalar with an empty vector gives us an empty vector: 1 / { } = { } */
    /* scalar is an int */
    CPPUNIT_ASSERT_NO_THROW(Int1 / emptyVector);
    const auto IntDivEmptyVector = Int1 / emptyVector;
    CPPUNIT_ASSERT(std::equal(IntDivEmptyVector.begin(), IntDivEmptyVector.end(), emptyVector.begin()));
    /* scalar is a double */
    CPPUNIT_ASSERT_NO_THROW(Double1 / emptyVector);
    const auto DoubleDivEmptyVector = Double1 / emptyVector;
    CPPUNIT_ASSERT(std::equal(DoubleDivEmptyVector.begin(), DoubleDivEmptyVector.end(), emptyVector.begin()));
    /* dividing any vector with a zero scalar throws an exception */
    /* vector is empty, scalar is an int */
    CPPUNIT_ASSERT_THROW(emptyVector / Int0, tthMEMexception);
    /* vector is empty, scalar is a double */
    CPPUNIT_ASSERT_THROW(emptyVector / Double0, tthMEMexception);
    /* vector is non-empty. scalar is an int */
    CPPUNIT_ASSERT_THROW(vectorWElem0 / Int0, tthMEMexception);
    /* vector is non-empty. scalar is an double */
    CPPUNIT_ASSERT_THROW(vectorWElem0 / Double0, tthMEMexception);
    /* dividing with a vector that contains a zero element throws an exception, regardless of the nature of scalar */
    /* scalar is an int */
    CPPUNIT_ASSERT_THROW(Int0 / vectorWElem0, tthMEMexception);
    /* scalar is an double */
    CPPUNIT_ASSERT_THROW(Double0 / vectorWElem0, tthMEMexception);
    /* complex division (int + reverse) */
    CPPUNIT_ASSERT_NO_THROW(vectorScalarDiv / IntM2);
    const auto vectorScalarDivResultInt = vectorScalarDiv / IntM2;
    CPPUNIT_ASSERT(std::equal(vectorScalarDivResultInt.begin(), vectorScalarDivResultInt.end(), vectorScalarDivResult.begin()));
    CPPUNIT_ASSERT_THROW(IntM2 / vectorScalarDiv, tthMEMexception);
    /* complex division (double + reverse) */
    CPPUNIT_ASSERT_NO_THROW(vectorScalarDiv / DoubleM2);
    const auto vectorScalarDivResultDouble = vectorScalarDiv / DoubleM2;
    CPPUNIT_ASSERT(std::equal(vectorScalarDivResultDouble.begin(), vectorScalarDivResultDouble.end(), vectorScalarDivResult.begin()));
    CPPUNIT_ASSERT_THROW(DoubleM2 / vectorScalarDiv, tthMEMexception);
    /* complex division 2 (int) */
    CPPUNIT_ASSERT_NO_THROW(vectorScalarDiv2 / IntM2);
    const auto vectorScalarDiv2ResultInt = vectorScalarDiv2 / IntM2;
    CPPUNIT_ASSERT(std::equal(vectorScalarDiv2ResultInt.begin(), vectorScalarDiv2ResultInt.end(), vectorScalarDiv2Result.begin()));
    /* complex division 2 (double) */
    CPPUNIT_ASSERT_NO_THROW(vectorScalarDiv2 / DoubleM2);
    const auto vectorScalarDiv2ResultDouble = vectorScalarDiv2 / DoubleM2;
    CPPUNIT_ASSERT(std::equal(vectorScalarDiv2ResultDouble.begin(), vectorScalarDiv2ResultDouble.end(), vectorScalarDiv2Result.begin()));
    /* complex division 2 (int, reverse) */
    CPPUNIT_ASSERT_NO_THROW(IntM2 / vectorScalarDiv2);
    const auto vectorScalarDiv2ResultIntR = IntM2 / vectorScalarDiv2;
    CPPUNIT_ASSERT(std::equal(vectorScalarDiv2ResultIntR.begin(), vectorScalarDiv2ResultIntR.end(), vectorScalarDiv2ResultR.begin()));
    /* complex division 2 (double, reverse) */
    CPPUNIT_ASSERT_NO_THROW(DoubleM2 / vectorScalarDiv2);
    const auto vectorScalarDiv2ResultDoubleR = DoubleM2 / vectorScalarDiv2;
    CPPUNIT_ASSERT(std::equal(vectorScalarDiv2ResultDoubleR.begin(), vectorScalarDiv2ResultDoubleR.end(), vectorScalarDiv2ResultR.begin()));
  }

  /**
   * @brief tests:
   *
   * std::vector<VectorElementType> &
   * operator+=(std::vector<VectorElementType> & lhs,
   *            ScalarType rhs);
   */
  void
  testVectorScalarAdditionAssignment()
  {
    /* adding to an empty vector does nothing */
    /* scalar is an int */
    {
      auto emptyVectorEmptyInt = emptyVector;
      CPPUNIT_ASSERT_NO_THROW(emptyVectorEmptyInt += Int0);
    }
    auto emptyVectorEmptyInt = emptyVector;
    emptyVectorEmptyInt += Int0;
    CPPUNIT_ASSERT(std::equal(emptyVectorEmptyInt.begin(), emptyVectorEmptyInt.end(), emptyVector.begin()));
    /* scalar is a double */
    {
      auto emptyVectorEmptyDouble = emptyVector;
      CPPUNIT_ASSERT_NO_THROW(emptyVectorEmptyDouble += Double0);
    }
    auto emptyVectorEmptyDouble = emptyVector;
    emptyVectorEmptyDouble += Double0;
    CPPUNIT_ASSERT(std::equal(emptyVectorEmptyDouble.begin(), emptyVectorEmptyDouble.end(), emptyVector.begin()));
    /* simple addition: { 0 } += 1 => { 1 } */
    /* scalar is an int */
    {
      auto zeroVectorToNonZeroInt = vectorWElem0;
      CPPUNIT_ASSERT_NO_THROW(zeroVectorToNonZeroInt += Int1);
    }
    auto zeroVectorToNonZeroInt = vectorWElem0;
    zeroVectorToNonZeroInt += Int1;
    CPPUNIT_ASSERT(std::equal(zeroVectorToNonZeroInt.begin(), zeroVectorToNonZeroInt.end(), vectorWElem1.begin()));
    /* scalar is a double */
    {
      auto zeroVectorToNonZeroDouble = vectorWElem0;
      CPPUNIT_ASSERT_NO_THROW(zeroVectorToNonZeroDouble += Double1);
    }
    auto zeroVectorToNonZeroDouble = vectorWElem0;
    zeroVectorToNonZeroDouble += Double1;
    CPPUNIT_ASSERT(std::equal(zeroVectorToNonZeroDouble.begin(), zeroVectorToNonZeroDouble.end(), vectorWElem1.begin()));
  }

  /**
   * @brief tests:
   *
   * std::vector<VectorElementType> &
   * operator-=(std::vector<VectorElementType> & lhs,
   *            ScalarType rhs);
   */
  void
  testVectorScalarSubtractionAssignment()
  {
    /* subtracting from an empty vector does nothing */
    /* scalar is an int */
    {
      auto emptyVectorEmptyInt = emptyVector;
      CPPUNIT_ASSERT_NO_THROW(emptyVectorEmptyInt -= Int0);
    }
    auto emptyVectorEmptyInt = emptyVector;
    emptyVectorEmptyInt -= Int0;
    CPPUNIT_ASSERT(std::equal(emptyVectorEmptyInt.begin(), emptyVectorEmptyInt.end(), emptyVector.begin()));
    /* scalar is a double */
    {
      auto emptyVectorEmptyDouble = emptyVector;
      CPPUNIT_ASSERT_NO_THROW(emptyVectorEmptyDouble -= Double0);
    }
    auto emptyVectorEmptyDouble = emptyVector;
    emptyVectorEmptyDouble -= Double0;
    CPPUNIT_ASSERT(std::equal(emptyVectorEmptyDouble.begin(), emptyVectorEmptyDouble.end(), emptyVector.begin()));
    /* simple subtraction: { 0 } -= 1 => { -1 } */
    /* scalar is an int */
    {
      auto zeroVectorToNonZeroInt = vectorWElem0;
      CPPUNIT_ASSERT_NO_THROW(zeroVectorToNonZeroInt -= Int1);
    }
    auto zeroVectorToNonZeroInt = vectorWElem0;
    zeroVectorToNonZeroInt -= Int1;
    CPPUNIT_ASSERT(std::equal(zeroVectorToNonZeroInt.begin(), zeroVectorToNonZeroInt.end(), vectorWElemM1.begin()));
    /* scalar is a double */
    {
      auto zeroVectorToNonZeroDouble = vectorWElem0;
      CPPUNIT_ASSERT_NO_THROW(zeroVectorToNonZeroDouble -= Double1);
    }
    auto zeroVectorToNonZeroDouble = vectorWElem0;
    zeroVectorToNonZeroDouble -= Double1;
    CPPUNIT_ASSERT(std::equal(zeroVectorToNonZeroDouble.begin(), zeroVectorToNonZeroDouble.end(), vectorWElemM1.begin()));
  }

  /**
   * @brief tests:
   *
   * std::vector<VectorElementType> &
   * operator*=(std::vector<VectorElementType> & lhs,
   *            ScalarType rhs);
   */
  void
  testVectorScalarMultiplicationAssignment()
  {
    /* multiplying an empty vector does nothing */
    /* scalar is an int */
    {
      auto emptyVectorEmptyInt = emptyVector;
      CPPUNIT_ASSERT_NO_THROW(emptyVectorEmptyInt *= Int0);
    }
    auto emptyVectorEmptyInt = emptyVector;
    emptyVectorEmptyInt *= Int0;
    CPPUNIT_ASSERT(std::equal(emptyVectorEmptyInt.begin(), emptyVectorEmptyInt.end(), emptyVector.begin()));
    /* scalar is a double */
    {
      auto emptyVectorEmptyDouble = emptyVector;
      CPPUNIT_ASSERT_NO_THROW(emptyVectorEmptyDouble *= Double0);
    }
    auto emptyVectorEmptyDouble = emptyVector;
    emptyVectorEmptyDouble *= Double0;
    CPPUNIT_ASSERT(std::equal(emptyVectorEmptyDouble.begin(), emptyVectorEmptyDouble.end(), emptyVector.begin()));
    /* complex multiplication & assignment */
    /* scalar is an int */
    {
      auto complexVectorInt = vectorScalarProd;
      CPPUNIT_ASSERT_NO_THROW(complexVectorInt *= IntM2);
    }
    auto complexVectorInt = vectorScalarProd;
    complexVectorInt *= IntM2;
    CPPUNIT_ASSERT(std::equal(complexVectorInt.begin(), complexVectorInt.end(), vectorScalarProdResult.begin()));
    /* scalar is a double */
    {
      auto complexVectorDouble = vectorScalarProd;
      CPPUNIT_ASSERT_NO_THROW(complexVectorDouble *= DoubleM2);
    }
    auto complexVectorDouble = vectorScalarProd;
    complexVectorDouble *= DoubleM2;
    CPPUNIT_ASSERT(std::equal(complexVectorDouble.begin(), complexVectorDouble.end(), vectorScalarProdResult.begin()));
  }

  /**
   * @brief tests:
   *
   * std::vector<VectorElementType> &
   * operator/=(std::vector<VectorElementType> & lhs,
   *            ScalarType rhs);
   */
  void
  testVectorScalarDivisionAssignment()
  {
    /* dividing an empty vector with non-zero scalar does nothing */
    /* scalar is an int */
    {
      auto emptyVectorEmptyNonZeroScalarInt = emptyVector;
      CPPUNIT_ASSERT_NO_THROW(emptyVectorEmptyNonZeroScalarInt /= Int1);
    }
    auto emptyVectorEmptyNonZeroScalarInt = emptyVector;
    emptyVectorEmptyNonZeroScalarInt /= Int1;
    CPPUNIT_ASSERT(std::equal(emptyVectorEmptyNonZeroScalarInt.begin(), emptyVectorEmptyNonZeroScalarInt.end(), emptyVector.begin()));
    /* scalar is a double */
    {
      auto emptyVectorEmptyNonZeroScalarDouble = emptyVector;
      CPPUNIT_ASSERT_NO_THROW(emptyVectorEmptyNonZeroScalarDouble /= Double1);
    }
    auto emptyVectorEmptyNonZeroScalarDouble = emptyVector;
    emptyVectorEmptyNonZeroScalarDouble /= Double1;
    CPPUNIT_ASSERT(std::equal(emptyVectorEmptyNonZeroScalarDouble.begin(), emptyVectorEmptyNonZeroScalarDouble.end(), emptyVector.begin()));
    /* dividing any vector with zero scalar throws an exception */
    /* scalar is an int, empty vector */
    auto emptyVectorEmptyZeroScalarInt = emptyVector;
    CPPUNIT_ASSERT_THROW(emptyVectorEmptyZeroScalarInt /= Int0, tthMEMexception);
    /* scalar is a double, empty vector */
    auto emptyVectorEmptyZeroScalarDouble = emptyVector;
    CPPUNIT_ASSERT_THROW(emptyVectorEmptyZeroScalarDouble /= Double0, tthMEMexception);
    /* scalar is an int, non-empty vector */
    auto nonEmptyVectorEmptyZeroScalarInt = vectorWElem0;
    CPPUNIT_ASSERT_THROW(nonEmptyVectorEmptyZeroScalarInt /= Int0, tthMEMexception);
    /* scalar is a double, non-empty vector */
    auto nonEmptyVectorEmptyZeroScalarDouble = vectorWElem0;
    CPPUNIT_ASSERT_THROW(nonEmptyVectorEmptyZeroScalarDouble /= Double0, tthMEMexception);
    /* complex division & assignment */
    /* scalar is a zero int */
    auto complexVectorScalarInt0 = vectorScalarDiv;
    CPPUNIT_ASSERT_THROW(complexVectorScalarInt0 /= Int0, tthMEMexception);
    /* scalar is a zero double */
    auto complexVectorScalarDouble0 = vectorScalarDiv;
    CPPUNIT_ASSERT_THROW(complexVectorScalarDouble0 /= Double0, tthMEMexception);
    /* scalar is a non-zero int */
    {
      auto complexVectorScalarInt = vectorScalarDiv;
      CPPUNIT_ASSERT_NO_THROW(complexVectorScalarInt /= IntM2);
    }
    auto complexVectorScalarInt = vectorScalarDiv;
    complexVectorScalarInt /= IntM2;
    CPPUNIT_ASSERT(std::equal(complexVectorScalarInt.begin(), complexVectorScalarInt.end(), vectorScalarDivResult.begin()));
    /* scalar is a non-zero double */
    {
      auto complexVectorScalarDouble = vectorScalarDiv;
      CPPUNIT_ASSERT_NO_THROW(complexVectorScalarDouble /= DoubleM2);
    }
    auto complexVectorScalarDouble = vectorScalarDiv;
    complexVectorScalarDouble /= DoubleM2;
    CPPUNIT_ASSERT(std::equal(complexVectorScalarDouble.begin(), complexVectorScalarDouble.end(), vectorScalarDivResult.begin()));
  }

  /**
   * @brief tests:
   *
   * std::vector<VectorElementType> &
   * operator+=(std::vector<VectorElementType> & lhs,
   *            const std::vector<VectorElementType> & rhs);
   */
  void
  testVectorAdditionAssignment()
  {
    /* adding an empty vector to an empty vector has no effect */
    auto emptyVectorToEmpty = emptyVector;
    CPPUNIT_ASSERT_NO_THROW(emptyVectorToEmpty += emptyVector);
    CPPUNIT_ASSERT(std::equal(emptyVectorToEmpty.begin(), emptyVectorToEmpty.end(), emptyVector.begin()));
    /* adding non-empty vector to an empty vector throws an exception */
    auto emptyVectorToNonEmpty = emptyVector;
    CPPUNIT_ASSERT_THROW(emptyVectorToNonEmpty += vectorWElem0, tthMEMexception);
    /* adding two vectors of the same size, simple case: { 0 } + { 1 } = { 1 } */
    auto vectorWElem0Copy = vectorWElem0;
    CPPUNIT_ASSERT_NO_THROW(vectorWElem0Copy += vectorWElem1);
    CPPUNIT_ASSERT(std::equal(vectorWElem0Copy.begin(), vectorWElem0Copy.end(), vectorWElem1.begin()));
    /* adding two vectors of the same size, complex case */
    auto vectorSubOperand0Copy = vectorAddOperand0;
    CPPUNIT_ASSERT_NO_THROW(vectorSubOperand0Copy += vectorAddOperand1);
    CPPUNIT_ASSERT(std::equal(vectorSubOperand0Copy.begin(), vectorSubOperand0Copy.end(), vectorAddResult.begin()));
  }

  /**
   * @brief tests:
   *
   * std::vector<VectorElementType> &
   * operator-=(std::vector<VectorElementType> & lhs,
   *            const std::vector<VectorElementType> & rhs);
   */
  void
  testVectorSubtractionAssignment()
  {
    /* subtracting empty vector from an empty vector has no effect */
    auto emptyVectorToEmpty = emptyVector;
    CPPUNIT_ASSERT_NO_THROW(emptyVectorToEmpty -= emptyVector);
    CPPUNIT_ASSERT(std::equal(emptyVectorToEmpty.begin(), emptyVectorToEmpty.end(), emptyVector.begin()));
    /* subtracting non-empty vector from an empty vector throws an exception */
    auto emptyVectorToNonEmpty = emptyVector;
    CPPUNIT_ASSERT_THROW(emptyVectorToNonEmpty -= vectorWElem0, tthMEMexception);
    /* subtracting an empty vector from a non-empty vector throws an exception */
    auto emptyVectorToNonEmpty2 = vectorWElem0;
    CPPUNIT_ASSERT_THROW(emptyVectorToNonEmpty2 -= emptyVector, tthMEMexception);
    /* subtracting two vectors of the same size, simple case: { 1 } - { 0 } = { 1 } */
    auto vectorWElem1Copy = vectorWElem1;
    CPPUNIT_ASSERT_NO_THROW(vectorWElem1Copy -= vectorWElem0);
    CPPUNIT_ASSERT(std::equal(vectorWElem1Copy.begin(), vectorWElem1Copy.end(), vectorWElem1.begin()));
    /* subtracting two vectors of the same size, complex case */
    auto vectorSubOperand0Copy = vectorSubOperand0;
    CPPUNIT_ASSERT_NO_THROW(vectorSubOperand0Copy -= vectorSubOperand1);
    CPPUNIT_ASSERT(std::equal(vectorSubOperand0Copy.begin(), vectorSubOperand0Copy.end(), vectorSubResult.begin()));
    /* subtracting two vectors of the same size, complex case (reverse) */
    auto vectorSubOperand1Copy = vectorSubOperand1;
    CPPUNIT_ASSERT_NO_THROW(vectorSubOperand1Copy -= vectorSubOperand0);
    CPPUNIT_ASSERT(std::equal(vectorSubOperand1Copy.begin(), vectorSubOperand1Copy.end(), vectorSubResultR.begin()));
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(Test_tthMEMvecFunctions);

