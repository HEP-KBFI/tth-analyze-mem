#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <FWCore/Utilities/interface/Exception.h> // cms::Exception

#if !defined(__FILENAME__)
#define __FILENAME__ (__builtin_strrchr(__FILE__, '/') ? __builtin_strrchr(__FILE__, '/') + 1 : __FILE__)
#endif

#define TTHEXCEPTION_ERR_CODE_DEFAULT                     1
#define TTHEXCEPTION_ERR_CODE_FILE_NOT_FOUND              2
#define TTHEXCEPTION_ERR_CODE_DIVISION_BY_ZERO            3
#define TTHEXCEPTION_ERR_CODE_NONZERO_CHARGE_SUM          4
#define TTHEXCEPTION_ERR_CODE_SLHA                        5
#define TTHEXCEPTION_ERR_CODE_MGME                        6
#define TTHEXCEPTION_ERR_CODE_VM                          7
#define TTHEXCEPTION_ERR_CODE_INSUFFICIENT_INTEGRAND      8
#define TTHEXCEPTION_ERR_CODE_HIGGS_WIDTH                 9
#define TTHEXCEPTION_ERR_CODE_UNDEFINED_LOGLEVEL          10
#define TTHEXCEPTION_ERR_CODE_UNDEFINED_FUNCTION          11
#define TTHEXCEPTION_ERR_CODE_UNDEFINED_INTEGRATION_MODE  12
#define TTHEXCEPTION_ERR_CODE_MISSING_JET_COMBINATIONS    13
#define TTHEXCEPTION_ERR_CODE_MISSING_INTEGRAND           14
#define TTHEXCEPTION_ERR_CODE_MISSING_PDF                 15
#define TTHEXCEPTION_ERR_CODE_MISSING_MG                  16
#define TTHEXCEPTION_ERR_CODE_MISSING_MISSING_PERMUTATION 17
#define TTHEXCEPTION_ERR_CODE_INVALID_RANGE               18
#define TTHEXCEPTION_ERR_CODE_INVALID_NOF_JETS            19
#define TTHEXCEPTION_ERR_CODE_INVALID_PERMUTATION         20
#define TTHEXCEPTION_ERR_CODE_INVALID_HELICITY            21
#define TTHEXCEPTION_ERR_CODE_INVALID_PARAMETERSET        22
#define TTHEXCEPTION_ERR_CODE_IDXW_ZERO_PERMUTATION       23
#define TTHEXCEPTION_ERR_CODE_IDXW_INVALID_PERMUTATION    24
#define TTHEXCEPTION_ERR_CODE_IDXW_INVALID_ARGUMENT       25
#define TTHEXCEPTION_ERR_CODE_IDXW_INVALID_RANGE          26
#define TTHEXCEPTION_ERR_CODE_LH_SIGNAL_ALREADY_ADDED     27
#define TTHEXCEPTION_ERR_CODE_LH_BKG_ALREADY_ADDED        28
#define TTHEXCEPTION_ERR_CODE_LH_INVALID_RANGE            29
#define TTHEXCEPTION_ERR_CODE_LH_OVERLAP                  30
#define TTHEXCEPTION_ERR_CODE_MXMC_INVALID_MODE           31
#define TTHEXCEPTION_ERR_CODE_MXMC_INVALID_ANNSUM         32
#define TTHEXCEPTION_ERR_CODE_MXMC_INVALID_ALPHA          33
#define TTHEXCEPTION_ERR_CODE_MXMC_INVALID_NOF_CHAINS     34
#define TTHEXCEPTION_ERR_CODE_MXMC_INVALID_NOF_BATCHES    35
#define TTHEXCEPTION_ERR_CODE_MXMC_INVALID_BATCH_DIV      36
#define TTHEXCEPTION_ERR_CODE_MXMC_INVALID_TEMP           37
#define TTHEXCEPTION_ERR_CODE_MXMC_INVALID_NU             38
#define TTHEXCEPTION_ERR_CODE_MXMC_MISSING_INTEGRAND      39
#define TTHEXCEPTION_ERR_CODE_MXMC_RUNTIME                40

/**
 * @brief The tthMEMexception class
 *
 * The class pretty-much replicates the functionality of cms::Exception class,
 * with an additional feature that the user can also supply an error code as an integer
 * to the constructor of the exception so that it could be caught in a try-catch block.
 * Weirdly enough, the cms::Exception class provides a returnCode() function
 * but it has no effect whatsoever since the class never stores or provides any
 * setters for the return code.
 */
class tthMEMexception
  : public cms::Exception
{
public:
  explicit tthMEMexception(const std::string & aCategory,
                           int errCode = TTHEXCEPTION_ERR_CODE_DEFAULT)
    : cms::Exception(aCategory)
    , errCode_(errCode)
  {}
  explicit tthMEMexception(char const *        aCategory,
                           int errCode = TTHEXCEPTION_ERR_CODE_DEFAULT)
    : cms::Exception(aCategory)
    , errCode_(errCode)
  {}

  tthMEMexception(std::string const& aCategory,
                  std::string const& message,
                  int errCode = TTHEXCEPTION_ERR_CODE_DEFAULT)
    : cms::Exception(aCategory, message)
    , errCode_(errCode)
  {}
  tthMEMexception(char const *        aCategory,
                  const std::string & message,
                  int errCode = TTHEXCEPTION_ERR_CODE_DEFAULT)
    : cms::Exception(aCategory, message)
    , errCode_(errCode)
  {}
  tthMEMexception(const std::string & aCategory,
                  char const *        message,
                  int errCode = TTHEXCEPTION_ERR_CODE_DEFAULT)
    : cms::Exception(aCategory, message)
    , errCode_(errCode)
  {}
  tthMEMexception(char const *        aCategory,
                  char const *        message,
                  int errCode = TTHEXCEPTION_ERR_CODE_DEFAULT)
    : cms::Exception(aCategory, message)
    , errCode_(errCode)
  {}

  tthMEMexception(const std::string & aCategory,
                  const std::string & message,
                  const Exception & another,
                  int errCode = TTHEXCEPTION_ERR_CODE_DEFAULT)
    : cms::Exception(aCategory, message, another)
    , errCode_(errCode)
  {}

  tthMEMexception(const Exception & other,
                  int errCode = TTHEXCEPTION_ERR_CODE_DEFAULT)
    : cms::Exception(other)
    , errCode_(errCode)
  {}

  int
  getErrCode() const
  {
    return errCode_;
  }

private:
  int errCode_;
};

// supply only category and throw message
#define throw_line(category) throw tthMEMexception(category) \
  << __FILENAME__ << ':' \
  << __FUNCTION__ << ':' \
  << __LINE__ << ": "

// supply category, throw message and the error code
// (cannot ,,overload'' macro definitions, hence new definition with a different name)
#define throw_line_ext(category, err_code) throw tthMEMexception(category, err_code) \
  << __FILENAME__ << ':' \
  << __FUNCTION__ << ':' \
  << __LINE__ << ": "

#endif // EXCEPTION_H
