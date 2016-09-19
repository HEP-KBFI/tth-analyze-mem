#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <FWCore/Utilities/interface/Exception.h> // cms::Exception

#if !defined(__FILENAME__)
#define __FILENAME__ (__builtin_strrchr(__FILE__, '/') ? __builtin_strrchr(__FILE__, '/') + 1 : __FILE__)
#endif

#define throw_line(category) throw cms::Exception(category) \
  << __FILENAME__ << ':' \
  << __FUNCTION__ << ':' \
  << __LINE__ << ": "

#endif // EXCEPTION_H
