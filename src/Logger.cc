#include "tthAnalysis/tthMEM/interface/Logger.h"

#include "tthAnalysis/tthMEM/interface/Exception.h" // throw_line()

#include <chrono> // std::chrono::
#include <ctime> // std::time_t, std::tm
#include <iomanip> // std::setprecision(), std::setw()

#include <boost/format.hpp> // boost::format()
#include <boost/algorithm/string/predicate.hpp> // boost::iequals()

namespace tthMEM
{
  std::string
  getTimeStamp()
  {
    std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
    std::chrono::system_clock::duration tp = now.time_since_epoch();
    tp -= std::chrono::duration_cast<std::chrono::seconds>(tp);
    std::time_t tt = std::chrono::system_clock::to_time_t(now);
    std::tm t = *std::localtime(&tt);
    boost::format fmt("[%04u-%02u-%02u %02u:%02u:%02u.%03u]");
    fmt % (t.tm_year + 1900) % (t.tm_mon + 1) % t.tm_mday % t.tm_hour % t.tm_min
        % t.tm_sec % static_cast<unsigned>(tp / std::chrono::milliseconds(1));
    return fmt.str();
  }

  Logger::Logger(Logger::LogLevel logLevel)
    : holder_(new Holder(logLevel, *this))
  {}

  Logger::~Logger()
  {
    holder_.reset();
    holder_ = nullptr;
  }

  void
  Logger::setLogLevel(Logger::LogLevel logLevel)
  {
    logLevel_ = logLevel;
  }

  void
  Logger::setLogLevel(const std::string & logLevelString)
  {
    if(logLevelString == "") return; // using the default
    else if(boost::iequals(logLevelString, "error"))   logLevel_ = Logger::LogLevel::kError;
    else if(boost::iequals(logLevelString, "warning")) logLevel_ = Logger::LogLevel::kWarning;
    else if(boost::iequals(logLevelString, "fixme"))   logLevel_ = Logger::LogLevel::kFixme;
    else if(boost::iequals(logLevelString, "info"))    logLevel_ = Logger::LogLevel::kInfo;
    else if(boost::iequals(logLevelString, "debug"))   logLevel_ = Logger::LogLevel::kDebug;
    else if(boost::iequals(logLevelString, "verbose")) logLevel_ = Logger::LogLevel::kVerbose;
    else if(boost::iequals(logLevelString, "trace"))   logLevel_ = Logger::LogLevel::kTrace;
    else if(boost::iequals(logLevelString, "all"))     logLevel_ = Logger::LogLevel::kAll;
  }

  std::string
  Logger::getLogLevel()
  {
    if(logLevel_ == Logger::LogLevel::kError)   return "error";
    if(logLevel_ == Logger::LogLevel::kWarning) return "warning";
    if(logLevel_ == Logger::LogLevel::kFixme)   return "fixme";
    if(logLevel_ == Logger::LogLevel::kInfo)    return "info";
    if(logLevel_ == Logger::LogLevel::kDebug)   return "debug";
    if(logLevel_ == Logger::LogLevel::kVerbose) return "verbose";
    if(logLevel_ == Logger::LogLevel::kTrace)   return "trace";
    if(logLevel_ == Logger::LogLevel::kAll)     return "all";
    throw_line("Logger") << "Unspecified logging level";
  }

  void
  Logger::enableLogging(bool enableLogging)
  {
    enableLogging_ = enableLogging;
  }

  void
  Logger::enableTimeStamp()
  {
    enableTimeStamp_ = true;
  }

  void
  Logger::disableTimeStamp()
  {
    enableTimeStamp_ = false;
  }

  void
  Logger::enableTimeStamp(bool enableTimeStamp)
  {
    enableTimeStamp_ = enableTimeStamp;
  }

  void
  Logger::setFloatPrecision(unsigned floatPrecision)
  {
    floatPrecision_ = floatPrecision;
  }

  unsigned
  Logger::getFloatPrecision()
  {
    return floatPrecision_;
  }

  void
  Logger::flush()
  {
    os_.flush();
  }

  Logger::Holder::Holder(Logger::LogLevel logLevel,
                         Logger & logger)
    : logger_(logger)
    , enableLogging_(logger_.logLevel_ >= logLevel && logger_.enableLogging_)
  {
    if(enableLogging_)
    {
      ss_ << std::fixed << std::setprecision(logger_.floatPrecision_)
          << std::setw(10) << std::left << logger_.logLevelStrings_[logLevel];
      if(logger_.enableTimeStamp_) ss_ << tthMEM::getTimeStamp() << ' ';
    }
  }
  Logger::Holder::~Holder()
  {
    if(enableLogging_)
    {
      ss_ << '\n';
      logger_.os_ << ss_.str();
    }
  }

  std::vector<std::string> Logger::logLevelStrings_ =
  {
    "[error]", "[warning]", "[fixme]", "[info]", "[debug]", "[verbose]", "[trace]", "[all]"
  };

  Logger::LogLevel Logger::logLevel_ = Logger::LogLevel::kAll;
  bool Logger::enableLogging_ = true;
  bool Logger::enableTimeStamp_ = true;
  std::ostream & Logger::os_ = std::cout;
  unsigned Logger::floatPrecision_ = 3;

  Logger
  wrap(Logger::LogLevel logLevel)
  {
    return Logger(logLevel);
  }
}
