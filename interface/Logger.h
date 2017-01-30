#ifndef LOGGER_H
#define LOGGER_H

#include <string> // std::string
#include <iostream> // std::ostream
#include <memory> // std::shared_ptr<>
#include <sstream> // std::stringstream
#include <vector> // std::vector<>

#if !defined(__FILENAME__)
#define __FILENAME__ (__builtin_strrchr(__FILE__, '/') ? __builtin_strrchr(__FILE__, '/') + 1 : __FILE__)
#endif

#define INFIXMSG '[' << __FILENAME__ << "][" << __FUNCTION__ << ':' << __LINE__ << "] "

#ifdef DISABLE_LOGGING
#define LOG_ENABLER if(0)
#else
#define LOG_ENABLER if(1)
#endif

#define LOGERR  LOG_ENABLER tthMEM::wrap(tthMEM::Logger::LogLevel::kError)   << INFIXMSG
#define LOGWARN LOG_ENABLER tthMEM::wrap(tthMEM::Logger::LogLevel::kWarning) << INFIXMSG
#define LOGFIX  LOG_ENABLER tthMEM::wrap(tthMEM::Logger::LogLevel::kFixme)   << INFIXMSG
#define LOGINFO LOG_ENABLER tthMEM::wrap(tthMEM::Logger::LogLevel::kInfo)    << INFIXMSG
#define LOGDBG  LOG_ENABLER tthMEM::wrap(tthMEM::Logger::LogLevel::kDebug)   << INFIXMSG
#define LOGVRB  LOG_ENABLER tthMEM::wrap(tthMEM::Logger::LogLevel::kVerbose) << INFIXMSG
#define LOGTRC  LOG_ENABLER tthMEM::wrap(tthMEM::Logger::LogLevel::kTrace)   << INFIXMSG
#define LOGALL  LOG_ENABLER tthMEM::wrap(tthMEM::Logger::LogLevel::kAll)     << INFIXMSG

#define LOGERR_S  LOGERR  << std::scientific
#define LOGWARN_S LOGWARN << std::scientific
#define LOGFIX_S  LOGFIX  << std::scientific
#define LOGINFO_S LOGINFO << std::scientific
#define LOGDBG_S  LOGDBG  << std::scientific
#define LOGVRB_S  LOGVRB  << std::scientific
#define LOGTRC_S  LOGTRC  << std::scientific
#define LOGALL_S  LOGALL  << std::scientific

namespace tthMEM
{
  std::string
  getTimeStamp();

  class Logger
  {
  public:
    enum LogLevel
    {
      kError   = 0,
      kWarning = 1,
      kFixme   = 2,
      kInfo    = 3,
      kDebug   = 4,
      kVerbose = 5,
      kTrace   = 6,
      kAll     = 7
    };

    Logger(Logger::LogLevel logLevel);
    ~Logger();

    template <typename T>
    friend std::ostream &
    operator<<(const Logger & l,
               const T & t)
    {
      return (l.holder_ -> ss_) << t;
    }

    static void
    enableLogging(bool enableLogging);

    static void
    setLogLevel(LogLevel logLevel);

    static void
    setLogLevel(const std::string & logLevelString);

    static std::string
    getLogLevel();

    static void
    enableTimeStamp();

    static void
    disableTimeStamp();

    static void
    enableTimeStamp(bool enableTimeStamp);

    static void
    setFloatPrecision(unsigned floatPrecision);

    static unsigned
    getFloatPrecision();

    static void
    flush();

  private:
    class Holder
    {
    public:
      Holder(Logger::LogLevel logLevel,
             Logger & logger);
      ~Holder();

      Logger & logger_;
      std::stringstream ss_;
      bool enableLogging_;
    };

    mutable std::shared_ptr<Holder> holder_;

    static LogLevel logLevel_;
    static bool enableLogging_;
    static bool enableTimeStamp_;
    static std::ostream & os_;
    static unsigned floatPrecision_;
    static std::vector<std::string> logLevelStrings_;
  };

  Logger
  wrap(Logger::LogLevel logLevel);
}

#endif // LOGGER_H
