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
#define INFIXMSG_S INFIXMSG << std::scientific

#define LOGERR  tthMEM::wrap(tthMEM::Logger::LogLevel::kError)   << INFIXMSG
#define LOGWARN tthMEM::wrap(tthMEM::Logger::LogLevel::kWarning) << INFIXMSG
#define LOGFIX  tthMEM::wrap(tthMEM::Logger::LogLevel::kFixme)   << INFIXMSG
#define LOGINFO tthMEM::wrap(tthMEM::Logger::LogLevel::kInfo)    << INFIXMSG
#define LOGDBG  tthMEM::wrap(tthMEM::Logger::LogLevel::kDebug)   << INFIXMSG
#define LOGVRB  tthMEM::wrap(tthMEM::Logger::LogLevel::kVerbose) << INFIXMSG
#define LOGTRC  tthMEM::wrap(tthMEM::Logger::LogLevel::kTrace)   << INFIXMSG
#define LOGALL  tthMEM::wrap(tthMEM::Logger::LogLevel::kAll)     << INFIXMSG

#define LOGERR_S  tthMEM::wrap(tthMEM::Logger::LogLevel::kError)   << INFIXMSG_S
#define LOGWARN_S tthMEM::wrap(tthMEM::Logger::LogLevel::kWarning) << INFIXMSG_S
#define LOGFIX_S  tthMEM::wrap(tthMEM::Logger::LogLevel::kFixme)   << INFIXMSG_S
#define LOGINFO_S tthMEM::wrap(tthMEM::Logger::LogLevel::kInfo)    << INFIXMSG_S
#define LOGDBG_S  tthMEM::wrap(tthMEM::Logger::LogLevel::kDebug)   << INFIXMSG_S
#define LOGVRB_S  tthMEM::wrap(tthMEM::Logger::LogLevel::kVerbose) << INFIXMSG_S
#define LOGTRC_S  tthMEM::wrap(tthMEM::Logger::LogLevel::kTrace)   << INFIXMSG_S
#define LOGALL_S  tthMEM::wrap(tthMEM::Logger::LogLevel::kAll)     << INFIXMSG_S

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
      bool enableLogging_;
      std::stringstream ss_;
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
