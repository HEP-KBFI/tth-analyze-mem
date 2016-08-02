#ifndef LOGGER_H
#define LOGGER_H

#include <string> // std::string
#include <iostream> // std::ostream, std::cout
#include <memory> // std::shared_ptr<>
#include <sstream> // std::stringstream
#include <vector> // std::vector<>
#include <iomanip> // std::setw(), std::left

#define LOGERR  tthMEM::wrap(tthMEM::Logger::LogLevel::kError)  <<"["<<__FUNCTION__<<":"<<__LINE__<<"] "
#define LOGWARN tthMEM::wrap(tthMEM::Logger::LogLevel::kWarning)<<"["<<__FUNCTION__<<":"<<__LINE__<<"] "
#define LOGINFO tthMEM::wrap(tthMEM::Logger::LogLevel::kInfo)   <<"["<<__FUNCTION__<<":"<<__LINE__<<"] "
#define LOGDBG  tthMEM::wrap(tthMEM::Logger::LogLevel::kDebug)  <<"["<<__FUNCTION__<<":"<<__LINE__<<"] "

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
      kInfo    = 2,
      kDebug   = 3
    };

    Logger(Logger::LogLevel logLevel);

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
    static std::ostream * os_;
    static std::vector<std::string> logLevelStrings_;
  };

  Logger
  wrap(Logger::LogLevel logLevel);
}

#endif // LOGGER_H
