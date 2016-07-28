#ifndef LOGGER_H
#define LOGGER_H

#include <string> // std::string
#include <iostream> // std::ostream, std::cout
#include <memory> // std::shared_ptr<>
#include <sstream> // std::stringstream

#define LOGERR  tthMEM::wrap(std::cout, "error",   true)
#define LOGWARN tthMEM::wrap(std::cout, "warning", true)
#define LOGINFO tthMEM::wrap(std::cout, "info",    true)

namespace tthMEM
{
  class Logger
  {
  public:
    Logger(std::ostream & os,
           const std::string & level,
           bool enable)
      : holder_(new Holder(os, level, enable))
    {}
    template <typename T>
    friend std::ostream & operator<<(const Logger & l,
                                     const T & t)
    {
      return (l.holder_ -> ss_) << t;
    }
  private:
    class Holder
    {
    public:
      Holder(std::ostream & os,
             const std::string level,
             bool enable)
        : os_(os)
        , enable_(enable)
      {
        ss_ << "[" << level << "] ";
      }
      ~Holder()
      {
        ss_ << "\n";
        if(enable_) os_ << ss_.str();
      }

      std::ostream & os_;
      std::stringstream ss_;
      bool enable_;
    };

    mutable std::shared_ptr<Holder> holder_;
  };

  Logger
  wrap(std::ostream & os,
       const std::string & level,
       bool enable);
}

#endif // LOGGER_H
