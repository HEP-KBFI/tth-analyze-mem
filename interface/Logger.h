#ifndef LOGGER_H
#define LOGGER_H

#include <string> // std::string
#include <iostream> // std::ostream, std::cout
#include <memory> // std::shared_ptr<>

#define LOGERR  tthMEM::wrap(std::cout, "error")
#define LOGWARN tthMEM::wrap(std::cout, "warning")
#define LOGINFO tthMEM::wrap(std::cout, "info")

namespace tthMEM
{
  class Logger
  {
  public:
    Logger(std::ostream & os,
           const std::string & level)
      : holder_(new Holder(os, level))
    {}
    template <typename T>
    friend std::ostream & operator<<(const Logger & l,
                                     const T & t)
    {
      return (l.holder_ -> os_) << t;
    }
  private:
    class Holder
    {
    public:
      Holder(std::ostream & os,
             const std::string level)
        : os_(os)
      {
        os_ << "[" << level << "] ";
      }
      ~Holder()
      {
        os_ << "\n";
      }

      std::ostream & os_;
    };

    mutable std::shared_ptr<Holder> holder_;
  };

  Logger
  wrap(std::ostream & os,
       const std::string & level);
}

#endif // LOGGER_H
