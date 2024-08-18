#ifndef SIGNAL_HANDLER_H
#define SIGNAL_HANDLER_H

#include <chrono>
#include <csignal>
#include <exception>
#include <iostream>
#include <thread>

class SignalException : public std::exception
{
public:
  SignalException(int signal) : signal_(signal)
  {
  }
  virtual const char *what() const noexcept override
  {
    return "Signal received";
  }
  int getSignal() const
  {
    return signal_;
  }

private:
  int signal_;
};

extern std::exception_ptr globalExceptionPtr;

void signalHandler(int signal);

#endif  // SIGNAL_HANDLER_H
