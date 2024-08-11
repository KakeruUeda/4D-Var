#include "SignalHandler.h"

std::exception_ptr globalExceptionPtr = nullptr;

void signalHandler(int signal)
{
  globalExceptionPtr = std::make_exception_ptr(SignalException(signal));
}