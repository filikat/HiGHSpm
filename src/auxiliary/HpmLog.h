#ifndef HIGHSPM_LOG_H
#define HIGHSPM_LOG_H

#include <string>

#include "io/HighsIO.h"

// Interface to Highs logging.
// Call Log::setOptions to set the HighsLogOptions.
// Call Log::printf for normal printing, same syntax as printf.
// Call Log::printw for warnings, Log::printe for errors.
// If log_options_ is null, nothing is printed.

namespace highspm {

class Log {
  static const HighsLogOptions* log_options_;

  // Private ctor and dtor
  Log();
  ~Log() = default;

 public:
  static void setOptions(const HighsLogOptions& log_options);

  // normal printing
  template <typename... Args>
  static void printf(const char* format, Args... args) {
    if (log_options_)
      highsLogUser(*log_options_, HighsLogType::kInfo, format, args...);
  }

  // print warnings
  template <typename... Args>
  static void printw(const char* format, Args... args) {
    if (log_options_)
      highsLogUser(*log_options_, HighsLogType::kWarning, format, args...);
  }

  // print errors
  template <typename... Args>
  static void printe(const char* format, Args... args) {
    if (log_options_)
      highsLogUser(*log_options_, HighsLogType::kError, format, args...);
  }
};

}  // namespace highspm

#endif