#ifndef HIGHSPM_LOG_H
#define HIGHSPM_LOG_H

#include <sstream>

#include "auxiliary/IntConfig.h"
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

  static void print(std::stringstream& ss);

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

// Functions to print using streams, taken from IPX.
std::string format(double d, Int width, Int prec,
                   std::ios_base::fmtflags floatfield);
std::string format(Int i, Int width);
inline std::string sci(double d, Int width, Int prec) {
  return format(d, width, prec, std::ios_base::scientific);
}
inline std::string fix(double d, Int width, Int prec) {
  return format(d, width, prec, std::ios_base::fixed);
}

}  // namespace highspm

#endif