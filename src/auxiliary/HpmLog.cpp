#include "HpmLog.h"

namespace highspm {

const HighsLogOptions* Log::log_options_ = nullptr;

void Log::setOptions(const HighsLogOptions& log_options) {
  log_options_ = &log_options;
}

void Log::print(std::stringstream& ss) {
  if (log_options_)
    highsLogUser(*log_options_, HighsLogType::kInfo, "%s", ss.str().c_str());
}

std::string format(double d, Int width, Int prec,
                   std::ios_base::fmtflags floatfield) {
  std::ostringstream s;
  s.precision(prec);
  s.width(width);
  s.setf(floatfield, std::ios_base::floatfield);
  s << d;
  return s.str();
}

std::string format(Int i, int width) {
  std::ostringstream s;
  s.width(width);
  s << i;
  return s.str();
}

}  // namespace highspm
