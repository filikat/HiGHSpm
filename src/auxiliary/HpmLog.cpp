#include "HpmLog.h"

namespace highspm {

const HighsLogOptions* Log::log_options_ = nullptr;

void Log::setOptions(const HighsLogOptions& log_options) {
  log_options_ = &log_options;
}

}  // namespace highspm
