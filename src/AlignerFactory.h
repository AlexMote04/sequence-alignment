#pragma once

#include "ISequenceAligner.h"
#include "NeedlemanWunsch.h"
#include <string>
#include <memory>
#include <stdexcept>

// Factory function to create an aligner object.
inline std::unique_ptr<ISequenceAligner> create_aligner(
    const std::string &algorithm_name,
    int gap, int match, int mismatch)
{
  if (algorithm_name == "nw" || algorithm_name == "needleman-wunsch")
  {
    return std::make_unique<NeedlemanWunsch>(gap, match, mismatch);
  }

  // If no algorithm matches, throw an exception
  throw std::runtime_error("Unknown alignment algorithm: " + algorithm_name);
}