#include "AlignerFactory.h"
#include "ISequenceAligner.h"
#include <iostream>
#include <string>
#include <memory>

int main()
{
  // Temp before implementing argv
  std::string s1 = "GACATACA";
  std::string s2 = "CAT";
  std::string algo_name = "nw";

  // Scoring
  int gap = -2;
  int match = 1;
  int mismatch = -1;

  try
  {
    // Use the factory to get the correct algorithm object.
    // 'aligner' is a smart pointer to the interface, not a concrete class.
    std::unique_ptr<ISequenceAligner> aligner =
        create_aligner(algo_name, gap, match, mismatch);

    // Perform the alignment
    aligner->align(s1, s2);

    // Print optimal alignments
    aligner->print_alignments();
  }
  catch (const std::exception &e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}