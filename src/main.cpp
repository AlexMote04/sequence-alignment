#include "AlignerFactory.h"
#include "ISequenceAligner.h"
#include <iostream>
#include <string>
#include <memory>

int main(int argc, char *argv[])
{
  // Scoring
  int GAP = -2;
  int MATCH = 1;
  int MISMATCH = -1;

  // Check that we have the correct number of arguments
  if (argc != 3)
  {
    std::cerr << "Usage: " << argv[0] << " <SEQUENCE1> <SEQUENCE2>" << std::endl;
    return 1;
  }

  // Get sequences from command-line arguments
  std::string s1 = argv[1];
  std::string s2 = argv[2];
  std::string algo_name = "nw";

  try
  {
    // Use the factory to get the correct algorithm object.
    // 'aligner' is a smart pointer to the interface, not a concrete class.
    std::unique_ptr<ISequenceAligner> aligner =
        create_aligner(algo_name, GAP, MATCH, MISMATCH);

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