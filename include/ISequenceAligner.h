#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <utility>

// Type alias for clarity
using AlignmentPair = std::pair<std::string, std::string>;

// Interface for all sequence alignment algorithms.
// Defines the common methods all alignment algorithms must implement
class ISequenceAligner
{
public:
  // Virtual destructor
  virtual ~ISequenceAligner() = default;

  // Performs the alignment on two sequences.
  virtual void align(const std::string &s1, const std::string &s2) = 0;

  // Gets the optimal alignment score.
  virtual int get_optimal_score() const = 0;

  // Gets the vector representing the optimal score matrix
  virtual std::vector<int> get_optimal_score_matrix() const = 0;

  // Gets the vector representing the traceback matrix
  virtual std::vector<unsigned char> get_traceback_matrix() const = 0;

  // Gets a vector of all optimal alignment pairs.
  virtual std::vector<AlignmentPair> get_formatted_alignments() const = 0;

  // Helper function to print all results.
  void print_alignments(std::ostream &out = std::cout) const
  {
    // Call the pure virtual methods
    int score = get_optimal_score();
    std::vector<AlignmentPair> alignments = get_formatted_alignments();

    out << "Optimal Score: " << score << "\n";
    out << "Found " << alignments.size() << " optimal alignment(s):\n";

    int count = 1;
    for (const auto &pair : alignments)
    {
      out << "\n--- Alignment " << count++ << " ---\n";
      out << "Seq1: " << pair.first << "\n";
      out << "Seq2: " << pair.second << "\n";
    }
  }
};