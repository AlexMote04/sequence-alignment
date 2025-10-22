#pragma once

#include "ISequenceAligner.h"
#include <vector>
#include <string>

class NeedlemanWunsch : public ISequenceAligner
{
private:
  // Scoring parameters
  int gap_penalty;
  int match_score;
  int mismatch_penalty;

  // Bitmasks for pointer matrix
  static constexpr unsigned char DIAG = 1; // 001
  static constexpr unsigned char UP = 2;   // 010
  static constexpr unsigned char LEFT = 4; // 100

  // Algorithm data
  std::string seq1;
  std::string seq2;
  int M, N;
  std::vector<int> score_matrix;
  std::vector<unsigned char> traceback_matrix;
  std::vector<std::vector<unsigned char>> alignments;

  // --- Helper Function Declarations ---
  int index_1D(int i, int j) const;
  int substitution_score(int i, int j) const;
  void initialize_matrices();
  void fill_matrices();
  void traceback(std::vector<unsigned char> &path, int i, int j);

public:
  // Constructor
  NeedlemanWunsch(int gap = -2, int match = 1, int mismatch = -1);

  // --- Interface Method Declarations ---
  void align(const std::string &s1, const std::string &s2) override;
  int get_optimal_score() const override;
  std::vector<AlignmentPair> get_formatted_alignments() const override;
};