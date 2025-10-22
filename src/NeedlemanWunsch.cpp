#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <sstream>

class NeedlemanWunsch
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

  // Algorithm data (cleared on each run)
  std::string seq1;
  std::string seq2;
  int M, N;
  std::vector<int> score_matrix;
  std::vector<unsigned char> traceback_matrix;
  std::vector<std::vector<unsigned char>> alignments;

  // Helper: Calculates the 1D index from 2D coordinates
  int index_1D(int i, int j) const
  {
    // Use N (Columns) to calculate the jump between rows
    return i * N + j;
  }

  // Helper: Calculates match/mismatch score
  int substitution_score(int i, int j) const
  {
    // Note: i and j are 1-based indices here, so we use i-1 and j-1
    if (seq1[i - 1] == seq2[j - 1])
      return match_score;
    return mismatch_penalty;
  }

public:
  // Constructor just sets up scoring
  NeedlemanWunsch(int gap = -2, int match = 1, int mismatch = -1)
      : gap_penalty(gap), match_score(match), mismatch_penalty(mismatch) {}

  // Public method to run the alignment
  void align(const std::string &s1, const std::string &s2)
  {
    // Set sequences and dimensions
    seq1 = s1;
    seq2 = s2;
    M = s1.length() + 1;
    N = s2.length() + 1;

    // Resize and clear old data
    score_matrix.assign(M * N, 0);
    traceback_matrix.assign(M * N, 0);
    alignments.clear();

    // Run the algorithm
    initialize_matrices();
    fill_matrices();
    traceback({}, M - 1, N - 1);
  }

private:
  void initialize_matrices()
  {
    // Initialization of Row 0 (Gaps in Seq1/Rows)
    for (int j = 0; j < N; j++)
    {
      score_matrix[index_1D(0, j)] = j * gap_penalty;
      traceback_matrix[index_1D(0, j)] = LEFT;
    }
    // Initialization of Column 0 (Gaps in Seq2/Columns)
    for (int i = 1; i < M; i++)
    {
      score_matrix[index_1D(i, 0)] = i * gap_penalty;
      traceback_matrix[index_1D(i, 0)] = UP;
    }
  }

  void fill_matrices()
  {
    for (int i = 1; i < M; i++)
    {
      for (int j = 1; j < N; j++)
      {
        // Calculate predecessor scores
        int S_diag = score_matrix[index_1D(i - 1, j - 1)] + substitution_score(i, j);
        int S_up = score_matrix[index_1D(i - 1, j)] + gap_penalty;
        int S_left = score_matrix[index_1D(i, j - 1)] + gap_penalty;

        // Find max score and update score matrix
        int S_max = std::max({S_diag, S_up, S_left});
        score_matrix[index_1D(i, j)] = S_max;

        // Determine and store the traceback bitmask
        unsigned char pointer = 0;
        if (S_diag == S_max)
          pointer |= DIAG;
        if (S_up == S_max)
          pointer |= UP;
        if (S_left == S_max)
          pointer |= LEFT;

        traceback_matrix[index_1D(i, j)] = pointer;
      }
    }
  }

  void traceback(std::vector<unsigned char> &path, int i, int j)
  {
    if (i == 0 && j == 0)
    {
      alignments.push_back(path);
    }
    else
    {
      int pointer = traceback_matrix[index_1D(i, j)];
      if (pointer & DIAG)
      {
        path.push_back(DIAG);
        traceback(path, i - 1, j - 1);
        path.pop_back(); // Backtrack by removing the last element
      }
      if (pointer & UP)
      {
        path.push_back(UP);
        traceback(path, i - 1, j);
        path.pop_back();
      }
      if (pointer & LEFT)
      {
        path.push_back(LEFT);
        traceback(path, i, j - 1);
        path.pop_back();
      }
    }
  }

public:
  int get_optimal_score_val() const
  {
    return score_matrix[index_1D(M - 1, N - 1)];
  }

  std::vector<std::pair<std::string, std::string>> get_formatted_alignments() const
  {
    std::vector<std::pair<std::string, std::string>> alignment_pairs;

    // Iterate over each path found by the traceback
    for (const auto &alignment_path : alignments)
    {
      std::stringstream sequence1;
      std::stringstream sequence2;

      // Initialize i and j inside the loop
      int i = 0;
      int j = 0;

      // Reconstruct the alignment strings by reading the path in reverse
      for (auto it = alignment_path.rbegin(); it != alignment_path.rend(); ++it)
      {
        if (*it == DIAG)
        {
          sequence1 << seq1[i];
          sequence2 << seq2[j];
          i++;
          j++;
        }
        else if (*it == UP)
        {
          sequence1 << seq1[i];
          sequence2 << "-";
          i++;
        }
        else // LEFT
        {
          sequence1 << "-";
          sequence2 << seq2[j];
          j++;
        }
      }

      // Add the constructed pair to the vector
      alignment_pairs.push_back({sequence1.str(), sequence2.str()});
    }

    // Return the completed vector
    return alignment_pairs;
  }
};