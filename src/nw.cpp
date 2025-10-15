#include <iostream>
#include <vector>
#include <string>
#include <max>

constexpr int GAP_PENALTY = -1;

int main(int argc, char **argv)
{
  if (argc != 3)
  {
    std::cerr << "Usage: " << argv[0] << " <sequence 1> <sequence 2>" << std::endl;
  }

  std::string seq1 = argv[1];
  std::string seq2 = argv[2];

  int M = seq1.length();
  int N = seq2.length();

  // Needleman-Wunsch Algorithm
  std::vector<int> grid((M + 1) * (N + 1), 0);

  // Set first row
  for (int i = 1; i <= M; i++)
  {
    grid[i * (M + 1)] = grid[(i - 1) * (M + 1)] + GAP_PENALTY;
  }

  // Set first column
  for (int j = 1; j < N; j++)
  {
    grid[j] = grid[j - 1] + GAP_PENALTY;
  }

  // Declare re-declaring variable in loop
  int index;
  int match;
  int row_offset;

  // Fill matrix
  for (int i = 1; i < M + 1; i++)
  {
    row_offset = i * (M + 1);
    for (int j = 1; j < N + 1; j++)
    {
      // Match
      index = row_offset + j;
      match = seq1[i - 1] == seq2[j - 1] ? 1 : -1;
      grid[index] = std::max({grid[index - M - 2] + match, grid[index - 1] + GAP_PENALTY, grid[index - M - 1] + GAP_PENALTY});
    }
  }

  // Traceback
}