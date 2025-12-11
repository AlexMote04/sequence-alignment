#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>

// Maximum query length supported by the thread stack
// If you increase this too much (e.g. > 4000), you might hit Local Memory limits
#define MAX_QUERY_LEN 2048

// Serial Needleman-Wunsch (Score Only)
// This runs entirely inside ONE thread
// Uses O(min(N,M)) space optimization to avoid allocating a full matrix.
__device__ int nw_score_only(
    const char *query, int q_len,
    const char *target, int t_len,
    int match, int mismatch, int gap)
{
  // We allocate the DP column on the thread's local stack.
  // 'column' acts as the current vertical slice of the DP matrix.
  int column[MAX_QUERY_LEN + 1];

  // 1. Initialize the first column (Boundary condition: gap penalties)
  for (int i = 0; i <= q_len; ++i)
  {
    column[i] = i * gap;
  }

  // 2. Iterate through Target (Outer Loop = Columns of DP matrix)
  for (int j = 1; j <= t_len; ++j)
  {

    // 'diag' holds the value of the cell at (i-1, j-1)
    // At the start of a new column j, (0, j-1) is the top boundary
    int diag = column[0];

    // Update the top boundary for the current column j
    column[0] = j * gap;

    // Iterate through Query (Inner Loop = Rows of DP matrix)
    for (int i = 1; i <= q_len; ++i)
    {

      // Calculate costs
      int match_cost = (query[i - 1] == target[j - 1]) ? match : mismatch;

      // Retrieve previous values:
      // val_north = column[i];      // This is effectively (i, j-1) or "Left" in matrix terms
      // val_west  = column[i-1];    // This is effectively (i-1, j) or "Up" (newly computed)
      // val_nw    = diag;           // This is effectively (i-1, j-1) or "Diagonal"

      // Note: Naming depends on if you visualize Query on X or Y axis.
      // Standard: Query=Vertical(i), Target=Horizontal(j)

      int score_diag = diag + match_cost;
      int score_up = column[i] + gap;
      int score_left = column[i - 1] + gap;

      // Find Maximum
      int max_score = score_diag;
      if (score_up > max_score)
        max_score = score_up;
      if (score_left > max_score)
        max_score = score_left;

      // Save current column[i] as 'diag' for the next iteration (i+1)
      // BEFORE we overwrite it with the new max_score
      diag = column[i];

      // Update the array with the new score
      column[i] = max_score;
    }
  }

  // The result is the bottom-right corner of the matrix
  return column[q_len];
}

// ------------------------------------------------------------------
// Global Kernel: Database Search
// ------------------------------------------------------------------
__global__ void db_search_kernel(
    const char *query, int query_len,
    const char *db_data, const int *db_offsets, const int *db_lengths,
    int *results, int num_entries,
    int match, int mismatch, int gap)
{
  // 1. Compute global thread ID
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  // 2. Boundary Check
  if (tid >= num_entries)
    return;

  // 3. Locate the specific target sequence in flattened memory
  int target_idx = db_offsets[tid];
  int target_len = db_lengths[tid];
  const char *target_seq = &db_data[target_idx];

  // 4. Safety Check
  if (query_len > MAX_QUERY_LEN)
  {
    // Error flag if query exceeds stack limit
    results[tid] = -999999;
    return;
  }

  // 5. Compute Score
  results[tid] = nw_score_only(query, query_len, target_seq, target_len, match, mismatch, gap);
}

// ------------------------------------------------------------------
// C Wrapper (Called from C++)
// ------------------------------------------------------------------
extern "C" void launch_db_search_wrapper(
    const char *d_query, int query_len,
    const char *d_db_data, const int *d_db_offsets, const int *d_db_lengths,
    int *d_results, int num_entries,
    int match, int mismatch, int gap)
{
  // Define Grid Dimensions
  int threadsPerBlock = 256;
  int blocksPerGrid = (num_entries + threadsPerBlock - 1) / threadsPerBlock;

  // Launch Kernel
  db_search_kernel<<<blocksPerGrid, threadsPerBlock>>>(
      d_query, query_len,
      d_db_data, d_db_offsets, d_db_lengths,
      d_results, num_entries,
      match, mismatch, gap);

  // Check for launch errors (good practice)
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    printf("Kernel Launch Error: %s\n", cudaGetErrorString(err));
  }

  cudaDeviceSynchronize();
}