#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>

// ------------------------------------------------------------------
// Kernel 1: Boundary Initialization
// ------------------------------------------------------------------
__global__ void init_boundaries(int *score_matrix, int lenA, int lenB, int gap_penalty)
{
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  int width = lenB + 1;

  // Initialize Top Row (Row 0, Columns 0 to lenB)
  if (idx <= lenB)
  {
    score_matrix[idx] = idx * gap_penalty;
  }

  // Initialize Left Column (Column 0, Rows 0 to lenA)
  if (idx <= lenA)
  {
    score_matrix[idx * width] = idx * gap_penalty;
  }
}

// ------------------------------------------------------------------
// Kernel 2: Wavefront Computation
// ------------------------------------------------------------------
__global__ void needleman_wunsch_block(const char *A, const char *B, int lenA, int lenB, int match, int mismatch, int gap_penalty, int *score_matrix)
{
  extern __shared__ char smem_chars[];

  char *s_A = smem_chars;
  char *s_B = smem_chars + lenA;

  int tid = threadIdx.x;
  int width = lenB + 1;

  if (tid < lenA)
    s_A[tid] = A[tid];
  if (tid < lenB)
    s_B[tid] = B[tid];

  __syncthreads();

  for (int k = 0; k < lenA + lenB - 1; ++k)
  {
    int i = tid + 1;
    int j = k - tid + 1;

    bool active = (i > 0 && i <= lenA && j > 0 && j <= lenB);

    if (active)
    {
      int idx = i * width + j;

      int val_north = score_matrix[(i - 1) * width + j];
      int val_west = score_matrix[i * width + (j - 1)];
      int val_nw = score_matrix[(i - 1) * width + (j - 1)];

      char char_A = s_A[i - 1];
      char char_B = s_B[j - 1];

      int match_score = (char_A == char_B) ? match : mismatch;

      int score_diag = val_nw + match_score;
      int score_up = val_north + gap_penalty;
      int score_left = val_west + gap_penalty;

      int max_score = score_diag;
      if (score_up > max_score)
        max_score = score_up;
      if (score_left > max_score)
        max_score = score_left;

      score_matrix[idx] = max_score;
    }
    __syncthreads();
  }
}

// ------------------------------------------------------------------
// CRITICAL: The Wrapper Function to be called from C++
// ------------------------------------------------------------------
extern "C" void launch_nw_kernel_wrapper(const char *d_A, const char *d_B, int *d_matrix, int lenA, int lenB, int match, int mismatch, int gap)
{
  // Constraint check for this specific basic kernel
  if (lenA > 1024)
  {
    printf("Error: Basic kernel supports max lenA=1024\n");
    return;
  }

  // 1. Initialize Boundaries (New Step)
  // We launch enough threads to cover the longer dimension
  int max_dim = (lenA > lenB) ? lenA : lenB;
  int initThreads = 256;
  int initBlocks = (max_dim + initThreads) / initThreads; // Ceiling division

  init_boundaries<<<initBlocks, initThreads>>>(d_matrix, lenA, lenB, gap);

  // Ensure initialization is done before starting the wavefront
  cudaDeviceSynchronize();

  // 2. Run Wavefront
  int threadsPerBlock = lenA;
  int blocksPerGrid = 1;
  size_t sharedMemSize = (lenA + lenB) * sizeof(char);

  needleman_wunsch_block<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(
      d_A, d_B, lenA, lenB, match, mismatch, gap, d_matrix);

  cudaDeviceSynchronize();
}