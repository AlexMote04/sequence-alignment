__global__ void needleman_wunsch_block(char *A, char *B, int lenA, int lenB, int *score_matrix)
{
  // 1. Allocate Shared Memory
  extern __shared__ int smem[];

  char *smem_chars = (char *)smem;

  char *s_A = smem_chars;        // Starts at byte 0
  char *s_B = smem_chars + lenA; // Starts at byte lenA

  // Ensure s_diag starts on a 4-byte boundary
  int vars_size = lenA + lenB;
  int aligned_offset = (vars_size + 3) / 4 * 4; // Round up to nearest multiple of 4

  // Note: We aren't using s_diag in this logic yet, but this is how you safely declare it
  int *s_diag = (int *)(smem_chars + aligned_offset);

  // 2. Cooperative Load (Coalesced)
  int tid = threadIdx.x;

  // Boundary check: Ensure we don't read past the real input arrays
  if (tid < lenA)
    s_A[tid] = A[tid];
  if (tid < lenB)
    s_B[tid] = B[tid];

  __syncthreads();

  // 3. Loop over diagonals
  // Max diagonals = lenA + lenB - 1
  for (int k = 0; k < lenA + lenB - 1; ++k)
  {
    // Mapping 1D threads to 2D diagonal slice
    int i = tid + 1;
    int j = k - i + 2;

    // 4. Check bounds
    // Valid cells must be within the matrix (1-based indices)
    bool active = (i > 0 && i <= lenA && j > 0 && j <= lenB);

    if (active)
    {
      // A. Calculate global memory indices
      int width = lenB + 1;
      int idx = i * width + j;

      // B. Read neighbors from Global Memory
      int val_north = score_matrix[(i - 1) * width + j];
      int val_west = score_matrix[i * width + (j - 1)];
      int val_nw = score_matrix[(i - 1) * width + (j - 1)];

      // C. Determine Match/Mismatch
      // Use Shared Memory for characters (Fast)
      // Note: sequence arrays are 0-indexed, matrix is 1-indexed
      char char_A = s_A[i - 1];
      char char_B = s_B[j - 1];

      int match_score = (char_A == char_B) ? 1 : -1;
      int gap_penalty = -2;

      // D. Calculate scores
      int score_diag = val_nw + match_score;
      int score_up = val_north + gap_penalty;
      int score_left = val_west + gap_penalty;

      // E. Find Max
      int max_score = score_diag;
      if (score_up > max_score)
        max_score = score_up;
      if (score_left > max_score)
        max_score = score_left;

      // F. Write result
      score_matrix[idx] = max_score;
    }

    // 7. Barrier Synchronization
    // Wait for diagonal k to finish before starting k+1
    __syncthreads();
  }
}