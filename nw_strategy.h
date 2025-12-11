#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <chrono>

// Include CUDA runtime for memory functions
#include <cuda_runtime.h>

// ==========================================
// 1. Base Class
// ==========================================
using AlignmentPair = std::pair<std::string, std::string>;

class NeedlemanWunschBase
{
protected:
  int match_score;
  int mismatch_score;
  int gap_penalty;

  int get_index(int i, int j, int width) const
  {
    return i * width + j;
  }

public:
  NeedlemanWunschBase(int match = 1, int mismatch = -1, int gap = -2)
      : match_score(match), mismatch_score(mismatch), gap_penalty(gap) {}

  virtual ~NeedlemanWunschBase() = default;

  // Standard Signature: Inputs (Strings, Lengths), then Output (Matrix)
  virtual void compute_score_matrix(const std::string &s1, const std::string &s2, int lenA, int lenB, std::vector<int> &scoreMatrix) = 0;

  virtual std::string get_name() const = 0;

  int get_optimal_score(int lenA, int lenB, const std::vector<int> &scoreMatrix) const
  {
    return scoreMatrix[get_index(lenA, lenB, lenB + 1)];
  }

  AlignmentPair run_traceback(const std::string &s1, const std::string &s2, int lenA, int lenB, const std::vector<int> &scoreMatrix)
  {
    std::string alignA = "";
    std::string alignB = "";
    int i = lenA;
    int j = lenB;
    int width = lenB + 1;

    while (i > 0 || j > 0)
    {
      int currentScore = scoreMatrix[get_index(i, j, width)];
      bool wentDiag = false;

      if (i > 0 && j > 0)
      {
        int diagScore = scoreMatrix[get_index(i - 1, j - 1, width)];
        int scoreMove = (s1[i - 1] == s2[j - 1]) ? match_score : mismatch_score;
        if (currentScore == diagScore + scoreMove)
        {
          alignA += s1[i - 1];
          alignB += s2[j - 1];
          i--;
          j--;
          wentDiag = true;
        }
      }

      if (!wentDiag)
      {
        if (i > 0)
        {
          int upScore = scoreMatrix[get_index(i - 1, j, width)];
          if (currentScore == upScore + gap_penalty)
          {
            alignA += s1[i - 1];
            alignB += '-';
            i--;
          }
          else
          {
            alignA += '-';
            alignB += s2[j - 1];
            j--;
          }
        }
        else
        {
          alignA += '-';
          alignB += s2[j - 1];
          j--;
        }
      }
    }
    std::reverse(alignA.begin(), alignA.end());
    std::reverse(alignB.begin(), alignB.end());
    return {alignA, alignB};
  }
};

// ==========================================
// 2. CPU Implementation
// ==========================================
class CpuNaiveNW : public NeedlemanWunschBase
{
public:
  using NeedlemanWunschBase::NeedlemanWunschBase; // Inherit constructor

  std::string get_name() const override { return "CPU Sequential"; }

  void compute_score_matrix(const std::string &s1, const std::string &s2, int lenA, int lenB, std::vector<int> &scoreMatrix) override
  {
    int width = lenB + 1;

    // Initialize Boundaries
    for (int i = 0; i <= lenA; i++)
      scoreMatrix[get_index(i, 0, width)] = i * gap_penalty;
    for (int j = 0; j <= lenB; j++)
      scoreMatrix[get_index(0, j, width)] = j * gap_penalty;

    // Main Loop
    for (int i = 1; i <= lenA; i++)
    {
      for (int j = 1; j <= lenB; j++)
      {
        int val_north = scoreMatrix[get_index(i - 1, j, width)];
        int val_west = scoreMatrix[get_index(i, j - 1, width)];
        int val_nw = scoreMatrix[get_index(i - 1, j - 1, width)];

        int match_cost = (s1[i - 1] == s2[j - 1]) ? match_score : mismatch_score;

        int score_diag = val_nw + match_cost;
        int score_up = val_north + gap_penalty;
        int score_left = val_west + gap_penalty;

        scoreMatrix[get_index(i, j, width)] = std::max({score_diag, score_up, score_left});
      }
    }
  }
};

// ==========================================
// 3. CUDA Implementation Class
// ==========================================
// Forward declaration of the C wrapper function defined in .cu file
extern "C" void launch_nw_kernel_wrapper(const char *d_A, const char *d_B, int *d_matrix, int lenA, int lenB, int match, int mismatch, int gap);

class CudaBasicNW : public NeedlemanWunschBase
{
public:
  using NeedlemanWunschBase::NeedlemanWunschBase;

  std::string get_name() const override { return "CUDA Basic Global Memory"; }

  void compute_score_matrix(const std::string &seqA, const std::string &seqB, int lenA, int lenB, std::vector<int> &scoreMatrix) override
  {
    int width = lenB + 1;
    size_t matrixSize = (lenA + 1) * (lenB + 1) * sizeof(int);

    // 1. Allocate Device Memory
    char *d_A, *d_B;
    int *d_matrix;
    cudaMalloc(&d_A, lenA * sizeof(char));
    cudaMalloc(&d_B, lenB * sizeof(char));
    cudaMalloc(&d_matrix, matrixSize);

    // 2. Copy Data (Host -> Device)
    cudaMemcpy(d_A, seqA.c_str(), lenA * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, seqB.c_str(), lenB * sizeof(char), cudaMemcpyHostToDevice);

    // 3. Launch Kernel
    launch_nw_kernel_wrapper(d_A, d_B, d_matrix, lenA, lenB, match_score, mismatch_score, gap_penalty);

    // 4. Copy Result (Device -> Host)
    cudaMemcpy(scoreMatrix.data(), d_matrix, matrixSize, cudaMemcpyDeviceToHost);

    // 5. Free Memory
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_matrix);
  }
};

// ==========================================
// 4. Benchmarker
// ==========================================
class Benchmarker
{
public:
  void run_comparison(std::vector<NeedlemanWunschBase *> strategies,
                      const std::string &seqA,
                      const std::string &seqB)
  {

    int lenA = seqA.length();
    int lenB = seqB.length();
    std::vector<int> resultMatrix((lenA + 1) * (lenB + 1));

    std::cout << "--- Benchmarking (Input sizes: " << lenA << "x" << lenB << ") ---\n";

    for (auto &strategy : strategies)
    {
      std::fill(resultMatrix.begin(), resultMatrix.end(), 0);
      auto start = std::chrono::high_resolution_clock::now();

      strategy->compute_score_matrix(seqA, seqB, lenA, lenB, resultMatrix);

      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::milli> elapsed = end - start;
      std::cout << strategy->get_name() << " Matrix Fill Time: " << elapsed.count() << " ms" << std::endl;
    }
  }
};