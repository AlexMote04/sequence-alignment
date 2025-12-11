#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <ctime>

#include "nw_strategy.h"

// Helper to generate random DNA sequences
std::string generate_random_sequence(int length)
{
  const char nucleotides[] = "ACGT";
  std::string seq;
  seq.reserve(length);
  for (int i = 0; i < length; ++i)
  {
    seq += nucleotides[rand() % 4];
  }
  return seq;
}

int main()
{
  srand(time(0));

  // 1. Setup Test Data
  // NOTE: Keep lenA <= 1024 for the basic kernel!
  int lenA = 1000;
  int lenB = 1000;

  std::string seqA = generate_random_sequence(lenA);
  std::string seqB = generate_random_sequence(lenB);

  std::cout << "Generated sequences of length " << lenA << " and " << lenB << ".\n";

  // 2. Instantiate Strategies
  // (You can pass custom scoring: match=1, mismatch=-1, gap=-2)
  CpuNaiveNW cpuStrategy(1, -1, -2);
  CudaBasicNW cudaStrategy(1, -1, -2);

  // 3. CORRECTNESS CHECK (Critical for Project)
  std::cout << "\n--- Verifying Correctness ---\n";
  std::vector<int> cpuResult((lenA + 1) * (lenB + 1));
  std::vector<int> gpuResult((lenA + 1) * (lenB + 1));

  cpuStrategy.compute_score_matrix(seqA, seqB, lenA, lenB, cpuResult);
  cudaStrategy.compute_score_matrix(seqA, seqB, lenA, lenB, gpuResult);

  int cpuScore = cpuStrategy.get_optimal_score(lenA, lenB, cpuResult);
  int gpuScore = cudaStrategy.get_optimal_score(lenA, lenB, gpuResult);

  if (cpuScore == gpuScore)
  {
    std::cout << "[PASS] Scores match! (" << cpuScore << ")\n";

    // Run traceback on one of them to see alignment
    auto alignment = cudaStrategy.run_traceback(seqA, seqB, lenA, lenB, gpuResult);
    std::cout << "Alignment Length: " << alignment.first.length() << "\n";
  }
  else
  {
    std::cout << "[FAIL] Score mismatch! CPU: " << cpuScore << ", GPU: " << gpuScore << "\n";
    return 1;
  }

  // 4. RUN BENCHMARK
  std::cout << "\n--- Starting Benchmark ---\n";
  Benchmarker benchmarker;

  std::vector<NeedlemanWunschBase *> strategies;
  strategies.push_back(&cpuStrategy);
  strategies.push_back(&cudaStrategy);

  // You can run multiple iterations or different sizes here
  benchmarker.run_comparison(strategies, seqA, seqB);

  return 0;
}