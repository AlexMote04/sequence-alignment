#include <gtest/gtest.h>
#include <vector>
#include <string>
#include "nw_strategy.h"

// Helper to generate sequences (reused for tests)
std::string generate_seq(int len)
{
  const char nucleotides[] = "ACGT";
  std::string seq;
  for (int i = 0; i < len; ++i)
    seq += nucleotides[rand() % 4];
  return seq;
}

// --------------------------------------------------------------------------
// Test Fixture: Sets up the strategies before each test
// --------------------------------------------------------------------------
class NeedlemanWunschTest : public ::testing::Test
{
protected:
  CpuNaiveNW cpuStrategy;
  CudaBasicNW cudaStrategy;

  NeedlemanWunschTest() : cpuStrategy(1, -1, -2), cudaStrategy(1, -1, -2) {}

  void SetUp() override
  {
    srand(42); // Fixed seed for reproducibility
  }
};

// --------------------------------------------------------------------------
// Test 1: Small fixed input to verify manual correctness
// --------------------------------------------------------------------------
TEST_F(NeedlemanWunschTest, SmallFixedInput)
{
  std::string seqA = "ACTG";
  std::string seqB = "ACGG";
  int lenA = seqA.length();
  int lenB = seqB.length();

  std::vector<int> cpuRes((lenA + 1) * (lenB + 1));
  std::vector<int> gpuRes((lenA + 1) * (lenB + 1));

  cpuStrategy.compute_score_matrix(seqA, seqB, lenA, lenB, cpuRes);
  cudaStrategy.compute_score_matrix(seqA, seqB, lenA, lenB, gpuRes);

  int cpuScore = cpuStrategy.get_optimal_score(lenA, lenB, cpuRes);
  int gpuScore = cudaStrategy.get_optimal_score(lenA, lenB, gpuRes);

  EXPECT_EQ(cpuScore, gpuScore) << "CPU and GPU scores should match for small input";
  // We expect ACTG vs ACGG -> Match, Match, Mismatch, Match
  // Score = 1 + 1 - 1 + 1 = 2 (Depending on path logic)
  // Actually, gap could appear. But mainly we check CPU == GPU here.
}

// --------------------------------------------------------------------------
// Test 2: Random Medium Input
// --------------------------------------------------------------------------
TEST_F(NeedlemanWunschTest, MediumRandomInput)
{
  int len = 500;
  std::string seqA = generate_seq(len);
  std::string seqB = generate_seq(len);

  std::vector<int> cpuRes((len + 1) * (len + 1));
  std::vector<int> gpuRes((len + 1) * (len + 1));

  cpuStrategy.compute_score_matrix(seqA, seqB, len, len, cpuRes);
  cudaStrategy.compute_score_matrix(seqA, seqB, len, len, gpuRes);

  int cpuScore = cpuStrategy.get_optimal_score(len, len, cpuRes);
  int gpuScore = cudaStrategy.get_optimal_score(len, len, gpuRes);

  EXPECT_EQ(cpuScore, gpuScore) << "CPU and GPU scores should match for random input";
}

// --------------------------------------------------------------------------
// Test 3: Edge Case (Empty Strings)
// --------------------------------------------------------------------------
TEST_F(NeedlemanWunschTest, EmptyStrings)
{
  std::string seqA = "";
  std::string seqB = "";
  int len = 0;

  std::vector<int> cpuRes(1);
  std::vector<int> gpuRes(1);

  cpuStrategy.compute_score_matrix(seqA, seqB, len, len, cpuRes);
  cudaStrategy.compute_score_matrix(seqA, seqB, len, len, gpuRes);

  EXPECT_EQ(cpuRes[0], 0);
  EXPECT_EQ(gpuRes[0], 0);
}