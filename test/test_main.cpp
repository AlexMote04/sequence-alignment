#include <gtest/gtest.h>
#include <gmock/gmock.h> // UnorderedElementsAre
#include "NeedlemanWunsch.h"
#include <memory>

class NeedlemanWunschTest : public ::testing::Test
{
protected:
  std::unique_ptr<NeedlemanWunsch> aligner;

  void SetUp() override
  {
    aligner = std::make_unique<NeedlemanWunsch>(-2, 1, -1);
  }
};

// Test case for two empty strings
TEST_F(NeedlemanWunschTest, EmptyStrings)
{
  aligner->align("", "");

  // Score should be 0
  EXPECT_EQ(aligner->get_optimal_score(), 0);

  // Should be one alignment: two empty strings
  auto alignments = aligner->get_formatted_alignments();
  EXPECT_EQ(alignments.size(), 1);
  EXPECT_EQ(alignments[0].first, "");
  EXPECT_EQ(alignments[0].second, "");
}

// Test case for one empty string
TEST_F(NeedlemanWunschTest, OneEmptyString)
{
  // Scoring: gap = -2, match = 1, mismatch = -1
  aligner->align("A", "");

  // Aligning "A" with "" requires one gap. Score = -2.
  EXPECT_EQ(aligner->get_optimal_score(), -2);

  auto alignments = aligner->get_formatted_alignments();
  EXPECT_EQ(alignments.size(), 1);
  EXPECT_EQ(alignments[0].first, "A");
  EXPECT_EQ(alignments[0].second, "-");
}

// Test case for a perfect match
TEST_F(NeedlemanWunschTest, PerfectMatch)
{
  // Scoring: gap = -2, match = 1, mismatch = -1
  aligner->align("A", "A");

  // One match: Score = 1
  EXPECT_EQ(aligner->get_optimal_score(), 1);

  auto alignments = aligner->get_formatted_alignments();
  EXPECT_EQ(alignments.size(), 1);
  EXPECT_EQ(alignments[0].first, "A");
  EXPECT_EQ(alignments[0].second, "A");
}

// Test case for a simple mismatch
TEST_F(NeedlemanWunschTest, SimpleMismatch)
{
  // Scoring: gap = -2, match = 1, mismatch = -1
  aligner->align("A", "G");

  // The best score is one mismatch (-1),
  // which is better than two gaps (-2 + -2 = -4).
  EXPECT_EQ(aligner->get_optimal_score(), -1);

  auto alignments = aligner->get_formatted_alignments();
  EXPECT_EQ(alignments.size(), 1);
  EXPECT_EQ(alignments[0].first, "A");
  EXPECT_EQ(alignments[0].second, "G");
}

TEST_F(NeedlemanWunschTest, MultipleOptimalAlignments)
{
  aligner->align("AA", "A");

  // The optimal score is -1
  EXPECT_EQ(aligner->get_optimal_score(), -1);

  // The expected alignments
  AlignmentPair align1 = {"AA", "A-"};
  AlignmentPair align2 = {"AA", "-A"};

  // Use gmock's matcher to check that the vector contains
  // both expected pairs, in any order.
  EXPECT_THAT(aligner->get_formatted_alignments(),
              ::testing::UnorderedElementsAre(align1, align2));
}