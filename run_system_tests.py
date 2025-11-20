"""
This module runs a series of system tests comparing alignment 
scores between a gold standard (Bio.Align) my own implementations of NeedlemanWunsch
"""

import subprocess
import random
import re
import sys
from Bio import Align

# --- Configuration ---
GAP_PENALTY = -2
MATCH_SCORE = 1
MISMATCH_PENALTY = -1

# Path to your compiled C++ program
EXECUTABLE_PATH = "./build/sequence_aligner"
NUM_TESTS = 100

# Use biopython to get the "correct" score
def get_oracle_score(seq1: str, seq2: str) -> int:
    """Aligns two sequences using biopython and returns the score."""
    # Configure the aligner with our exact parameters
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global' # Use 'global' for Needleman-Wunsch
    aligner.match_score = MATCH_SCORE
    aligner.mismatch_score = MISMATCH_PENALTY

    # Biopython handles gap penalties differently.
    # We set a single "open" and "extend" penalty.
    aligner.open_gap_score = GAP_PENALTY 
    aligner.extend_gap_score = GAP_PENALTY

    # Biopython only aligns non-empty strings
    if not seq1 or not seq2:
        return (len(seq1) + len(seq2)) * GAP_PENALTY

    # Calculate the alignment score
    # We just get the score, as comparing formatted alignments is complex
    score = aligner.score(seq1, seq2)
    return int(score)

def get_my_program_score(seq1: str, seq2: str) -> int:
    """Runs the C++ executable and parses its stdout for the score."""

    # Run the C++ program as a subprocess
    result = subprocess.run(
        [EXECUTABLE_PATH, seq1, seq2],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print("--- C++ Program Error ---")
        print(result.stderr)
        raise RuntimeError("Your C++ program crashed.")

    # Parse the stdout
    # We look for the line "Optimal Score: X"
    match = re.search(r"Optimal Score: (-?\d+)", result.stdout)
    if not match:
        print("--- C++ Program Output ---")
        print(result.stdout)
        raise ValueError("Could not parse the 'Optimal Score' from C++ output.")

    return int(match.group(1))

# Test Helper
def generate_random_sequence(max_len:int = 50) -> str:
    """Generates a random DNA sequence."""
    length = random.randint(0, max_len)
    return "".join(random.choice("ATGC") for _ in range(length))

# Main Test Runner
def run_tests():
    print(f"--- Running {NUM_TESTS} Random System Tests ---")
    for i in range(NUM_TESTS):
        # Generate two random sequences
        s1 = generate_random_sequence()
        s2 = generate_random_sequence()

        print(f"Test {i+1}/{NUM_TESTS}: Aligning '{s1}' and '{s2}'...", end="")

        try:
            # Get scores from both implementations
            oracle_score = get_oracle_score(s1, s2)
            my_score = get_my_program_score(s1, s2)

            # Compare
            assert my_score == oracle_score
            print(" PASSED")

        except AssertionError as e:
            print(" FAILED!")
            print(f"  Sequences: '{s1}', '{s2}'")
            print(f"  Oracle Score (Biopython): {oracle_score}")
            print(f"  My Score (C++):         {my_score}")
            print(f"  Error: {e}")
            sys.exit(1) # Stop on first failure

    print(f"\n--- All {NUM_TESTS} System Tests Passed! ---")

if __name__ == "__main__":
    run_tests()
