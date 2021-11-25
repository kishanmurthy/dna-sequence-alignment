from basic import SequenceAlignmentBasic
from efficient import SequenceAlignmentEfficient
from utils import compute_cost_bio_stat, compute_cost, generate_random_sequence

SEQUENCE_LENGTH = 100
TESTS = 100

ALPHA =  {'A': {'A':0,'C':110,'G':48,'T':94},
          'C': {'A':110,'C':0,'G':118,'T':48},
          'G':  {'A':48,'C':118,'G':0,'T':110},
          'T': {'A':94,'C':48,'G':110,'T':0}
          }

DELTA = 30

def generate_X_and_Y():
    return (generate_random_sequence(SEQUENCE_LENGTH), generate_random_sequence(SEQUENCE_LENGTH))


def run_basic(X,Y):
    sequence_alignment = SequenceAlignmentBasic(X,Y,ALPHA,DELTA)
    cost_dp = sequence_alignment.calculate_alignment_cost()
    X_align, Y_align = sequence_alignment.find_alignment()
    cost_alignment = compute_cost(X_align,Y_align,ALPHA,DELTA)
    return cost_dp , cost_alignment


def run_efficient(X,Y):
    sequence_alignment = SequenceAlignmentEfficient(X,Y,ALPHA,DELTA)
    cost_dp = sequence_alignment.calculate_alignment_cost()
    X_align, Y_align = sequence_alignment.find_alignment()
    cost_alignment = compute_cost(X_align,Y_align,ALPHA,DELTA)
    return cost_dp[-1] , cost_alignment


def test_sequence_alignment(X,Y):
    cost_dp_basic, cost_alignment_basic = run_basic(X,Y)
    cost_bio_stats = compute_cost_bio_stat(X,Y,ALPHA,DELTA)
    cost_dp_efficient, cost_alignment_efficient = run_efficient(X,Y)
    basic_pass = cost_dp_basic == cost_alignment_basic == int(cost_bio_stats)
    efficient_pass = cost_dp_efficient == cost_alignment_efficient == int(cost_bio_stats)

    return basic_pass, efficient_pass

def run_tests():
    failed_basic_cases = []
    failed_efficient_cases = []
    passed = 0
    for i in range(TESTS):
        X,Y = generate_X_and_Y()
        basic_pass, efficient_pass = test_sequence_alignment(X,Y)

        if basic_pass and efficient_pass:
            print(f"Test Passed: {i}")
            passed+=1
        else:
            if not basic_pass:
                print(f"Test Failed Basic: {i}")
                failed_basic_cases.append((X,Y))
            if not efficient_pass:
                print(f"Test Failed Efficient: {i}")
                failed_efficient_cases.append((X,Y))
    if passed == TESTS:
        print("All Tests Passed!")
    return failed_basic_cases, failed_efficient_cases


failed_basic_cases, failed_efficient_cases= run_tests()

print("Tests Completed!")