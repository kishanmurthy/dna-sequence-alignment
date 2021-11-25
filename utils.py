from Bio import pairwise2
import random

def compute_cost_bio_stat(X,Y, alpha,delta):
    
    bio_dict = {}
    for k,v in alpha.items():
        for k1,v1 in v.items():
            bio_dict[(k,k1)] = -v1

    alignments = pairwise2.align.globalds(X, Y,bio_dict,-delta,-delta)
    return -1 * alignments[0].score


def compute_cost(X_align, Y_align, alpha, delta):
    i = 0
    cost = 0
    while i<len(X_align):
        if X_align[i] != '_' and Y_align[i] != '_':
            cost+= alpha[X_align[i]][Y_align[i]]
        else:
            cost+= delta
        i+=1
    return cost


def generate_random_sequence(length):
    letters = 'ACGT'
    result_str = ''.join(random.choice(letters) for i in range(length))
    return result_str