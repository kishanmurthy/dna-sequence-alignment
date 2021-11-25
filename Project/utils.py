from Bio import pairwise2

def check_bio_stats(X,Y, alpha):
    
    bio_dict = {}
    for k,v in alpha.items():
        for k1,v1 in v.items():
            bio_dict[(k,k1)] = -v1

    alignments = pairwise2.align.globalds(X, Y,bio_dict,-30,-30)
    return -1 * alignments[0].score