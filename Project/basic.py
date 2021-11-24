import numpy as np
import pandas as pd

def read_and_generate_strings(input_file):
    X=''
    Y=''
    with open(f'{input_file}') as file:
        lines = file.readlines()
        X = lines[0].strip()

        Y = ''
        flag = 0

        for i in range(1,len(lines)):
            val = lines[i].strip()
            try:
                val = int(val)

                if flag == 0:
                    X = X[:val+1] + X + X[val+1:]
                else:
                    Y = Y[:val+1] + Y + Y[val+1:]
            except:
                Y = val
                flag = 1
    return X,Y



alpha =  {'A': {'A':0,'C':110,'G':48,'T':94},
          'C': {'A':110,'C':0,'G':118,'T':48},
          'G':  {'A':48,'C':118,'G':0,'T':110},
          'T': {'A':94,'C':48,'G':110,'T':0}
          }


delta = 30



def calculate_alignment_cost(X,Y):
    len_x = len(X)
    len_y = len(Y)
    
    dp = [[0 for j in range(len_y+1)] for i in range(len_x+1)]
    
    for i in range(len_x+1):
        dp[i][0] = i*delta
        
    for i in range(len_y+1):
        dp[0][i] = i*delta
        
        
    for i in range(1,len_x+1):
        for j in range(1,len_y+1):
            
            matching = dp[i-1][j-1] + alpha[X[i-1]][Y[j-1]]
            
            mismatch_x = dp[i-1][j] + delta
            mismatch_y = dp[i][j-1] + delta
            
            dp[i][j] = min(matching,mismatch_x,mismatch_y)
    return dp[len_x-1][len_y-1], dp



def find_alignment(X,Y,dp):
    X_aligned = ""
    Y_aligned = ""
    i,j = len(X)-1, len(Y)-1
    while i >= 0 and j >= 0:
        
        min_cost = min(dp[i][j] + alpha[X[i]][Y[j]] ,dp[i][j+1] + delta ,dp[i+1][j] + delta)
        if min_cost==(dp[i][j]+ alpha[X[i]][Y[j]]):
            X_aligned = X[i] + X_aligned
            Y_aligned = Y[j] + Y_aligned
            i-=1
            j-=1
        elif min_cost== (dp[i][j+1] + delta):
            X_aligned = X[i] + X_aligned
            Y_aligned = '_' + Y_aligned
            i-=1
        elif min_cost== (dp[i+1][j] + delta):
            X_aligned = '_' + X_aligned
            Y_aligned = Y[j] + Y_aligned
            j-=1
    
    while i >= 0:
        X_aligned = X[i] + X_aligned
        Y_aligned = '_' + Y_aligned
        i-=1

    while j >= 0:
        X_aligned = '_' + X_aligned
        Y_aligned = Y[j] + Y_aligned
        j-=1
    return X_aligned, Y_aligned



X,Y =read_and_generate_strings("input1.txt")

cost, dp = calculate_alignment_cost(X,Y)
pd.DataFrame(np.array(dp)).to_csv("tmp.csv",index=False)
X_align, Y_align = find_alignment(X,Y,dp)

print(X_align)
print(Y_align)