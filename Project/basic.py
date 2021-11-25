import tracemalloc
import time

def read_and_generate_strings(input_file):
    X=''
    Y=''
    with open(f'{input_file}') as file:
        lines = file.readlines()
        X = lines[0].strip()
        X_len = len(X)
        #print('X:',X_len)
        Y = ''
        flag,j,k = 0,0,0

        for i in range(1,len(lines)):
            val = lines[i].strip()
            try:
                val = int(val)

                if flag == 0:
                    j += 1
                    X = X[:val+1] + X + X[val+1:]
                else:
                    k += 1
                    Y = Y[:val+1] + Y + Y[val+1:]
            except:
                Y = val
                Y_len = len(Y)
                #print('Y:',len(Y))
                flag = 1
        
        #print(X,len(X),j,X_len,(2**j)*X_len)
        #print(Y,len(Y),k,Y_len,(2**k)*Y_len)

    return X,Y


class SequenceAlignmentBasic():
    
    def __init__(self,X,Y, alpha, delta):
        self.X = X
        self.Y = Y
        self.alpha = alpha
        self.delta = delta
        self.len_x = len(self.X)
        self.len_y = len(self.Y)
        self.dp = [[0 for j in range(self.len_y+1)] for i in range(self.len_x+1)]
        
    def calculate_alignment_cost(self):       

        for i in range(self.len_x+1):
            self.dp[i][0] = i*delta   

        for i in range(self.len_y+1):
            self.dp[0][i] = i*delta

        for i in range(1,self.len_x+1):
            for j in range(1,self.len_y+1):

                matching = self.dp[i-1][j-1] + self.alpha[self.X[i-1]][self.Y[j-1]]
                mismatch_x = self.dp[i-1][j] + self.delta
                mismatch_y = self.dp[i][j-1] + self.delta

                self.dp[i][j] = min(matching,mismatch_x,mismatch_y)

        return self.dp[self.len_x-1][self.len_y-1]



    def find_alignment(self):
        X_aligned = ""
        Y_aligned = ""
        i,j = self.len_x-1, self.len_y-1
        while i >= 0 and j >= 0:
            
            matching = self.dp[i][j] + self.alpha[self.X[i]][self.Y[j]]
            mismatch_x = self.dp[i][j+1] + self.delta
            mismatch_y = self.dp[i+1][j] + self.delta
            
            min_cost = min(matching,mismatch_x,mismatch_y)
            
            if min_cost==matching:
                X_aligned = self.X[i] + X_aligned
                Y_aligned = self.Y[j] + Y_aligned
                i-=1
                j-=1
            elif min_cost== mismatch_x:
                X_aligned = self.X[i] + X_aligned
                Y_aligned = '_' + Y_aligned
                i-=1
            elif min_cost== mismatch_y:
                X_aligned = '_' + X_aligned
                Y_aligned = self.Y[j] + Y_aligned
                j-=1

        while i >= 0:
            X_aligned = self.X[i] + X_aligned
            Y_aligned = '_' + Y_aligned
            i-=1

        while j >= 0:
            X_aligned = '_' + X_aligned
            Y_aligned = self.Y[j] + Y_aligned
            j-=1
        
        return X_aligned, Y_aligned

tracemalloc.start()
start = time.time()

alpha =  {'A': {'A':0,'C':110,'G':48,'T':94},
          'C': {'A':110,'C':0,'G':118,'T':48},
          'G':  {'A':48,'C':118,'G':0,'T':110},
          'T': {'A':94,'C':48,'G':110,'T':0}
          }


delta = 30

X,Y =read_and_generate_strings("basecase/input1.txt")
sequence_alignment = SequenceAlignmentBasic(X,Y,alpha,delta)
cost = sequence_alignment.calculate_alignment_cost()
X_align, Y_align = sequence_alignment.find_alignment()

output_file = open("output_basic.txt", "w")
        

output_file.write('First 50 Elements X: ' + str(X_align[:50])  + ' Y: ' + str(Y_align[:50] )+'\n')
output_file.write('Last  50 Elements X: ' + str(X_align[-50:]) + ' Y: ' + str(Y_align[-50:])+'\n')


current, peak = tracemalloc.get_traced_memory()
output_file.write(f"Time Taken: {time.time()-start}"+'\n')
output_file.write(f"Current memory usage is {current / 10**3} KB; Peak was {peak / 10**3} KB"+'\n')

output_file.close()
tracemalloc.stop()