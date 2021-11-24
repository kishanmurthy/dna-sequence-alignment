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

def check_end_case_and_add_gaps(x,y):
    if x=="" or y=="":
        if x=="":
            return True,"_" * len(y), y
        elif y=="":
            return True, x, "_" * len(x)
        else:
            return True, "",""
    elif len(x)==1 and len(y)==1:
        return True, x,y
    return False, x, y

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

class SequenceAlignmentEfficient():
    
    def __init__(self,X,Y, alpha, delta):
        self.X = X
        self.Y = Y
        self.alpha = alpha
        self.delta = delta
        
    def calculate_alignment_cost(self,X=None,Y=None):
        if X==None and Y==None:
            X = self.X
            Y = self.Y
        len_x = len(X)
        len_y = len(Y)

        dp = [0 for j in range(len_y+1)]

        for j in range(len_y+1):
            dp[j] = j*self.delta

        for i in range(1,len_x+1):
            
            previous = dp[0]
            dp[0] = i* self.delta
            
            for j in range(1,len_y+1):

                matching =  previous + self.alpha[X[i-1]][Y[j-1]]
                mismatch_x = dp[j] + delta
                mismatch_y = dp[j-1] + delta
                previous = dp[j]
                dp[j] = min(matching,mismatch_x,mismatch_y)
        return dp[-1]
    
    
    def divide_one_element(self, X,Y):
        min_pair_cost = float("inf")
        min_pair_y = 0
        y_len = len(Y)
        
        for i in range(y_len):
            pair_cost = self.calculate_alignment_cost(X[0],Y[i]) + (y_len-1) * self.delta
            if pair_cost < min_pair_cost:
                min_pair_cost = pair_cost
                min_pair_y = i
        
        if min_pair_y>0:
            return "", X, Y[:min_pair_y], Y[min_pair_y:]
        else:
            return X, "", Y[0], Y[1:]
        
    def divide(self,X,Y):
        if len(X)==1:
            return self.divide_one_element(X,Y)

        mid = len(X) // 2
        x1, x2 = X[:mid] , X[mid:]
        x2_reverse = x2[::-1]
        y_reverse = Y[::-1]
        
        result1 = []
        result2 = []
        
        for i in range(0,len(Y)+1):
            result1.append(self.calculate_alignment_cost(x1, Y[:i]))
        
        for i in range(0,len(Y)+1):
            result2.append(self.calculate_alignment_cost(x2_reverse, y_reverse[:i]))

        result2.reverse()

        min_divide_cost = float("inf")
        divide_split = 0
        
        for i in range(0,len(Y)+1):
            if (result1[i] + result2[i]) < min_divide_cost:
                min_divide_cost = result1[i] + result2[i]
                divide_split = i
                
        return x1,x2,Y[:divide_split],Y[divide_split:]
    
    def divide_and_conquer(self,X,Y):
        x1, x2, y1, y2 = self.divide(X,Y)
    
        complete1,x1,y1 = check_end_case_and_add_gaps(x1,y1)
        complete2,x2,y2 = check_end_case_and_add_gaps(x2,y2)

        if not complete1 and not complete2:
            X1,Y1 = self.divide_and_conquer(x1,y1)
            X2,Y2 = self.divide_and_conquer(x2,y2)
            return (X1 + X2, Y1+Y2)
        elif not complete1:
            X1,Y1 = self.divide_and_conquer(x1,y1)
            return (X1+x2, Y1 +y2)
        elif not complete2:
            X2,Y2 = self.divide_and_conquer(x2,y2)
            return (x1+X2, y1 +Y2)
        return (x1+x2,y1+y2)
    
    def find_alignment(self):
        complete,X,Y = check_end_case_and_add_gaps(self.X,self.Y)
        if complete:
            return (X,Y)
        else:
            return self.divide_and_conquer(X,Y)


alpha =  {'A': {'A':0,'C':110,'G':48,'T':94},
          'C': {'A':110,'C':0,'G':118,'T':48},
          'G':  {'A':48,'C':118,'G':0,'T':110},
          'T': {'A':94,'C':48,'G':110,'T':0}
          }

delta = 30

X,Y = read_and_generate_strings("input1.txt")

sequence_alignment = SequenceAlignmentEfficient(X,Y,alpha,delta)
cost = sequence_alignment.calculate_alignment_cost()
X_align, Y_align = sequence_alignment.find_alignment()


print(X_align)
print(Y_align)