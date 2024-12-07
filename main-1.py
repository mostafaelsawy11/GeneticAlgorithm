import random
import sys
import numpy as np


ITERATIONS = 50
POP_SIZE = 200
PM = 0.1
ELTISM_RATIO = 0.9
pc = 0.7
dependancyFactor = 3


class Task_allocation:
    def __init__(self,tasks,Time,boundries ,cost):
        self.pop_size = POP_SIZE
        self.n = int(self.pop_size * ELTISM_RATIO)
        if self.n % 2 == 1:
            self.n -= 1
        self.nEltism = self.pop_size - self.n
        self.max_time = Time
        self.nTasks = tasks
        self.boundries = boundries
        self.array = []
        self.generation = []
        self.selection = []
        self.cost = cost
        self.scores=[]
        self.t = 0
        self.RUN()

    def init_population(self):
        i = 0
        while i < self.pop_size:
            for x in range(self.nTasks):
                r1 = round(random.uniform(self.boundries[x*2], self.boundries[x*2+1]), 3)
                self.array.append(r1)
            if self.calc_fit(self.array) == self.max_time:
                self.generation.append(self.array)
            else:
                self.generation.append(self.updateArray(self.array))
            self.array = []
            i+=1
            
    def updateArray(self,upArr):
        res = self.calc_fit(upArr)
        while(res > self.max_time):
            # print("res if it more ==> ", res)
            for i in range(self.nTasks):
                if upArr[i] > self.boundries[i*2]:
                    plus = res - self.max_time
                    delta = upArr[i] - self.boundries[i*2]
                    r = round(random.uniform(0,delta), 3)
                    res -= min(r,plus)
                    res = round(res,3)
                    upArr[i] -= min(r,plus)
                    upArr[i] = round(upArr[i],3)
            if res == self.max_time + 0.001:
                break
    
        while(res < self.max_time):
            # print("res if it less ==> ", res)
            for i in range(self.nTasks):
                if upArr[i] < self.boundries[i*2+1]:
                    minus= self.max_time - res;
                    delta = self.boundries[i*2+1] - upArr[i]
                    r = round(random.uniform(0,delta), 3)
                    res += min(r,minus)
                    res = round(res,3)
                    upArr[i] += min(r,minus)
                    upArr[i] = round(upArr[i],3)
            if res == self.max_time - 0.001:
                break
                
        return upArr
    
    def SelectTasks(self,r1,r2):
         return self.calc_real_fitness(self.generation[r1],self.generation[r2])
    
    def calculate_Scores(self,arr):
        sum1 = 0
        for i in range(len(arr)):
            sum1 += (arr[i]*self.cost[i])
        return sum1
       
    def calc_real_fitness(self,arr1,arr2):
        sum1 = 0 
        sum2 = 0 
        for i in range(len(arr1)):
            sum1 += (arr1[i]*self.cost[i])
        for i in range(len(arr2)):
            sum2 += (arr2[i]*self.cost[i])
        if sum1<sum2:
            return arr1
        else:
            return arr2
    
    def Two_Point_CrossOverTwo(self,arr1,arr2,r1,r2):
        return arr1[:r1] + arr2[r1:r2]+arr1[r2:] , arr2[:r1]+arr1[r1:r2]+arr2[r2:]
        
    def Two_Point_CrossOver(self,arr1,arr2):
        r1 = random.randint(1,self.nTasks//2)
        r2 = random.randint(self.nTasks//2+1,self.nTasks-1)
        r = random.random()
        if(r <= pc):
            off1,off2 = self.Two_Point_CrossOverTwo(arr1,arr2,r1,r2)
            return off1,off2
        else:
            return arr1,arr2
        
    def NonUniformMutation(self,chr):
        for x in range(len(chr)):
            r = random.random()
            if r <= PM:
                r1 = random.random()
                delta = 0
                if r1 <= 0.5:
                    delta = chr[x] - self.boundries[x*2]
                else:
                    delta =  self.boundries[x*2+1]-chr[x] 
                r2 = random.random()
                dd = (1 - self.t / ITERATIONS)**dependancyFactor
                y = delta *(1-(r2**dd))
                if r1 < 0.5:
                    chr[x] -= y
                else:
                    chr[x] += y
        return self.updateArray(chr)

    def create_generation(self):
        for i in range(len(self.generation)):
            self.scores.append(self.calculate_Scores(self.generation[i]))
        min_indices = np.argsort(self.scores)[:self.nEltism]
        for i in min_indices:
            self.selection.append(self.generation[i])
        self.generation = self.selection
        self.selection = []
        self.scores = []

    def calc_fit(self, chemical):
        sum = 0
        for c in chemical:
            sum += c
        return sum
       
    def RUN(self):
        self.init_population()
        for i in range(ITERATIONS):
            while(len(self.selection)<self.pop_size - self.nEltism):
                r1 = random.randint(0,len(self.generation)-1)
                r2 = random.randint(0,len(self.generation)-1)
                off1 = self.SelectTasks(r1,r2)
                r1 = random.randint(0,len(self.generation)-1)
                r2 = random.randint(0,len(self.generation)-1)
                off2 = self.SelectTasks(r1,r2)
                off3,off4 = self.Two_Point_CrossOver(off1,off2)
                off1 = self.NonUniformMutation(off3)
                off2 = self.NonUniformMutation(off4)
                self.selection.append(off1)
                self.selection.append(off2)
            self.create_generation()
            self.t = 1
        result =  sys.maxsize
        var = []
        for i in range(len(self.generation)):
            sum = self.calculate_Scores(self.generation[i])
            if result > sum:
                result = sum
                var = self.generation[i]
        print("Chemical Proportions: ", var)
        print("Total Cost: ", result)

input_file = "input.txt"
with open(input_file, "r") as file:
    lines = file.readlines()

number_of_testcases = int(lines[0].strip())
index = 1

for x in range(number_of_testcases):
    number_of_tasks = int(lines[index].strip())
    index += 1
    limit = float(lines[index].strip())
    index += 1

    bound = []
    for i in range(number_of_tasks * 2):
        bound.append(float(lines[index].strip()))
        index += 1

    co = []
    for i in range(number_of_tasks):
        co.append(float(lines[index].strip()))
        index += 1

    print("Dataset ", x + 1)

    Task_allocation(
        tasks = number_of_tasks,
        Time = limit,
        boundries = bound,
        cost = co
    )
    print("\n\n----------------------------\n\n")