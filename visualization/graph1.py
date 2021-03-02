import csv
import matplotlib.pyplot as plt
import numpy as np
import math

x1 = []
p1 = []
p2 = []
q1 = []
q2 = []
r1 = []
r2 = []
s1 = []
s2 = []



with open('../results/exp1.csv') as file:
    csv_reader = csv.reader(file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if (line_count==0):
            numOfRuns = float(row[9])
        else:
            x1.append(float(row[0]))
            p1.append(float(row[1]))
            p2.append(float(row[2]))
            q1.append(float(row[3]))
            q2.append(float(row[4]))
            r1.append(float(row[5]))
            r2.append(float(row[6]))
            s1.append(float(row[7]))
            s2.append(float(row[8]))
        line_count = line_count + 1

print("Num of runs is ",numOfRuns)

plt.errorbar(x1, p1, [x/math.sqrt(numOfRuns) for x in p2], label = 'CG')
plt.errorbar(x1, q1, [x/math.sqrt(numOfRuns) for x in q2], label = 'G')
plt.errorbar(x1, r1, [x/math.sqrt(numOfRuns) for x in r2], label = 'DPG')
plt.errorbar(x1, s1, [x/math.sqrt(numOfRuns) for x in s2], label = 'R')
plt.xlabel("Rank")
plt.ylabel("Average utility")
plt.legend()
plt.title("Utility vs Rank, 20 grid 80 spurious")
plt.savefig("Cardinality constraint")
plt.show()
