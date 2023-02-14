import numpy as np
import sys
args = sys.argv

input = []
with open(args[1], 'r') as file:
  for line in file:
    l=line.split('\t')
    input.append(l)

input = input[1:]

zero_num = 0
one_num = 0
current = []
for row in input:
    if row[len(row)-1][0] == '1':
        one_num += 1
        current.append(1)
    else:
        zero_num +=1
        current.append(0)

pred = []
if zero_num > one_num:
    for i in range(len(input)):
        pred.append(0)
else:
    for i in range(len(input)):
        pred.append(1)

e_input = 0
for i in range(len(input)):
    if current[i] != pred[i]:
        e_input += 1
er_input = e_input/len(input)

p = []
p.append(zero_num/len(input))
p.append(one_num/len(input))
entro = -sum(p*np.log2(p))


f= open(args[2],"w+")
f.write("entropy: " + str(entro))
f.write('\n')
f.write("error: " + str(er_input))
f.close()