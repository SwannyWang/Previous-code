import numpy as np
import sys


args = sys.argv

train = []
with open(args[1], 'r') as file:
  for line in file:
    l=line.split('\t')
    train.append(l)

train = train[1:]

test = []
with open(args[2], 'r') as file:
  for line in file:
    l=line.split('\t')
    test.append(l)

test = test[1:]


zero_num_train = 0
one_num_train = 0
hd_train = []
for row in train:
    if row[len(row)-1][0] == '1':
        one_num_train += 1
        hd_train.append(1)
    else:
        zero_num_train +=1
        hd_train.append(0)

hd_test = []
for line in test:
    if line[len(line)-1][0] == '1':
        hd_test.append(1)
    else:
        hd_test.append(0)

pred_train = []
pred_test = []
if zero_num_train > one_num_train:
    for i in range(len(train)):
        pred_train.append(0)
    for j in range(len(test)):
        pred_test.append(0)
else:
    for i in range(len(train)):
        pred_train.append(1)
    for j in range(len(test)):
        pred_test.append(1)

e_train = 0
for i in range(len(hd_train)):
    if hd_train[i] != pred_train[i]:
        e_train += 1
er_train = e_train/len(train)

e_test = 0
for i in range(len(hd_test)):
    if hd_test[i] != pred_test[i]:
        e_test += 1
er_test = e_test/len(test)

#output txt files
f= open(args[3],"w+")
for result in pred_train:
        f.write(str(result))
        f.write('\n')
f.close()

f= open(args[4],"w+")
for result in pred_test:
        f.write(str(result))
        f.write('\n')
f.close()

f= open(args[5],"w+")
f.write('error(train): ' + str(er_train))
f.write('\n')
f.write('error(test): ' + str(er_test))
f.close()

