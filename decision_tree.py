import numpy as np
import sys

class Node:
    '''
    Here is an arbitrary Node class that will form the basis of your decision
    tree. 
    Note:
        - the attributes provided are not exhaustive: you may add and remove
        attributes as needed, and you may allow the Node to take in initial
        arguments as well
        - you may add any methods to the Node class if desired 
    '''
    def __init__(self, data_new, depth, used_attr_index, v, left = None, right = None, attr = None):
        self.data_new = data_new
        self.depth = depth
        self.used_attr_index = used_attr_index
        self.vote = v
        self.left = left
        self.right = right
        self.attr = attr




args = sys.argv

train = []
with open(args[1], 'r') as file:
  for line in file:
    l=line.split('\t')
    train.append(l)

train = np.array(train)
for i in range(train.shape[0]):
    train[i,-1] = train[i,-1][0:-1]
train = train[1:,]
attributes = train[0,:]

test = []
with open(args[2], 'r') as file:
  for line in file:
    l=line.split('\t')
    test.append(l)

test = np.array(test)
for i in range(test.shape[0]):
    test[i,-1] = test[i,-1][0:-1]
test = test[1:,]

max_depth = int(args[3])

def entro_calculator(input):
    zero_num = 0
    one_num = 0
    p = []
    a = 0
    uni = sorted(np.unique(input))
    if len(uni) == 1:
        return a
    else:
        for y in input:
            if y[0] == uni[1]:
                one_num += 1
            else:
                zero_num += 1
        p.append(zero_num / len(input))
        p.append(one_num / len(input))
        for i in p:
            if i != 0:
                a += -i * np.log2(i)
    return a

def mi_calculator(predictor,response):
    HY = entro_calculator(response)
    uni = sorted(np.unique(predictor))
    HYX = []
    prob_x = []
    for x in uni:
        prob_y = []
        y_zero_count = 0
        y_one_count = 0
        x_count = 0
        for i in range(len(predictor)):
            if predictor[i] == x:
                x_count += 1
                if response[i][0] == '0':
                    y_zero_count += 1
                else:
                    y_one_count +=1
        prob_y.append(y_zero_count/x_count)
        prob_y.append(y_one_count/x_count)
        prob_x.append(x_count/len(predictor))
        sum_yx = 0
        for i in range(len(prob_y)):
            if prob_y[i] != 0:
                sum_yx += -prob_y[i] * np.log2(prob_y[i])
            else:
                sum_yx += 0
        HYX.append(sum_yx)
    sum_hyx = 0
    for i in range(len(HYX)):
        sum_hyx += prob_x[i]*HYX[i]
    return HY-sum_hyx, HYX

def majority_vote(y):
    zero_count = 0
    one_count = 0
    uni = sorted(np.unique(y))
    for i in y:
        if i == uni[0]:
            zero_count +=1
        else:
            one_count += 1
    if zero_count > one_count:
        return uni[0]
    else:
        return uni[1]

class dctr():
    #Calculate every attribute's MI
    def MI_selector(self, data):
        mi = []
        for i in range(data.shape[1]-1):
            mi.append(mi_calculator(data[:,i], data[:,-1])[0])
        mi = np.array(mi)
        max_mi_index = np.argmax(mi)
        return max_mi_index
    #Split data given an attribute
    def split(self, data, used_attr_index):
        data_new_0 = []
        data_new_1 = []
        uni = sorted(np.unique(data[:,used_attr_index]))
        #print(uni)
        for i in range(data.shape[0]):
            if data[i,used_attr_index] == uni[0]:
                data_new_0.append(data[i,:])
            else:
                data_new_1.append(data[i,:])
        data_new_0 = np.array(data_new_0)
        data_new_1 = np.array(data_new_1)
        return data_new_0, data_new_1
    #Construct tree
    def build_tr(self, data, max_depth, depth=0, used_attr_index = 0):
        H = []
        used_attr_index = self.MI_selector(data)
        uni = sorted(np.unique(data[:,used_attr_index]))
        for i in range(data.shape[1]-1):
            H.append(mi_calculator(data[:,i], data[:,-1])[1])
        if depth == max_depth or max_depth == 0 or depth > data.shape[1]-1:
            mv = majority_vote(data[:,-1])
            node = Node(data, depth, used_attr_index, mv)
        elif 0.0 in H[used_attr_index]:
            depth += 1
            data_0 = self.split(data,used_attr_index)[0]
            data_1 = self.split(data,used_attr_index)[1]
            node = Node(data, depth, used_attr_index, data[0,-1], None, None, attributes[used_attr_index])
            node.left = Node(data_0, max_depth, depth, data_0[0,-1], None, None, attributes[self.MI_selector(data_0)])
            node.right = Node(data_1, max_depth, depth, data_1[0,-1], None, None, attributes[self.MI_selector(data_1)])
        else:
            node = Node(data, depth, used_attr_index, 0, None, None, attributes[used_attr_index])
            depth += 1
            node.left = self.build_tr(self.split(data,used_attr_index)[0], max_depth, depth, used_attr_index)
            node.right = self.build_tr(self.split(data,used_attr_index)[1], max_depth, depth, used_attr_index)
        return node

def dictionarize(set1, set2):
    key = []
    value = []
    d = dict(enumerate(set2.flatten(), 0))
    for i in range(len(set1)):
        key.append(set1[i])
        value.append(set2[i])
    for j in range(len(key)):
        d[key[j]] = d[j]
        del d[j]
    return d


def predict(node, dictionary, dataset):
    uni = sorted(np.unique(dataset[1,:]))
    if node.attr == None:
        return node.vote
    elif dictionary[node.attr] == uni[0] and node.left != None:
        return predict(node.left, dictionary, dataset)
    elif dictionary[node.attr] == uni[1] and node.right != None:
        return predict(node.right, dictionary, dataset)

def error_rate(node, data):
    pred = []
    true = []
    for i in range(data.shape[0]):
        pred.append(predict(node, dictionarize(attributes,train[i,:]),train))
        true.append(data[i,-1])
    error_number = 0
    for i in range(len(pred)):
        if pred[i] != true[i]:
            error_number += 1
    return pred, error_number / len(pred)

tree1 = dctr()
root = tree1.build_tr(train, max_depth)

def pretty_print(data, node, max_depth):
    uni = sorted(np.unique(data[:,1]))
    tree = dctr()
    print('[' + str(np.count_nonzero(data[:,-1] == uni[0])) + ' ' + uni[0] + '/' + str(np.count_nonzero(data[:,-1] == uni[1])) + ' ' + uni[1] + ']')
    print_part_two(data, node, max_depth)
def print_part_two(data, node, max_depth):
    for i in range(max_depth):
        if i == node.depth:
            print('|'*(i+1) + ' ' + node.attr + ' ' + '=' + ' ' + uni[0]+':'
                  + '[' + str(np.count_nonzero(node.left.data_new[:,-1] == uni[0]))
                  + ' ' + uni[0]+'/'
                  + str(np.count_nonzero(node.left.data_new[:,-1] == uni[1])) + ' ' + uni[1]+']')
            print_part_two(node.left.data_new, node.left, max_depth)
            print('|'*(i+1) + ' ' + node.attr + ' ' + '=' + ' ' + uni[1]+':'
                  + '[' + str(np.count_nonzero(node.right.data_new[:,-1] == uni[0]))
                  + ' ' + uni[0]+'/'
                  + str(np.count_nonzero(node.right.data_new[:,-1] == uni[1])) + ' ' + uni[1]+']')
            print_part_two(node.right.data_new, node.right, max_depth)


f = open(args[4],"w+")
for result in error_rate(root, train)[0]:
    f.write(str(result))
    f.write('\n')
f.close()

f = open(args[5],"w+")
for result in error_rate(root, test)[0]:
    f.write(str(result))
    f.write('\n')
f.close()

f= open(args[6],"w+")
f.write('error(train):' + ' ' + str(error_rate(root, train)[1]))
f.write('\n')
f.write('error(test):' + ' ' + str(error_rate(root, test)[1]))
f.close()





if __name__ == '__main__':
    pass
