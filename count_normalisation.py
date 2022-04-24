import os
import math

dir = '/home/bc234/test_dir'
count = 0
for file in os.listdir(dir):
    if file.endswith('.txt'):
        count = count + 1
        fo = open(os.path.join(dir, file))
        rl = fo.read().splitlines()
        
        # sum of col
        def col_sum():
            sum = 0
            for val in rl:
                sum = float(sum) + float(val)
                sum = str(sum) 
            return sum
        
        # each count normalised (CPM)
        c = open('/home/bc234/test_dir/file_' + str(count) + 'NORM.txt', 'w')
        for i in rl:
            if i == '0':
              continue
            s = float(i) / float(col_sum())
            s = s * 1000000
            d = math.log(s, 10)
            e = str(s) + '\n'
            c.write(e)
