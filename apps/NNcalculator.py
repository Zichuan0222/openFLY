import numpy as np
import math
import csv


filename = "/home/zichuan/openFLY/build/data/NNlist.csv"

coor1 = np.zeros(10, 10, 10)
coor2 = np.zeros(10, 10, 10)

for i in range(10):
    x1 = i
    x2 = i + 0.5

    for j in range(10):
        y1 = j
        y2 = j + 0.5

        for k in range(10):
            z1 = k
            z2 = k + 0.5

            coor1[i][j][k] = [x1, y1, z1]
            coor2[i][j][k] = [x2, y2, z2]


                    
# writing to csv file 
with open(filename, 'w') as csvfile: 
    # creating a csv writer object 
    csvwriter = csv.writer(csvfile) 

    # writing the data rows 
    csvwriter.writerows([str(x1), str(y1), str(z1)])
    csvwriter.writerows([str(x2), str(y2), str(z2)])

