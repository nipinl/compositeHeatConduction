import matplotlib.pyplot as plt
import numpy as np
geom = []
temp = []

fp = open("geometry","r+")
for i, line in enumerate(fp):
    if i == 0:
        geom.append(line.split())
fp.close()

#getting data
shells = geom[0][0]
Length = geom[0][1]
Width = geom[0][2]
N = int(geom[0][3])
M = int(geom[0][4])
times = int(geom[0][5])

x = np.linspace(0, 4, 3)
y = np.linspace(0, 2, 5)
[X, Y] = np.meshgrid(x,y)
enableColorBar=0

for n in range(times-1):
    fp = open("temperature","r+")
    for i, line in enumerate(fp):
        if(i>n*5 and i < (n+1)*5+1):
            temp.append(line.split())
            #print(i)
    #fig, ax = plt.subplots(1, 1)
    plt.contourf(X,Y,temp)
    """ plt.set_title('Filled Contour Plot')
    plt.set_xlabel('feature_x')
    plt.set_ylabel('feature_y') """
    if enableColorBar==0:
        enableColorBar=1
        plt.colorbar()
    plt.savefig(str(n).zfill(5) )
    temp.clear()
    fp.close()
#print(temp)
#temp.clear()
