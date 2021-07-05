import matplotlib.pyplot as plt
import numpy as np
geom = []
temp = []
fp = open("temperature","r+")
for i, line in enumerate(fp):
    if i == 0:
        geom.append(line.split())
fp.close()
Length = geom[0][0]
Width = geom[0][1]
N = int(geom[0][2])
M = int(geom[0][3])
times = int(geom[0][4])
print(geom)
x = np.linspace(0, 4, 3)
y = np.linspace(0, 2, 5)
[X, Y] = np.meshgrid(x,y)
for n in range(times):
    fp = open("temperature","r+")
    for i, line in enumerate(fp):
        if(i>n*5 and i < (n+1)*5+1):
            temp.append(line.split())
            print(i)
    print(n)
    print(temp)
    fig, ax = plt.subplots(1, 1)
    ax.contourf(X,Y,temp)
    ax.set_title('Filled Contour Plot')
    ax.set_xlabel('feature_x')
    ax.set_ylabel('feature_y')
    plt.savefig(str(n).zfill(5) )
    temp.clear()
    fp.close()
#print(temp)
#temp.clear()
