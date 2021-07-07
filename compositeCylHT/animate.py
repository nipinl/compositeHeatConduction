import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
geom = []
temp = []

fp = open("s.geom","r+")
for i, line in enumerate(fp):
    if i == 0:
        geom.append(line.split())
fp.close()

#getting data
shells = geom[0][0]
Length = float(geom[0][1])
Width = float(geom[0][2])
N = int(geom[0][3])
M = int(geom[0][4])
times = int(geom[0][5])

x = np.linspace(0, Length, N)
y = np.linspace(0, Width, M)
[X, Y] = np.meshgrid(x,y)
enableColorBar=0
n=1

data = np.genfromtxt('s.temp', delimiter=',')
minTemp = data.min() 
maxTemp = data.max() 
def getTempForTime(i):   # function returns temp for time=i
    return data[i*M:(i+1)*M,:]
T = getTempForTime(100)
fig = plt.figure()
ax = plt.axes(xlim=(0, Length), ylim=(0, Width), xlabel='Length', ylabel='Width')
cvals = np.linspace(minTemp,maxTemp,20)      # set contour values 
cont = plt.contourf(X, Y, getTempForTime(0), cvals)    # first image on screen
plt.colorbar()

# animation function
def animate(i):
    global cont
    T = getTempForTime(i)
    for c in cont.collections:
        c.remove()  # removes only the contours, leaves the rest intact
    cont = plt.contourf(X, Y, T, cvals)
    plt.title('t = %i:  %.2f' % (i,times-1))
    return cont

anim = animation.FuncAnimation(fig, animate, frames=times-1, repeat=False)
anim.save('animation.mp4', writer=animation.FFMpegWriter())
#for saving the last frame
T = getTempForTime(times-2)
cont = plt.contourf(X, Y, T, cvals)
plt.title('t = %i:  %.2f' % (i,times-1))
plt.savefig(str(times-1).zfill(5) )

""" for n in range(times-1):
    fp = open("s.temp","r+")
    for i, line in enumerate(fp):
        if(i>=n*M and i < (n+1)*M):
            temp.append(line.split(','))
            #print(i) """
""" plt.contourf(X,Y,T)

plt.colorbar()
plt.savefig(str(n).zfill(5) )

fp.close() """

