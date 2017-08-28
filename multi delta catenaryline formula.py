import matplotlib.pylab as plt
import math
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
##delta=[10,20,40,60]
delta=50
delta_seq=[2,5,10,30,50,70]
gama=0.05
e=math.e
def ch(x):
    return((e**x+e**(-x))/2)
def sh(x):
    return((e**x-e**(-x))/2)

#L=20/math.cos(math.atan(0.8/60))
L=30
x=np.arange(0,L+0.1,0.1)
for delta in delta_seq:
    t=delta/gama
    h=13
    #悬链线
##    Loa=L/2-t*sh(h/2/t/sh(L/2/t))
##    print('Loa',Loa)
##    y=t*(ch((Loa-x)/t)-ch(Loa/t))

    #斜抛线
    beta=math.atan(h/L)
    Loa=L/2-t*math.sin(beta)
    y=x*math.tan(beta)-x*(L-x)/2/t/math.cos(beta)

    print(y[-1])
    plt.plot(x,y)
plt.show()
