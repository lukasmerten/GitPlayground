import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

def clip(x, a, b):
    if x>=b:
        return b
    if x<=a:
        return a
    else:
        return x

def changeHeight(X, h):
    """Change the height of a random section
    a 2-dimensional array. The section is defined by a line which 
    conects two random points of the landscape array."""
    try:
        a, b = np.shape(X)
    except ValueError, e:
        print 'The shape of X has to be a 2-tuple \nInstead X has a shape of '+str(np.shape(X))
    
    x1, y1 = np.random.rand(2)*a
    x2, y2 = np.random.rand(2)*b
    
    m = (y2-y1)/abs(x2-x1)
    y0 = x2-m*x1
    
    for i, c in enumerate(X.transpose()[:,]):
        yi = int(clip(i*m+y0, 0, b))
          
        cAdd = [1]*yi+[0]*(b-yi)
        
        c += np.array(cAdd)*h
        
    return X

if __name__ == "__main__":    
    A = np.zeros([100, 100])
    N = 2500
    for i in range(N):
        # Change height by +1 or -1
        A = changeHeight(A, (0.5-np.random.rand())*2.)
        # use a normal distributed height change for realistic mountain
        A = changeHeight(A, np.random.normal(0,1))
    if i%(N/100)==0:
        print i/(N/100)
        
    # Plot the result
    X, Y = np.meshgrid(range(100), range(100))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(X, Y, A, cmap='gist_earth', linewidth=0., rstride=2, 
    cstride=2, edgecolor=None)
    plt.show()
