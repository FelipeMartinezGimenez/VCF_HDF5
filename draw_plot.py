import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d.axes3d import Axes3D


def crear_grafico3D_angulo(X, Y, Z, fhand, ang1=0, ang2=90, ejes=[1,2,3]):
    
    
    eje1, eje2, eje3=ejes
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X, Y, Z, color='black', marker='+')
    ax.view_init(ang1, ang2)
    canvas = FigureCanvas(fig)
    xn='componente '+str(ejes[0])
    yn='componente '+str(ejes[1])
    zn='componente '+str(ejes[2])
    ax.set_xlabel(xn)
    ax.set_ylabel(yn)
    ax.set_zlabel(zn)
    canvas.print_figure(fhand)
    del ax
    del fig



if __name__=='__main__':
   
    crear_grafico3D_angulo(X, Y, Z, fhand, ang1, ang2)
