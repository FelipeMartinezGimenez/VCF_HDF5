
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from mpl_toolkits.mplot3d import Axes3D
import os




def draw_pca(projections, color='blue', dir_ ='.'):
    X, Y, Z = projections
    
    fig = Figure()
    
    axes = fig.add_subplot(111, projection='3d', xticks=[], yticks=[], zticks=[])
    axes.scatter3D(X, Y, Z, color=color)
    axes.axis("off")
       
    for angle in range(0, 360, 10):
        fname = 'image%04d.png' % (angle)
        fpath = os.path.join(dir_, fname)
        fhand = open(fpath, 'w')
        
        axes.view_init(azim=angle)
        canvas = FigureCanvas(fig)
        canvas.print_figure(fhand)
        fhand.close()    
    
