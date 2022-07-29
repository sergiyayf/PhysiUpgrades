import matplotlib as mpl
mpl.use("TkAgg")
import matplotlib.pyplot as plt
from pyMCDS import pyMCDS
from matplotlib.widgets import Button
from matplotlib.text import Annotation
from matplotlib.patches import ConnectionPatch
import numpy as np
from adjustText import adjust_text

def annotate(axis, text, x, y, parentx, parenty):
    text_annotation = Annotation(text, xy=(x, y), xycoords='data')
    axis.add_artist(text_annotation)
    if(parentx != -1):
        xyA = (x, y)
        xyB = (parentx, parenty)
        coordsA = "data"
        coordsB = "data"
        con = ConnectionPatch(xyA, xyB, coordsA, coordsB,
                              arrowstyle="-|>", shrinkA=0, shrinkB=0,
                              mutation_scale=10, fc="w")
        axis.add_artist(con)
# define the behaviour -> what happens when you pick a dot on the scatterplot by clicking close to it
def onpick(event):
    # step 1: take the index of the dot which was picked
    ind = event.ind

    # step 2: save the actual coordinates of the click, so we can position the text label properly
    label_pos_x = event.mouseevent.xdata
    label_pos_y = event.mouseevent.ydata

    # just in case two dots are very close, this offset will help the labels not appear one on top of each other
    offset = 0

    # if the dots are to close one to another, a list of dots clicked is returned by the matplotlib library
    for i in ind:
        id = ids[int(i)]
        parentid = parents[int(i)]
        if(parentid >= 0):
            parentpos = ids.index(int(parentid))
        else:
            parentpos = 0
        # step 3: take the label for the corresponding instance of the data
        if parentid != -1:
            label = "me: "+str(id)+"\nparent: "+str(parentid)
            annotate(
                ax,
                label,
                label_pos_x + offset,
                label_pos_y + offset,
                positionx[parentpos],
                positiony[parentpos]
            )
        else:
            label = "me: "+str(id)+"\nparent: none"
            annotate(
                ax,
                label,
                label_pos_x + offset,
                label_pos_y + offset,
                -1,
                -1
            )
        # step 6: force re-draw
        # adjust_text(label, expand_points=(1, 1), expand_text=(1, 1),
        #             arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
        ax.figure.canvas.draw_idle()

        # alter the offset just in case there are more than one dots affected by the click
        offset += 0.01



def test_saving_parent_ID():
    mcds1 = pyMCDS('output00000020.xml', '../PhysiCell_V_1.10.4/output')
    #print(mcds1.get_cell_df())
    print()
    cell_df = mcds1.get_cell_df()
    ax.scatter(cell_df['position_x'], cell_df['position_y'], color='r', picker=True, pickradius=5)
    print(mcds1.get_2D_mesh())
    xmin = mcds1.get_2D_mesh()[0][0][0] - mcds1.get_mesh_spacing()/2
    xmax = mcds1.get_2D_mesh()[0][0][-1] + mcds1.get_mesh_spacing()/2
    ymin = mcds1.get_2D_mesh()[1][0][0] - mcds1.get_mesh_spacing() / 2
    ymax = mcds1.get_2D_mesh()[1][-1][-1] + mcds1.get_mesh_spacing() / 2

    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])

# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_xlim([0, 10])
# ax.set_ylim([0, 10])
xs = []
ys = []
points = []
# cid = fig.canvas.mpl_connect('button_press_event', onclick)
# plt.show()
fig,ax = plt.subplots()
cid = fig.canvas.mpl_connect('pick_event', onpick)
mcds1 = pyMCDS('output00000020.xml', '../PhysiCell_V_1.10.4/output')
cell_df = mcds1.get_cell_df()
points = []
parents=list(cell_df["parent_ID"])
ids =list(cell_df["ID"])
positionx = list(cell_df["position_x"])
positiony = list(cell_df["position_y"])
test_saving_parent_ID()
plt.show()



