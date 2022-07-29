import matplotlib as mpl
mpl.use("TkAgg")
import matplotlib.pyplot as plt
import numpy as np
from pyMCDS import pyMCDS

def onclick(event):
    print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          (event.button, event.x, event.y, event.xdata, event.ydata))
    plt.plot(event.xdata, event.ydata, marker='o',color='k',ms=6,alpha=0.4)
    xs.append(event.xdata)
    ys.append(event.ydata)
    print(xs)
    nearest_cell = find_nearest_cell(cell_df, xs[-1], ys[-1])
    connect_to_parents(nearest_cell, cell_df, ax=ax)

    fig.canvas.draw()


def plot_one_cell_type(cell_df, ax=None, c_type=0, color='grey', edgecolor='dimgrey'):
    x_pos = cell_df['position_x']
    y_pos = cell_df['position_y']
    cell_type = cell_df['cell_type']
    rad = (cell_df['total_volume'] * 3 / 4 / np.pi) ** (1 / 3)
    circles = [plt.Circle((xi, yi), radius=ri, linewidth=0) for xi, yi, ri in
               zip(x_pos[cell_type == c_type], y_pos[cell_type == c_type], rad[cell_type == c_type])]
    c = mpl.collections.PatchCollection(circles, color=color, alpha=0.75, edgecolor=edgecolor, linewidth=0.01)
    ax.scatter(x_pos[cell_type == c_type], y_pos[cell_type == c_type], color=color, s=1.5, alpha=0.75, linewidth=0.00,
               edgecolor=edgecolor)

    ax.add_collection(c)
    return

def plot_all_cell_types(ax,mcds):
    cell_df = mcds.get_cell_df()
    cell_types = np.array(cell_df['cell_type'])

    for idx,c_type in enumerate(np.unique(cell_types)):
        plot_one_cell_type(cell_df,ax,c_type,color=plt.get_cmap("tab20").colors[idx])

    xmin = mcds1.get_2D_mesh()[0][0][0] - mcds1.get_mesh_spacing() / 2
    xmax = mcds1.get_2D_mesh()[0][0][-1] + mcds1.get_mesh_spacing() / 2
    ymin = mcds1.get_2D_mesh()[1][0][0] - mcds1.get_mesh_spacing() / 2
    ymax = mcds1.get_2D_mesh()[1][-1][-1] + mcds1.get_mesh_spacing() / 2

    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.set_aspect('equal')
    return

def find_nearest_cell(cell_df,x,y):
    xs = cell_df['position_x']
    ys = cell_df['position_y']
    dist = np.sqrt((xs-x)**2+(ys-y)**2)

    return cell_df.loc[np.argmin(dist)]

def connect_to_parents(cell,cell_df,ax = None):
    parent_ID = cell['parent_ID']

    while parent_ID >= 0:
        parent_idx = np.where(cell_df['ID']==cell['parent_ID'])[0][0]
        parent = cell_df.loc[parent_idx]
        x = cell['position_x']
        y = cell['position_y']
        px = parent['position_x']
        py = parent['position_y']
        ax.plot([x,px],[y,py],color='k')
        cell = parent
        parent_ID = cell['parent_ID']
        # if parent_ID < 0:
        #     break
    return

def test_saving_parent_ID(ax):
    mcds1 = pyMCDS('output00000107.xml', '../PhysiCell_V_1.10.4/output')
    print(mcds1.get_cell_df())
    print()
    cell_df = mcds1.get_cell_df()

    ax.scatter(cell_df['position_x'],cell_df['position_y'],color=cell_coloring_function(cell_df['cell_type'].values()))
    print(mcds1.get_2D_mesh())
    xmin = mcds1.get_2D_mesh()[0][0][0] - mcds1.get_mesh_spacing()/2
    xmax = mcds1.get_2D_mesh()[0][0][-1] + mcds1.get_mesh_spacing()/2
    ymin = mcds1.get_2D_mesh()[1][0][0] - mcds1.get_mesh_spacing() / 2
    ymax = mcds1.get_2D_mesh()[1][-1][-1] + mcds1.get_mesh_spacing() / 2

    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
#

xs = []
ys = []
# test_saving_parent_ID(ax)
# cid = fig.canvas.mpl_connect('button_press_event', onclick)
mcds1 = pyMCDS('output00000107.xml', '../PhysiCell_V_1.10.4/output')
cell_df = mcds1.get_cell_df()
# a = find_nearest_cell(cell_df,1.0,1.0)
# connect_to_parents(a,cell_df,ax=ax)


fig, ax = plt.subplots(figsize=(10,10),dpi=150)
plot_all_cell_types(ax,mcds1)
cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show()


