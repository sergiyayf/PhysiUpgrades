import matplotlib as mpl
mpl.use('TkAgg')
from pyMCDS import pyMCDS
from pyMCDS_timeseries import pyMCDS_timeseries
#from my_clickablempl import find_nearest_cell as get_my_cell
import numpy as np
import matplotlib.pyplot as plt

def onclick(event):
    print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          (event.button, event.x, event.y, event.xdata, event.ydata))
    plt.plot(event.xdata, event.ydata, marker='o',color='k',ms=6,alpha=0.4)
    xs.append(event.xdata)
    ys.append(event.ydata)
    tracking(xs[-1], ys[-1])

    fig.canvas.draw()
def tracking(clicked_x, clicked_y):
    clicked_x = clicked_x
    clicked_y = clicked_y
    number_of_files = len(timeseries) - 1
    cell_df = timeseries[number_of_files].get_cell_df()
    my_cell = get_my_cell(cell_df, clicked_x, clicked_y)
    # Get number of my daughters
    max_ID = np.max(cell_df['ID'])


    # Loop through other times to recreate a tree + positional progression
    while number_of_files > 0:
        # load older cell data
        number_of_files -= 1
        older_mcds = timeseries[number_of_files]

        older_cell_df = older_mcds.get_cell_df()
        number_of_daughters = get_number_of_daughters(cell_df, my_cell['ID'], max_ID)

        # check if the cell of interest is in the older frame

        if my_cell['ID'] in older_cell_df['ID']:
            # get the cell of interest in the older frame

            my_cell_older = older_cell_df.loc[np.where(older_cell_df['ID'] == my_cell['ID'])[0][0]]
            my_new_number_of_daughters = get_number_of_daughters(older_cell_df, my_cell_older['ID'], max_ID)

            # connect my cell with itself older
            connect(ax, my_cell, my_cell_older)
            # check if number of daughter cells decreased, if yes, mark it as a division point

            if my_new_number_of_daughters < number_of_daughters:

                mark_as_division(ax, my_cell)

        # If the cell of interest is not in the older frame
        else:

            my_cell_older = older_cell_df.loc[np.where(older_cell_df['ID'] == my_cell['parent_ID'])[0][0]]
            connect(ax, my_cell, my_cell_older)
            mark_as_division(ax, my_cell)
            max_ID = my_cell['ID']

        my_cell = my_cell_older
        cell_df = older_cell_df
def get_my_cell(cell_df, x, y):
    xs = cell_df['position_x']
    ys = cell_df['position_y']
    dist = np.sqrt((xs - x) ** 2 + (ys - y) ** 2)

    return cell_df.loc[np.argmin(dist)]
def get_number_of_daughters(cell_df,me,max):

    cells = cell_df.iloc[np.where(cell_df['parent_ID']==me)]

    number_of_daughters = len(np.where(cells['ID']<max)[0])
    return number_of_daughters
def connect(ax,cell_a,cell_b):
    ax.plot([cell_a['position_x'],cell_b['position_x']],[cell_a['position_y'],cell_b['position_y']],color=plt.get_cmap("Set1").colors[int(cell_a['cell_type'])])
    ax2.plot([cell_a['position_x'],cell_b['position_x']],[cell_a['position_y'],cell_b['position_y']],color=plt.get_cmap("Set1").colors[int(cell_a['cell_type'])])

    return
def mark_as_division(ax,cell):
    ax.scatter(cell['position_x'],cell['position_y'],color = 'k',s=7,alpha=0.5)
    ax2.scatter(cell['position_x'], cell['position_y'], color='k',s=7, alpha=0.5)
    #ax2.text(cell['position_x'], cell['position_y'], 'ID: '+str(cell['ID']))
    #ax2.text(cell['position_x'], cell['position_y']-10, 'pID: '+ str(cell['parent_ID']))
    return
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
        plot_one_cell_type(cell_df,ax,c_type,color=plt.get_cmap("tab10").colors[idx])

    xmin = mcds.get_2D_mesh()[0][0][0] - mcds.get_mesh_spacing() / 2
    xmax = mcds.get_2D_mesh()[0][0][-1] + mcds.get_mesh_spacing() / 2
    ymin = mcds.get_2D_mesh()[1][0][0] - mcds.get_mesh_spacing() / 2
    ymax = mcds.get_2D_mesh()[1][-1][-1] + mcds.get_mesh_spacing() / 2

    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.set_aspect('equal')
    ax2.set_xlim([xmin, xmax])
    ax2.set_ylim([ymin, ymax])
    ax2.set_aspect('equal')
    return

if __name__ == '__main__':
    # Load the output
    fig, (ax, ax2) = plt.subplots(1, 2)

    timeseries = pyMCDS_timeseries('..\PhysiCell_V_1.10.4\output').timeseries
    number_of_files = len(timeseries) - 1
    mcds = timeseries[number_of_files]
    cell_df = mcds.get_cell_df()
    plot_all_cell_types(ax, mcds)
    xs = []
    ys = []
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    # Click on cell of interest and get it
