import matplotlib as mpl
mpl.use("TkAgg")
import matplotlib.pyplot as plt
from pyMCDS import pyMCDS

def onclick(event):
    print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          (event.button, event.x, event.y, event.xdata, event.ydata))
    plt.plot(event.xdata, event.ydata, 'o',ms=4)
    xs.append(event.xdata)
    ys.append(event.ydata)
    fig.canvas.draw()

def test_saving_parent_ID():
    mcds1 = pyMCDS('output00000001.xml', '..\PhysiCell_V_1.10.4\output')
    #print(mcds1.get_cell_df())
    print()
    cell_df = mcds1.get_cell_df()
    fig,ax = plt.subplots()
    ax.scatter(cell_df['position_x'],cell_df['position_y'],color='r')
    print(mcds1.get_2D_mesh())
    xmin = mcds1.get_2D_mesh()[0][0][0] - mcds1.get_mesh_spacing()/2
    xmax = mcds1.get_2D_mesh()[0][0][-1] + mcds1.get_mesh_spacing()/2
    ymin = mcds1.get_2D_mesh()[1][0][0] - mcds1.get_mesh_spacing() / 2
    ymax = mcds1.get_2D_mesh()[1][-1][-1] + mcds1.get_mesh_spacing() / 2

    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_xlim([0, 10])
# ax.set_ylim([0, 10])
# xs = []
# ys = []
#
# cid = fig.canvas.mpl_connect('button_press_event', onclick)
# plt.show()

test_saving_parent_ID()
plt.show()



