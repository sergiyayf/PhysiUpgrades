import matplotlib as mpl
mpl.use('TkAgg')
from pyMCDS import pyMCDS
from pyMCDS_timeseries import pyMCDS_timeseries
#from my_clickablempl import find_nearest_cell as get_my_cell
import numpy as np

def get_my_cell(cell_df,x,y):
    return cell_df.loc[10]
def get_number_of_daughters(me,max):
    return 0
def connect(a,b):
    return
def mark_as_division():
    return
# Load the output

timeseries = pyMCDS_timeseries('..\PhysiCell_V_1.10.4\output').timeseries
number_of_files = len(timeseries) -1
mcds = timeseries[number_of_files]
cell_df = mcds.get_cell_df()

# Click on cell of interest and get it
clicked_x = 110
clicked_y = 110
my_cell = get_my_cell(cell_df,clicked_x,clicked_y)
# Get number of my daughters
max_ID = 100
number_of_daughters = get_number_of_daughters(my_cell['ID'], max_ID)

# Loop through other times to recreate a tree + positional progression
while number_of_files>=0:
    # load older cell data
    number_of_files-=1
    older_mcds = timeseries[number_of_files]
    older_cell_df = older_mcds.get_cell_df()

    # check if the cell of interest is in the older frame
    if my_cell['ID'] in older_cell_df['ID']:
        # get the cell of interest in the older frame

        my_cell_older = older_cell_df.loc[np.where(older_cell_df['ID']==my_cell['ID'])[0][0]]
        my_new_number_of_daughters = get_number_of_daughters(my_cell_older['ID'],my_cell['ID'])

        # connect my cell with itself older
        connect(my_cell, my_cell_older)
        # check if number of daughter cells decreased, if yes, mark it as a division point
        if my_new_number_of_daughters < number_of_daughters:
            mark_as_division()

    # If the cell of interest is not in the older frame
    else:
        my_cell_older = older_cell_df.loc[np.where(older_cell_df['ID']==my_cell['parent_ID'])[0][0]]
        connect(my_cell, my_cell_older)
        mark_as_division()

    my_cell = my_cell_older
    print(my_cell['parent_ID'])


