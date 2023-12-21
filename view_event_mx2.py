import uproot
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib.widgets import Button, TextBox
import os
import subprocess
import sys



    
class Mx2Data:
    def __init__(self, filename):
        self.file = uproot.open(filename)

        self.hits_id_per_mod = self.file["minerva"]["hits_id_per_mod"].array()
        self.hit_strip = self.file["minerva"]["hit_strip"].array()
        self.hit_plane = self.file["minerva"]["hit_plane"].array()
        self.hit_module = self.file["minerva"]["hit_module"].array()
        self.hit_view = self.file["minerva"]["hit_view"].array()
        self.n_odhits = self.file["minerva"]["n_odhits"].array()
        self.hits_od_per_mod = self.file["minerva"]["hits_od_per_mod"].array()
        self.hit_bar = self.file["minerva"]["hit_bar"].array()
        self.hit_pe = self.file["minerva"]["hit_pe"].array()
        self.hit_time = self.file["minerva"]["hit_time"].array()
        self.hit_time_slice = self.file["minerva"]["hit_time_slice"].array()

        self.offsetX = self.file["minerva"]["offsetX"].array()
        self.offsetY = self.file["minerva"]["offsetY"].array()
        self.offsetZ = self.file["minerva"]["offsetZ"].array()

        self.ev_gps_time_sec = self.file["minerva"]["ev_gps_time_sec"].array()
        self.ev_gps_time_usec = self.file["minerva"]["ev_gps_time_usec"].array()

        self.n_slices = self.file["minerva"]["n_slices"].array()

        self.clus_id_coord = self.file["minerva"]["clus_id_coord"].array()
        self.clus_id_z = self.file["minerva"]["clus_id_z"].array()
        self.clus_id_view = self.file["minerva"]["clus_id_view"].array()
        self.clus_id_pe = self.file["minerva"]["clus_id_pe"].array()

        self.minerva_time = (self.ev_gps_time_sec - self.ev_gps_time_sec[0])*1e6+(self.ev_gps_time_usec)
        self.minerva_trigger = np.array((self.minerva_time/1.2e6)).astype(int)

        self.hits_id_per_mod = self.file["minerva"]["hits_id_per_mod"].array()
        self.hit_strip = self.file["minerva"]["hit_strip"].array()

class NdData:
    def __init__(self, filename):
        self.flow_file = h5.File(filename, 'r')

        self.data = self.flow_file['charge']['calib_prompt_hits']['data']
        self.nd_trigger = (self.data.fields("ts_pps")[:]/10/1.2e6).astype(int)  

# Mx2 conversion 
def strip_to_x(x, offsetY = 0, offsetX = 0, view = 1):
    tot_offset = offsetX
    if (view==2):
        tot_offset = calcUfromXY(offsetX, offsetY)
    if (view==3):
        tot_offset = calcVfromXY(offsetX, offsetY)
    return( 16.738776 * x  -1071.2776 - tot_offset) 

def module_to_z(x, offset = 0):
    if (x < 13):
        return(4001.4545 + x * 43.545455) - offset
    return(44.783333  * x + 8015.4947) - offset


strip_to_x = np.vectorize(strip_to_x)
module_to_z = np.vectorize(module_to_z)

# U and V positions
def calcUfromXY( x, y ) :
    return 0.5*( x - np.sqrt(3)*y )

def calcVfromXY( x, y ) :
    return 0.5*( x + np.sqrt(3)*y )

calcUfromXY = np.vectorize(calcUfromXY)
calcVfromXY = np.vectorize(calcVfromXY)


def plot_view_nd(Mx2Hits, entry, view, ax, first_draw = True,  min_pe=0, slice = 0):
    
    ax.clear()
    mask = (Mx2Hits.hit_module[entry]>0) & (Mx2Hits.hit_view[entry]==view) & (Mx2Hits.hit_pe[entry]>min_pe)
    if (slice>0):
        mask = (Mx2Hits.hit_module[entry]>0) & (Mx2Hits.hit_view[entry]==view) & (Mx2Hits.hit_pe[entry]>min_pe) & (Mx2Hits.hit_time_slice[entry]==slice) 
    

    
            
    if (len(Mx2Hits.hit_module[entry][mask])>0):
        x_bins = module_to_z(Mx2Hits.hit_module[entry][mask], Mx2Hits.offsetZ[entry])
        y_bins = strip_to_x(Mx2Hits.hit_strip[entry][mask], Mx2Hits.offsetY[entry], Mx2Hits.offsetX[entry], view)
        weights = Mx2Hits.hit_pe[entry][mask]
    x_bins = np.array(x_bins)
    y_bins = np.array(y_bins)
    weights = np.array(weights)
    mask = np.array(mask)     
    hist = ax.hist2d( x_bins, y_bins , weights=weights, bins=[module_to_z(np.concatenate((np.concatenate((np.arange(0,13,0.5),[12.5])),np.arange(13,44,.5))), Mx2Hits.offsetZ[entry]),strip_to_x(np.arange(-4, 130,.5), Mx2Hits.offsetY[entry], Mx2Hits.offsetX[entry], view)], cmap='magma_r', cmin=1e-4)
    
    title = "X view"
    y_title = "x [mm]"
    if (view==2):
        title = "U view"
        y_title = "u [mm]"
    if (view==3):
        title = "V view"
        y_title = "v [mm]"
        ax.set_xlabel("z [mm]", fontsize=20)
            
    ax.set_ylabel(y_title, fontsize=20)

    ax.set_title(title, fontsize=20, y=.95, pad=-14)
    
    
    
    hist[3].set_clim(vmin=0, vmax=35)
    
    if (first_draw):
        cl = plt.colorbar(hist[3], ax= ax)
        cl.set_label('Mx2 pe')
        hist[3].set_clim(vmin=0, vmax=35)
    ax.grid(True, axis='y')
    
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)
# Plotting Mx2 hit time    

def plot_mx2_time(Mx2Hits, entry, ax, min_pe=0, slice = 0):
   
    ax.clear()
    time_bin = np.linspace(0,16,1600)
    mask = np.array((Mx2Hits.hit_time_slice[entry] == 0)  & (Mx2Hits.hit_pe[entry]>min_pe))
    # plt.hist(hit_time[i]/1000, bins=time_bin, log=True, histtype='step')4

    if (slice >0): 
        mask =  np.array(Mx2Hits.hit_pe[entry]>min_pe)
        x_bins = np.array(Mx2Hits.hit_time[entry][mask]/1000)

        h0 = ax.hist(x_bins, bins=time_bin, log=True, histtype='stepfilled', color='black', alpha=.9)
        mask = np.array((Mx2Hits.hit_time_slice[entry] == slice) & (Mx2Hits.hit_pe[entry]>min_pe))
        x_bins = np.array(Mx2Hits.hit_time[entry][mask]/1000)

        h = ax.hist(x_bins, bins=time_bin, log=True, histtype='stepfilled', color='red')
        max_bin = h0[0].max()

    else:
        x_bins = np.array(Mx2Hits.hit_time[entry][mask]/1000)
        h0 = ax.hist(x_bins, bins=time_bin, log=True, histtype='stepfilled', color='black', alpha=.9)
        max_bin = h0[0].max()
        for ts in range(1,Mx2Hits.n_slices[entry]+1):
            mask = np.array((Mx2Hits.hit_time_slice[entry] == ts) & (Mx2Hits.hit_pe[entry]>min_pe))
            x_bins = np.array(Mx2Hits.hit_time[entry][mask]/1000)
            if (len(x_bins)):
                h = ax.hist(x_bins, bins=time_bin, log=True, histtype='stepfilled')
                max_bin = max(max_bin, h[0].max())
    
    ax.set_xlim(0,16)
    
    ax.set_xlabel("time [µs]", fontsize=15)
    ax.set_ylabel("hits", fontsize=15)
    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=15)
    ax.xaxis.set_label_coords(0.9, .2)
    
    # ax.set_xticks(20)
    
    ax.set_title("Mx2 timing", fontsize=15)
    ax.text(13, max_bin*12, f'Trigger number: {entry}', dict(size=15))
    ax.text(13, max_bin*2, f'Mx2 Slice number: {slice}', dict(size=15))
    ax.text(0, max_bin*2, f'#Run: {run:05}', dict(size=15))    
    
    ax.grid(which='both', axis="y",linestyle='--', linewidth='0.5', color='black')
    ax.grid(which='major', axis="x",linestyle='--', linewidth='0.5', color='black')
    
    # ax.yaxis.get_label().set_fontsize(20)


    
#main plotter function
def view_event(Mx2Hits, trig, first_draw=True, mx2_min_pe=0, slice=0):
    plot_mx2_time(Mx2Hits, trig,axs[0], mx2_min_pe, slice)
    plot_view_nd(Mx2Hits, trig,1,axs[1], first_draw, mx2_min_pe, slice)
    plot_view_nd(Mx2Hits, trig,2,axs[2], first_draw, mx2_min_pe, slice)
    plot_view_nd(Mx2Hits, trig,3,axs[3], first_draw, mx2_min_pe, slice)
    # f.savefig(f'plots/MnV_only_MR5_trigg_{trig}.png', facecolor='white')    1
    # plt.show()

    return f


def on_next_button_click(event):
    global current_event
    current_event+=1
    first_draw = False

    global current_slice 
    current_slice = 0

    global current_mx2_pe_cut
    current_mx2_pe_cut = 0

    if (current_event >= len(Mx2Hits.ev_gps_time_sec)):
        current_event = len(Mx2Hits.ev_gps_time_sec)-1
    view_event(Mx2Hits, current_event, first_draw)

def on_previous_button_click(event):
    global current_event
    current_event = max(0, current_event-1)
    
    global current_slice 
    current_slice = 0

    global current_mx2_pe_cut
    current_mx2_pe_cut = 0

    first_draw = False
    view_event(Mx2Hits, current_event, first_draw)
    

def on_next_slice_click(event):
    global current_event
    first_draw = False

    global current_mx2_pe_cut
    current_mx2_pe_cut = 0

    global current_slice
    current_slice +=1
    if (current_slice >= Mx2Hits.n_slices[current_event]):
        current_slice = Mx2Hits.n_slices[current_event]

    view_event(Mx2Hits, current_event, first_draw, current_mx2_pe_cut, current_slice)


def on_previous_slice_click(event):
    global current_event

    first_draw = False

    global current_mx2_pe_cut
    current_mx2_pe_cut = 0

    global current_slice
    current_slice -=1

    if (current_slice <=0):
        current_slice = 0

    view_event(Mx2Hits, current_event, first_draw, current_mx2_pe_cut, current_slice)
    


def on_text_entered(event):
    global current_event
    entered_event = int(event)
    current_event = entered_event
    first_draw = False
    global current_mx2_pe_cut
    current_mx2_pe_cut = 0

    global current_slice
    current_slice =0

    if (current_event < 0):
        current_event = 0
    if (current_event >= len(Mx2Hits.ev_gps_time_sec)):
        current_event = len(Mx2Hits.ev_gps_time_sec)-1
    view_event(Mx2Hits, current_event, first_draw, current_mx2_pe_cut, current_slice)

    

if len(sys.argv) != 3:
        print("Usage: python script_name.py <run> <trigger>")
        sys.exit(1)
try:
    # Attempt to convert the command-line argument to an integer
    run = int(sys.argv[1])
    trig = int(sys.argv[2])
except ValueError:
    print("Error: The argument must be an integer.")
    sys.exit(1)


current_event = trig
current_slice = 0
f, axs = plt.subplots(4, 1, figsize=(15,10), gridspec_kw={'height_ratios' : [.5,1,1,1]})
# Create a Button widget
next_button_ax = plt.axes([0.9, 0.09, 0.07, 0.05])  # [left, bottom, width, height]
next_button = Button(next_button_ax, 'Next Event')

before_button_ax = plt.axes([0.82, 0.09, 0.07, 0.05])  # [left, bottom, width, height]
before_button = Button(before_button_ax, 'Previous Event')

next_slice_ax = plt.axes([0.9, 0.01, 0.07, 0.05])  # [left, bottom, width, height]
next_slice = Button(next_slice_ax, 'Next slice')

before_slice_ax = plt.axes([0.82, 0.01, 0.07, 0.05])  # [left, bottom, width, height]
before_slice = Button(before_slice_ax, 'Previous slice')


# Connect the button to the function
next_button.on_clicked(on_next_button_click)
before_button.on_clicked(on_previous_button_click)

next_slice.on_clicked(on_next_slice_click)
before_slice.on_clicked(on_previous_slice_click)

# Create a TextBox widget for entering event numbers
text_box_ax = plt.axes([0.86, 0.17, 0.06, 0.05])  # [left, bottom, width, height]
event_number_input = TextBox(text_box_ax, 'Event:')
event_number_input.on_submit(on_text_entered)


Mx2Hits = Mx2Data(f"output_dir/run_{run:05}/dst/MiniRun5_1E19_RHC.minerva.{run:05}.dst.root")

f = view_event(Mx2Hits, current_event)
directory_path=f'output_dir/plots/{run:05}'

if not os.path.exists(directory_path):
    os.makedirs(directory_path)


f.savefig(f'output_dir/plots/{run:05}/MiniRun5_1E19_RHC.minerva{run:05}.{trig}.png', dpi=60, facecolor='white')
