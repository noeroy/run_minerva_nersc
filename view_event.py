import uproot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib.widgets import Button, TextBox
from ipywidgets import interact, Dropdown

import h5py as h5
import time



    
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


    
def plot_view_nd(Mx2Hits, NdFlow, trig_num, view, ax, first_draw = True,  min_pe=0, slice = 0):
    
    ax.clear()
    entry = Mx2Hits.minerva_trigger[trig_num]
    mask = (Mx2Hits.hit_module[entry]>0) & (Mx2Hits.hit_view[entry]==view) & (Mx2Hits.hit_pe[entry]>min_pe)
    if (slice>0):
        mask = (Mx2Hits.hit_module[entry]>0) & (Mx2Hits.hit_view[entry]==view) & (Mx2Hits.hit_pe[entry]>min_pe) & (Mx2Hits.hit_time_slice[entry]==slice)
    
    x_nd = NdFlow.data.fields("x")[NdFlow.nd_trigger==trig_num]*10
    y_nd = NdFlow.data.fields("y")[NdFlow.nd_trigger==trig_num]*10
    z_nd = NdFlow.data.fields("z")[NdFlow.nd_trigger==trig_num]*10
        
    Q_nd = NdFlow.data.fields("Q")[NdFlow.nd_trigger==trig_num]
    time = NdFlow.data.fields("ts_pps")[NdFlow.nd_trigger==trig_num]/10 - 1.2e6 * NdFlow.nd_trigger[NdFlow.nd_trigger==trig_num]
    mask1 = (time<10e3)
    mask2 = (time>10e3) 

    
    
    coord =  x_nd
    title = "X view"
    y_title = "x [mm]"
    if (x_nd.size>0):
        if (view==2):
            coord =  calcUfromXY(x_nd, y_nd)
            title = "U view"
            y_title = "u [mm]"
        if (view==3):
            coord =  calcVfromXY(x_nd, y_nd)
            title = "V view"
            y_title = "v [mm]"
            ax.set_xlabel("z [mm]", fontsize=20)

        hist_nd = ax.hist2d(z_nd[mask1], coord[mask1], weights=Q_nd[mask1],bins=[300,300], cmap='magma_r', cmin=1e-4)
        hist_nd2 = ax.hist2d(z_nd[mask2], coord[mask2], weights=Q_nd[mask2],bins=[300,300], cmap='viridis_r', cmin=1e-4)
        hist_nd[3].set_clim(vmin=0, vmax=50)
        hist_nd2[3].set_clim(vmin=0, vmax=50)
        if(view==1 and first_draw):
            cbar = plt.colorbar(hist_nd[3], ax=ax)
            cbar.set_label('ND Q')
            hist_nd[3].set_clim(vmin=0, vmax=50)

        if (view==2 and first_draw):
            cbar2 = plt.colorbar(hist_nd2[3], ax=ax)
            cbar2.set_label('ND Q delayed')
            hist_nd2[3].set_clim(vmin=0, vmax=30)
        
    x_bins = []
    y_bins = []
    weights = []
    if (len(Mx2Hits.hit_module[entry][mask])>0):
        x_bins = module_to_z(Mx2Hits.hit_module[entry][mask], Mx2Hits.offsetZ[entry])
        y_bins = strip_to_x(Mx2Hits.hit_strip[entry][mask], Mx2Hits.offsetY[entry], Mx2Hits.offsetX[entry], view)
        weights = Mx2Hits.hit_pe[entry][mask]
        
    hist = ax.hist2d( x_bins, y_bins , weights=weights, bins=[module_to_z(np.concatenate((np.concatenate((np.arange(0,13,0.5),[12.5])),np.arange(13,44,.5))), Mx2Hits.offsetZ[entry]),strip_to_x(np.arange(-4, 130,.5), Mx2Hits.offsetY[entry], Mx2Hits.offsetX[entry], view)], cmap='magma_r', cmin=1e-4)
    
    
    ax.set_ylabel(y_title, fontsize=20)
    ax.set_title(title, fontsize=20, y=.95, pad=-14)
    
    
    
    hist[3].set_clim(vmin=0, vmax=35)
    
    if (view==3 and first_draw):
        cl = plt.colorbar(hist[3], ax= ax)
        cl.set_label('Mx2 pe')
        hist[3].set_clim(vmin=0, vmax=35)
    ax.grid(True, axis='y')
    
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)
    # plt.draw()

# Plotting Mx2 hit time    

def plot_mx2_time(Mx2Hits, entry, ax, min_pe=0, slice = 0):
   
    ax.clear()
    time_bin = np.linspace(0,16,1600)
    mask = (Mx2Hits.hit_time_slice[entry] == 0)  & (Mx2Hits.hit_pe[entry]>min_pe)
    # plt.hist(hit_time[i]/1000, bins=time_bin, log=True, histtype='step')4

    if (slice >0): 
        mask =  (Mx2Hits.hit_pe[entry]>min_pe)
        h0 = ax.hist(Mx2Hits.hit_time[entry][mask]/1000, bins=time_bin, log=True, histtype='stepfilled', color='black', alpha=.9)
        mask = (Mx2Hits.hit_time_slice[entry] == slice) & (Mx2Hits.hit_pe[entry]>min_pe)
        h = ax.hist(Mx2Hits.hit_time[entry][mask]/1000, bins=time_bin, log=True, histtype='stepfilled', color='red')
        max_bin = h0[0].max()

    else:
        h0 = ax.hist(Mx2Hits.hit_time[entry][mask]/1000, bins=time_bin, log=True, histtype='stepfilled', color='black', alpha=.9)
        max_bin = h0[0].max()
        for ts in range(1,Mx2Hits.n_slices[entry]+1):
            mask = (Mx2Hits.hit_time_slice[entry] == ts) & (Mx2Hits.hit_pe[entry]>min_pe)
            if (len(Mx2Hits.hit_time[entry][mask])>0):
                h = ax.hist(Mx2Hits.hit_time[entry][mask]/1000, bins=time_bin, log=True, histtype='stepfilled')
                max_bin = max(max_bin, h[0].max())
    
    ax.set_xlim(0,16)
    
    ax.set_xlabel("time [µs]", fontsize=15)
    ax.set_ylabel("hits", fontsize=15)
    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=15)
    
    # ax.set_xticks(20)
    
    ax.set_title("Mx2 timing", fontsize=15)
    ax.text(13, max_bin*12, f'Trigger number: {entry}', dict(size=15))
    ax.text(13, max_bin*2, f'Mx2 Slice number: {slice}', dict(size=15))     
    ax.text(0, max_bin*2, f'#Run: {0:05}', dict(size=15))      
    
    ax.grid(which='both', axis="y",linestyle='--', linewidth='0.5', color='black')
    ax.grid(which='major', axis="x",linestyle='--', linewidth='0.5', color='black')
    
    # ax.yaxis.get_label().set_fontsize(20)

# Plotting 2x2 hit time    

def plot_2x2_time(NdFlow, entry, ax):
    ax.clear()
    time = NdFlow.data.fields("ts_pps")[:]/10 - 1.2e6 * NdFlow.nd_trigger
    mask1 = (time<10e3) & (NdFlow.nd_trigger == entry)
    
    ax.hist(time[mask1] ,histtype="stepfilled",bins=500, color="red",alpha=.7)

    ax.set_title("2x2 time", y=.9, pad=-14, fontsize=15)
    ax.set_xlabel("time [µs]", fontsize=15, loc='right')
    ax.set_ylabel("hits", fontsize=15)
    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=15)
    
    ax.grid(which='both', axis="y",linestyle='--', linewidth='0.5', color='black')
    ax.grid(which='major', axis="x",linestyle='--', linewidth='0.5', color='black')
    ax.set_xlim(0,200)


    
#main plotter function
def view_event(Mx2Hits, NdFlow, trig, first_draw=True, mx2_min_pe=0, slice=0):
    
    start = time.time()
    plot_mx2_time(Mx2Hits, trig,axs[0], mx2_min_pe, slice)
    end = time.time()
    print(1, end - start)

    start = time.time()
    plot_2x2_time(NdFlow, trig,axs[1])
    # plt.hist(time ,histtype="stepfilled",bins=100000, log=True)
    end = time.time()
    print(2, end - start)

    start = time.time()
    plot_view_nd(Mx2Hits, NdFlow, trig,1,axs[2], first_draw, mx2_min_pe, slice)
    end = time.time()
    print(3, end - start)

    start = time.time()
    plot_view_nd(Mx2Hits, NdFlow, trig,2,axs[3], first_draw, mx2_min_pe, slice)
    end = time.time()
    print(4, end - start)

    start = time.time()
    plot_view_nd(Mx2Hits, NdFlow, trig,3,axs[4], first_draw, mx2_min_pe, slice)
    end = time.time()
    print(5, end - start)
    print()
    # f.savefig(f'plots/trigg_{trig}.png', facecolor='white')    
    plt.show()
    end = time.time()

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
    view_event(Mx2Hits, NdFlow, current_event, first_draw)

def on_previous_button_click(event):
    global current_event
    current_event = max(0, current_event-1)
    
    global current_slice 
    current_slice = 0

    global current_mx2_pe_cut
    current_mx2_pe_cut = 0

    first_draw = False
    view_event(Mx2Hits, NdFlow, current_event, first_draw)
    

def on_next_slice_click(event):
    global current_event
    first_draw = False

    global current_mx2_pe_cut
    current_mx2_pe_cut = 0

    global current_slice
    current_slice +=1
    if (current_slice >= Mx2Hits.n_slices[current_event]):
        current_slice = Mx2Hits.n_slices[current_event]

    view_event(Mx2Hits, NdFlow, current_event, first_draw, current_mx2_pe_cut, current_slice)


def on_previous_slice_click(event):
    global current_event

    first_draw = False

    global current_mx2_pe_cut
    current_mx2_pe_cut = 0

    global current_slice
    current_slice -=1

    if (current_slice <=0):
        current_slice = 0

    view_event(Mx2Hits, NdFlow, current_event, first_draw, current_mx2_pe_cut, current_slice)
    


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
    view_event(Mx2Hits, NdFlow, current_event, first_draw, current_mx2_pe_cut, current_slice)

    



Mx2Hits = Mx2Data("spill_SIM_MINERVA_00000_volume.root")
NdFlow = NdData('PicoRun4.1_1E17_RHC.flow.00000.FLOW.hdf5')
current_event = 0
current_slice = 0
f, axs = plt.subplots(5, 1, figsize=(18,10), gridspec_kw={'height_ratios' : [.5,.5,1,1,1]})

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






view_event(Mx2Hits, NdFlow, current_event)
