import uproot
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
import subprocess
import sys


class Mx2Data:
    def __init__(self, filename):
        self.file = uproot.open(filename)

        self.n_idhits = self.file["minerva"]["n_idhits"].array()
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
        
        self.n_tracks = self.file["minerva"]["n_tracks"].array() 
        self.n_blobs_id = self.file["minerva"]["n_blobs_id"].array()
        self.trk_vis_energy = self.file["minerva"]["trk_vis_energy"].array() 

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


def draw_view(my_view=0):
    figs = []
    fig, axs = plt.subplots(1, 4, figsize=(15,5),  gridspec_kw={'width_ratios' : np.array([12,10,10,12])/12}, dpi=100)

    view_title = ["ALL views", "X view", "U view", "V view"]
    mod = np.array(np.concatenate(Mx2Hits.hit_module, axis=0))
    strip = np.array(np.concatenate(Mx2Hits.hit_strip, axis=0))
    pe = np.array(np.concatenate(Mx2Hits.hit_pe, axis=0) / len(Mx2Hits.hit_pe))
    view = np.array(np.concatenate(Mx2Hits.hit_view, axis=0)) 
    
    us = (mod>0) & (mod < 13) & (view==my_view)
    ds = (mod>=13) & (mod < 23) & (view==my_view)
    ecal = (mod>=23) & (mod < 33) & (view==my_view)
    hcal = (mod>=33)  & (view==my_view)

    
    if (my_view ==0):
        us = (mod>0) & (mod < 13) 
        ds = (mod>=13) & (mod < 23) 
        ecal = (mod>=23) & (mod < 33) 
        hcal = (mod>=33) 
        
    
    j=0
    axs[j].hist2d(mod[us], strip[us], weights=pe[us], bins=[np.arange(.5,13.5,1),np.arange(.5,128.5,1)])
    axs[j].set_title("Upstream")
    axs[j].set_ylabel("strip")
    axs[j].set_xlabel("module")

    j+=1
    axs[j].hist2d(mod[ds], strip[ds], weights=pe[ds], bins=[np.arange(12.5,23.5,1),np.arange(.5,128.5,1)])
    axs[j].set_title("Downstream")
    axs[j].set_xlabel("module")

    j+=1
    axs[j].hist2d(mod[ecal], strip[ecal], weights=pe[ecal], bins=[np.arange(22.5,33.5,1),np.arange(.5,128.5,1)])
    axs[j].set_title("ECal")
    axs[j].set_xlabel("module")

    j+=1
    hhcal = axs[j].hist2d(mod[hcal], strip[hcal], weights=pe[hcal], bins=[np.arange(32.5,45.5,1),np.arange(.5,128.5,1)])
    axs[j].set_title("HCal")
    axs[j].set_xlabel("module")
    axs[j].set_xticks(np.arange(34,46, 2))
    cl = plt.colorbar(hhcal[3], ax= axs[j])
    cl.set_label('Average pe/trigger')
    
    axs[0].text(0, 140, view_title[my_view], dict(size=15))
    # plt.show()
    figs.append(fig)
    return figs

def draw_time():
    figs = []
    fig, ax  = plt.subplots(1,1, figsize=(15,3))

    time = np.concatenate(Mx2Hits.hit_time, axis=0)/1000
    time_bin = np.linspace(0,16,1600)

    hist, edges, _  = ax.hist(time[(time>0)], bins=time_bin, histtype='stepfilled', log=True)
    ax.clear()
    scaled_hist = hist / len(Mx2Hits.hit_time)
    
    ax.hist(edges[:-1], weights=scaled_hist,bins=time_bin, histtype='stepfilled', log=True)

    # ax.step(edges[:-1], scaled_hist)
    # ax.set_yscale('log')
    # ax.fill_between(edges[:-1], scaled_hist)  # Filled steps


    ax.set_xlim(0,16)

    ax.set_xlabel("time [µs]", fontsize=15)
    ax.xaxis.set_label_coords(0.5, -0.1)
    ax.set_ylabel("hits", fontsize=15)
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)

    # ax.set_xticks(20)

    ax.set_title("Average hit time distribution", fontsize=20)
    # ax.text(13, max_bin*.7, f'Trigger number: {entry}', dict(size=15))
    
    ax.text(0, scaled_hist.max()*2, f'#Run: {0:05}', dict(size=15))

    ax.grid(which='both', axis="y",linestyle='--', linewidth='0.5', color='black')
    ax.grid(which='major', axis="x",linestyle='--', linewidth='0.5', color='black')
    # plt.show()
    figs.append(fig)
    return figs


def draw_ntracks():
    nbins=50
    x_bins = np.array(Mx2Hits.minerva_time/1e6)
    ntracks = np.array(Mx2Hits.n_tracks)
    
    nshowers = np.array(Mx2Hits.n_blobs_id)
    
    figs = []
    figs.append(draw_1D_hist(x_bins, nbins, ntracks, "Time [s]", "# reco tracks"))
    
    fig, ax = plt.subplots(1,1)
    fig.set_facecolor('white')
    labels = ["tracks","shower"]
    ax.pie([ntracks.sum(), nshowers.sum()], labels=labels,autopct='%1.1f%%')
    figs.append(fig)
    return figs

def draw_nhits():
    figs=[]
    nbins=50
    x_bins = np.array(Mx2Hits.minerva_time/1e6)
    weights = np.array(Mx2Hits.n_idhits)
    f = draw_1D_hist(x_bins, nbins, weights, "Time [s]", "# reco digits")
    figs.append(f)
    return figs

def draw_pe_time():  
    figs = []
    hpe_trig = []
    hE_trig = []
    for h in Mx2Hits.hit_pe:
        h2 = np.array(h)
        hpe_trig.append(h2[h2>0].sum())

    for h2 in Mx2Hits.hit_norm_energy:
        h2 = np.array(h2)
        hE_trig.append(h2[h2>0].sum())

    x_bins = np.array(Mx2Hits.minerva_time/1e6)
    hpe_trig = np.array(hpe_trig)
    hE_trig = np.array(hE_trig)
    
    nbins=50
    figs.append(draw_1D_hist(x_bins, nbins, hpe_trig, "Time [s]", "pe"))
    figs.append(draw_1D_hist(x_bins, nbins, hE_trig, "Time [s]", "Reco energy [MeV]"))
    
    mod = np.array(np.concatenate(Mx2Hits.hit_module, axis=0))
    
    av_energy = np.array(np.concatenate(Mx2Hits.hit_norm_energy, axis=0)) / len(Mx2Hits.hit_module)
    
    figs.append(draw_1D_hist(mod[mod>0], 44, av_energy[mod>0], "Module", "Average Energy [MeV]"))
    return figs 

def draw_1D_hist(x_bins, nbins, wheights, xtitle, ytitle):
    fig, ax  = plt.subplots(1,1, figsize=(15,3.5))
    hist, edges, _ = ax.hist(x_bins, weights=wheights, bins=nbins, histtype="step")
    ax.clear()
    ax.hist(edges[:-1], weights=hist, histtype="step", color='C0',  bins=nbins, align='mid')
    ax.hist(edges[:-1], weights=hist, histtype="stepfilled", color='C0', alpha=.7, bins=nbins, align='mid')
    
    
    ax.set_xlabel(xtitle,fontsize=15)
    ax.set_ylabel(ytitle, fontsize=15)
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)
    ax.grid(which='major', axis="x",linestyle='--', linewidth='0.5', color='black')
    ax.grid(which='major', axis="y",linestyle='--', linewidth='0.5', color='black')
    
    ax.set_xlim(edges[0], edges[-1])
    return(fig)


if len(sys.argv) != 2:
        print("Usage: python script_name.py <integer_argument>")
        sys.exit(1)
try:
    # Attempt to convert the command-line argument to an integer
    run = int(sys.argv[1])
except ValueError:
    print("Error: The argument must be an integer.")
    sys.exit(1)

Mx2Hits = Mx2Data(f"output_dir/run_{run:05}/dst/MiniRun5_1E19_RHC.minerva.{run:05}.dst.root")


figs = []
figs.append(draw_time())
figs.append(draw_nhits())
figs.append(draw_ntracks())
figs.append(draw_view())
figs.append(draw_view(1))
figs.append(draw_view(2))
figs.append(draw_view(3))

# # Save each figure as a separate PNG file, merge them in a pdf and remove the initial pngs.

directory_path=f'output_dir/plots/validation_plots/'

if not os.path.exists(directory_path):
    os.makedirs(directory_path)

i=0
png_filenames = []
for f_array in figs:
    for fig in f_array:
        png_filename = os.path.join(f'output_dir/plots/validation_plots/figure_{i + 1}.png')
        fig.savefig(png_filename, format='png', dpi=60, bbox_inches='tight')
        png_filenames.append(png_filename)
        i+=1

subprocess.run(['convert', 'output_dir/plots/validation_plots/figure_*', f'output_dir/plots/validation_plots/mx2_run_{run:05}.pdf'])
for png_filename in png_filenames:
    os.remove(png_filename)


