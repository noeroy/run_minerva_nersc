import uproot
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import os


    
class Mx2Data:
    def __init__(self, filename):
        self.file = uproot.open(filename)

        self.n_idhits = self.file["minerva"]["n_idhits"].array(library="np")
        self.hits_id_per_mod = self.file["minerva"]["hits_id_per_mod"].array(library="np")
        self.hit_strip = self.file["minerva"]["hit_strip"].array(library="np")
        self.hit_plane = self.file["minerva"]["hit_plane"].array(library="np")
        self.hit_module = self.file["minerva"]["hit_module"].array(library="np")
        self.hit_view = self.file["minerva"]["hit_view"].array(library="np")
        self.n_odhits = self.file["minerva"]["n_odhits"].array(library="np")
        self.hits_od_per_mod = self.file["minerva"]["hits_od_per_mod"].array(library="np")
        self.hit_bar = self.file["minerva"]["hit_bar"].array(library="np")
        self.hit_pe = self.file["minerva"]["hit_pe"].array(library="np")
        self.hit_time = self.file["minerva"]["hit_time"].array(library="np")
        self.hit_time_slice = self.file["minerva"]["hit_time_slice"].array(library="np")
        self.hit_norm_energy = self.file["minerva"]["hit_norm_energy"].array(library="np")

        self.offsetX = self.file["minerva"]["offsetX"].array(library="np")
        self.offsetY = self.file["minerva"]["offsetY"].array(library="np")
        self.offsetZ = self.file["minerva"]["offsetZ"].array(library="np")

        self.ev_gps_time_sec = self.file["minerva"]["ev_gps_time_sec"].array(library="np")
        self.ev_gps_time_usec = self.file["minerva"]["ev_gps_time_usec"].array(library="np")

        self.n_slices = self.file["minerva"]["n_slices"].array(library="np")

        self.clus_id_coord = self.file["minerva"]["clus_id_coord"].array(library="np")
        self.clus_id_z = self.file["minerva"]["clus_id_z"].array(library="np")
        self.clus_id_view = self.file["minerva"]["clus_id_view"].array(library="np")
        self.clus_id_pe = self.file["minerva"]["clus_id_pe"].array(library="np")
        self.clus_id_time_slice = self.file["minerva"]["clus_id_time_slice"].array(library="np")
        self.clus_id_type = self.file["minerva"]["clus_id_type"].array(library="np")

        self.minerva_time = (self.ev_gps_time_sec - self.ev_gps_time_sec[0])*1e6+(self.ev_gps_time_usec)
        self.minerva_trigger = np.array((self.minerva_time/1.2e6)).astype(int)

        self.hits_id_per_mod = self.file["minerva"]["hits_id_per_mod"].array(library="np")
        
        self.n_tracks = self.file["minerva"]["n_tracks"].array(library="np") 
        self.n_blobs_id = self.file["minerva"]["n_blobs_id"].array(library="np")
        self.trk_vis_energy = self.file["minerva"]["trk_vis_energy"].array(library="np")
        self.trk_time_slice = self.file["minerva"]["trk_time_slice"].array(library="np")


tot_bad_crossing = 0
tot_good_crossing = 0
path = "/pscratch/sd/m/mkramer/output/MegaRun5/run-minerva/MegaRun5_1E20_RHC.minerva/DST/0000000/"
i_f = 0
tot_bad_crossing = 0
tot_good_crossing = 0
max_cl = 40

tot_bad_crossing_ar = np.zeros(max_cl-1)
tot_good_crossing_ar = np.zeros(max_cl-1)
tot_good_crossing_ar2 = np.zeros(max_cl-1)
tot_good_crossing_ar1 = np.zeros(max_cl-1)

tot_bad_crossing_ar_trackable = np.zeros(max_cl-1)
tot_good_crossing_ar_trackable = np.zeros(max_cl-1)
tot_good_crossing_ar_trackable1 = np.zeros(max_cl-1)
tot_good_crossing_ar_trackable2 = np.zeros(max_cl-1)

for filename in os.listdir(path):
    if (i_f%100 == 0):
      print(i_f)
    i_f+=1
    if (i_f > 200):
      break
    file = path+filename
    Mx2Hits = Mx2Data(file)
    n_good_crossing = 0
    n_bad_crossing = 0
    for ievent in range(len(Mx2Hits.clus_id_time_slice)):
        clus_id_time_slice = (Mx2Hits.clus_id_time_slice[ievent])
        trk_time_slice = (Mx2Hits.trk_time_slice[ievent])
        clus_id_z = (Mx2Hits.clus_id_z[ievent])
        nslices = Mx2Hits.n_slices[ievent]
        clus_id_type = Mx2Hits.clus_id_type[ievent]

        for islice in range(1, nslices):
            clus_mask_trackable = (clus_id_time_slice == islice) & (clus_id_type == 1)
            clus_mask = (clus_id_time_slice == islice)
            trk_mask = (trk_time_slice == islice)
            us_sum = (clus_id_z[clus_mask_trackable] < 6000).sum() 
            ds_sum = (clus_id_z[clus_mask_trackable] > 6000).sum() 
            for ncl in range(max_cl-1):
                if (us_sum > ncl+1 and ds_sum > ncl+1):
                    if( trk_mask.sum() ==0):
                        tot_bad_crossing_ar_trackable[ncl] +=1
                    else:
                        tot_good_crossing_ar_trackable[ncl] +=1
                    if (trk_mask.sum() >=2):
                        tot_good_crossing_ar_trackable2[ncl] +=1
                    if (trk_mask.sum() ==1):
                        tot_good_crossing_ar_trackable1[ncl] +=1                        
            us_sum = (clus_id_z[clus_mask] < 6000).sum() 
            ds_sum = (clus_id_z[clus_mask] > 6000).sum() 
            for ncl in range(max_cl-1):
                if (us_sum > ncl+1 and ds_sum > ncl+1):
                    if( trk_mask.sum() ==0):
                        tot_bad_crossing_ar[ncl] +=1
                    else:
                        tot_good_crossing_ar[ncl] +=1
                    if (trk_mask.sum() >=2):
                        tot_good_crossing_ar2[ncl] +=1
                    if (trk_mask.sum() ==1):
                        tot_good_crossing_ar1[ncl] +=1    



trk_ratio = tot_good_crossing_ar / (tot_good_crossing_ar+tot_bad_crossing_ar)
trk_ratio_trackable = tot_good_crossing_ar_trackable / (tot_good_crossing_ar_trackable+tot_bad_crossing_ar_trackable)
nclus_array = np.arange(1,max_cl)
f1,axs = plt.subplots(2,1,figsize=(10,6))
axs[0].plot(nclus_array, trk_ratio, label="All clusters")
axs[0].plot(nclus_array, trk_ratio_trackable, label="Trackable clusters")
axs[0].set_xlabel("#of clusters upstream and downstream")
axs[0].set_ylabel("Tracking efficiency")
axs[0].legend()


axs[1].plot(nclus_array, tot_good_crossing_ar_trackable1/tot_good_crossing_ar_trackable, label="1 track in the slice")
axs[1].plot(nclus_array, tot_good_crossing_ar_trackable2/tot_good_crossing_ar_trackable, label=">1 tracks in the slice")
axs[1].set_xlabel("#of clusters upstream and downstream")
axs[1].set_ylabel("proportion")
axs[1].legend()
plt.savefig("tracking_eff_tot.pdf")
plt.savefig("tracking_eff_tot.png")
