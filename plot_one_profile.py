import os, sys
sys.path.append('/Users/apple/Desktop/Winwork/Seismicity-Visualization-master')
import matplotlib.pyplot as plt
from reader import read_fault, read_ctlg, slice_ctlg
import numpy as np
from obspy.imaging.beachball import beach
import matplotlib.ticker as ticker

# config
fig_fname = 'jz_cluster_pro_1.pdf'
fctlg = 'Jz_dtmsms_catalog11.txt'
ffault = 'Tibet_Faults.txt'
j = 2 ##cluster number 
min_num = 50
max_dist = 3
lon_rng = [103.6, 104]
lat_rng = [33, 33.4]
x_ratio = np.cos(lat_rng[0] * np.pi / 180)
dep_rng = [0, 30]
prof_wd = 6 / 111
alpha = 0.6
fsize_label = 14
fsize_title = 16
fig_size = (10, 10)
mark_size = 8
mag_corr = 0.6
line_wid = 1.
annot_x = lon_rng[0] + (lon_rng[1]-lon_rng[0]) * 0.75
annot_y = lat_rng[0] + (lat_rng[1]-lat_rng[0]) * 0.75
idx_x = lon_rng[0]
idx_y = lat_rng[1]
prof_pnt  = np.array([[103.75,33.21], [103.80,33.21]])
#colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:olive','tab:cyan']
colors = ['blue','orange','green','red','purple','brown','olive','cyan','slateblue','lawngreen','darkred','orchid','black','sienna','salmon','pink','orchid','black','skyblue','lawngreen','fuchsia','lime','gold']
#colormap = plt.cm.Gist_ncar
#colors = [colormap(i) for i in np.linspace(0,0.9,len(ax.collections))]

num_colors = len(colors)
print(num_colors)

events = read_ctlg(fctlg)
events = slice_ctlg(events, lat_rng=lat_rng, lon_rng=lon_rng, dep_rng=dep_rng)
num_events = len(events)
lat = np.array(list(events['lat']))
lon = np.array(list(events['lon']))
dep = np.array(list(events['dep']))
mag = (np.array(list(events['mag'])) + mag_corr) * mark_size

# i/o paths
plt.style.use('seaborn-ticks')


def sort_clust(clusts):
    dtype=[('lat','O'),('clust','O')]
    clusters = []
    for clust in clusts:
        lat = np.median(clust[:,1])
        clusters.append((lat, clust))
    clusters = np.array(clusters,dtype=dtype)
    clusters = np.sort(clusters, order='lat')
    return clusters['clust']

vec_ab = prof_pnt[1] - prof_pnt[0]
vec_ab[0] *= x_ratio
abs_ab = np.linalg.norm(vec_ab)
plt.figure(figsize=(4,8))
#plt.figure()
ccs = [0.5]
for i,cc in enumerate(ccs):       
    plt.subplot(1,1,i+1)
    fclust = 'jz_clusters_cc%s.txt'%cc
    # read file
    f=open(fclust); lines=f.readlines(); f.close()
    clusts = []
    for line in lines:
      if line[0]=='#': clusts.append([])
      else: 
        lat, lon, dep, mag = [float(code) for code in line.split(',')[1:]]
        if lat>lat_rng[1] or lat<lat_rng[0]: continue
        if lon>lon_rng[1] or lon<lon_rng[0]: continue
        clusts[-1].append([lon, lat, (mag+mag_corr)*mark_size,dep]) ##add dep
    bg_clust = []
    for clust in clusts:
        if len(clust)<min_num and len(clust)>2: bg_clust += clust
    bg_clust = np.array(bg_clust)
    clusts = [np.array(clust) for clust in clusts if len(clust)>=min_num]
    clusts = sort_clust(clusts)
    faults = read_fault(ffault, lat_rng, lon_rng)
    events = read_ctlg(fctlg)
    events = slice_ctlg(events, lat_rng=lat_rng, lon_rng=lon_rng, dep_rng=dep_rng)
    lat = np.array(list(events['lat']))
    lon = np.array(list(events['lon']))
    dep = np.array(list(events['dep']))
    mag = (np.array(list(events['mag'])) + mag_corr) * mark_size

    # fill up edge
    edgex = [0,0,abs_ab*111,abs_ab*111]
    edgey = [dep_rng[0],dep_rng[1],dep_rng[0],dep_rng[1]]
    plt.scatter(edgex, edgey, alpha=0)
    # plot catalog & clusters
   
    prof_dist, prof_dep, prof_mag = [], [], []
    for ii in range(len(lat)):
            loc_c = np.array([lon[ii], lat[ii]])
            vec_ac = loc_c - prof_pnt[0]
            vec_ac[0] *= x_ratio
            abs_ac = np.linalg.norm(vec_ac)
            cos = vec_ac.dot(vec_ab) / abs_ab / abs_ac
            if abs_ac * (1-cos**2)**0.5 > prof_wd: continue
            if cos<0 or abs_ac*cos>abs_ab: continue
            prof_dist.append(abs_ac * cos * 111)
            prof_dep.append(dep[ii])
            prof_mag.append(mag[ii])
#    plt.scatter(prof_dist, prof_dep, prof_mag, alpha=alpha, color='gray')


    clust = clusts[j]
    # calc along profile dist
    num_events = len(clust)
    lon = clust[:,0]
    lat = clust[:,1]
    mag = clust[:,2]
    dep = clust[:,3]
    prof_dist, prof_dep, prof_mag = [], [], []
    for ii in range(num_events):
            loc_c = np.array([lon[ii], lat[ii]])
            vec_ac = loc_c - prof_pnt[0]
            vec_ac[0] *= x_ratio
            abs_ac = np.linalg.norm(vec_ac)
            cos = vec_ac.dot(vec_ab) / abs_ab / abs_ac
            if abs_ac * (1-cos**2)**0.5 > prof_wd: continue
            if cos<0 or abs_ac*cos>abs_ab: continue
            prof_dist.append(abs_ac * cos * 111)
            prof_dep.append(dep[ii])
            prof_mag.append(mag[ii])
          #  print(prof_dep)
    plt.scatter(prof_dist, prof_dep, prof_mag, alpha=alpha, color=colors[j%num_colors])
    for fault in faults: plt.plot(fault[:,0], fault[:,1], color='k', linewidth=line_wid)
    plt.annotate('$CC_{min} = %s$ \n$N_{min} = %s$ \n$\Delta$$L_{max} = %skm$'%(cc, min_num, max_dist),
             (annot_x,annot_y), fontsize=fsize_label, fontweight='bold')
     
    ax = plt.gca()
    ax.grid(True, linestyle='-.')
    plt.axis([0, 6, 0, 12])
    ax.invert_yaxis()
   # ax.add_collection(beach1)
    tick_spacing = 0.1
    #ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    #ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]   
    plt.setp(ax.xaxis.get_majorticklabels(), fontsize=fsize_label)
    plt.setp(ax.yaxis.get_majorticklabels(), fontsize=fsize_label)
   # if i>0: plt.setp(ax.get_yticklabels(), visible=False)
font2 = {'family' : 'Times New Roman',
'weight' : 'normal',
'size'   : 15,
}
plt.xlabel('Along-Profile Distance (km)',font2)
plt.ylabel('Depth (km)',font2)

plt.tight_layout(True)
plt.savefig(fig_fname)
#plt.show()

