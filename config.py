""" Configure file for seismicity visualization
"""
import os
import numpy as np

# 1. plot event location
class Config_Loc(object):
  def __init__(self):

    # catalog info
    self.ctlg_path = '/home/zhouyj/Xiaojiang/run_pad/hypoinverse/output/zsy.csv'
    self.lon_rng = [102.5, 103.9] #[102.25, 103.8]
    self.lat_rng = [24.3, 26.5] #[25.9, 27.5]
    self.dep_rng = [0, 30.]
    self.prof_pnt = np.array([[103,25],[103.5,26]]) # ref points for profile
    self.prof_wd = 10. # width of profile km
    # plot params
    self.fig_title = 'PAD HypoInverse: ZSY Network'
    self.fig_fname = 'zsy_pad_hyp.pdf'
    self.fsize_label = 14
    self.fsize_title = 16
    self.mark_size = 3.
    self.alpha = 0.6
    self.fig_xsize = 5
    self.cmap = 'hot'
    self.cbar_pos = [0.76,0.6,0.03,0.25]
    self.cbar_ticks = np.arange(0,1.1,0.25)


# 2. plot FMD & M-t
class Config_Mag(object):
  def __init__(self):

    # catalog info
    self.ctlg_path = '/home/zhouyj/Xiaojiang/run_pad/hypoinverse/output/zsy.csv'
    self.mag_rng = [-1, 8.]
    # plot params
    self.fig_title_fmd = 'FMD: ZSY Network'
    self.fig_title_mt = 'M-t: ZSY Network'
    self.fsize_label = 14
    self.fsize_title = 16
    self.mark_size = 10.
    self.alpha = 0.6


