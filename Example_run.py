#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 11:56:10 2026

@author: Alexandria McPherson

Python implementation of multi-SSS preprocessing as in https://pubmed.ncbi.nlm.nih.gov/40541227/
Combined Python implementation with Fosters Inverse with mSSS as in https://www.biorxiv.org/content/10.64898/2026.07.21.739883v1

"""

import numpy as np
import os
import mne

## import Foster's inverse with mSSS, and mSSS functions
from fit_spheres_to_mri import fit_spheres_to_mri
from run_fosters_mSSS import apply_preprocessing

########################################################################
## Example Run ####
if __name__ == '__main__':
    #################################################
    ### 1. DEFINE PATHS
    ## define raw file directory and raw file names
    sample_dir ='~/path/to/raw/data'
    sample_files = ['~raw.fif']
    n_subjects = len(sample_files)
    trans_files=['~raw_trans.fif']
                 
    ## define subject MRI files
    subjects_dir = '~/freesurfer/subjects/'
    
    #################################################
    ## define constants
    trigger_chan = 'di2' 
    # epoch duration
    tmin = -0.05
    tmax = 0.3
    ##specify plotting args
    ts_args = ts_args = dict(
        time_unit="s",
        #ylim=dict(mag=(-400, 400)),
        gfp=True
    )
    topomap_args = dict(
        time_unit="s",
        #vlim=(-400,400)
    )
    baseline=None
    
    ## --- generate results per subject -------------------------------------------
    freq_min = 0.5
    freq_max = 80
    
    for file,trans in zip(sample_files,trans_files):
        ## read subject from file name
        subject = file[4:8]
        
        trans_path = os.path.join(sample_dir, trans)
        
        ###################################
        # 3. load data and find events
        raw = mne.io.read_raw_fif(os.path.join(sample_dir, file),preload=False, allow_maxshield='no')
        events = mne.find_events(raw, stim_channel=trigger_chan, shortest_event=1)

        ###################################
        # 4. Prepare data for preprocessing
        ## drop bad channels -- to make sure dimensions of data match info later on
        bad_indices = [raw.ch_names.index(ch) for ch in raw.info['bads']]
        bads = raw.info["bads"]
        raw.drop_channels(bads)
        raw.drop_channels(['di2'])
    
        ## high and low - pass, notch filter raw data
        raw.load_data().filter(l_freq=freq_min, h_freq=freq_max)
        meg_picks = mne.pick_types(raw.info, meg=True)
        raw.notch_filter(freqs=60, picks=meg_picks)
       
        ###################################
        # 5. Calculate expansion centers for mSSS. These are in HEAD COORDS
        conductivity = conductivity = (0.3, 0.006, 0.3)
        bem_model = mne.make_bem_model(subject=subject, ico=4, conductivity=conductivity, subjects_dir=subjects_dir)
        n_centers = 2
        centers = fit_spheres_to_mri(subjects_dir, 
                                     subject, 
                                     bem_model, 
                                     trans_path, 
                                     n_centers, 
                                     show_spheres=True)

        ###################################
        # 6. Run preprocessing
        # do_fos=TRUE and do_msss=TRUE will execute Fosters Inverse with mSSS and empirical N 
        # do_fos=FALSE and do_msss=TRUE will execute mSSS
        # do_fos=TRUE and do_msss=FALSE will execute Fosters Inverse with SSS and empirical N 
        ######################3
        do_fos=True
        do_msss=True
        Lin=8
        Lout=3
        ch_types=np.ones(raw.info["nchan"]) # 1's for magnetometers, 0's for gradiometers
        raw_msss = apply_preprocessing(np.transpose(centers[0]), 
                                   np.transpose(centers[1]), 
                                   raw, 
                                   do_fos, 
                                   ch_types, 
                                   Lin, Lout)
        
        ###################################
        # 7. make epochs, visualize evokeds
        epochs_msss = mne.Epochs(raw_msss, events, tmin=tmin, tmax=tmax, baseline=None, preload=True)
        evoked_msss = epochs_msss.average()
        fig = evoked_msss.plot_joint(ts_args=ts_args, topomap_args=topomap_args, title=subject+"Foster's with mSSS")






