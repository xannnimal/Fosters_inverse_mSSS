# Fosters_inverse_mSSS

These functions and Example implement the mSSS method as described in [Refined signal space separation methods for on-scalp MEG systems](https://pubmed.ncbi.nlm.nih.gov/40541227/) McPherson et al. 2025 from MATLAB into Python. See linked manuscript and this [Repository](https://github.com/xannnimal/refined_SSS_methods_onscalp_MEG) for original mSSS MATLAB implementation. All other functions for calculating the vector spherical harmonic expansions are from [Presentation of electromagnetic multichannel data: The signal space separation method](https://pubs.aip.org/aip/jap/article-abstract/97/12/124905/893620/Presentation-of-electromagnetic-multichannel-data?redirectedFrom=fulltext) Taulu and Kajola, 2005.

Next, these functions run the preprocessing method Foster's Inverse using either the SSS or the mSSS basis as derived and verified in [Noise Optimization of Basic Signal Component Extraction for Cryogenic and On-Scalp Magnetoencephalography (MEG)](https://www.biorxiv.org/content/10.64898/2026.07.21.739883v1) McPherson et al. July 2026. 

Foster's Inverse with SSS and/or mSSS is a novel MEG preprocessing method modified specifically for on-scalp MEG, such as OPM-MEG, data preprocessing, where the internal brain activity is isolated using multiple overlapping and optimized vector spherical harmonic basis sets (mSSS) such that the whole brain is spanned without encroaching on the sensors, and estimated utilizing information about the sensor noise profile and artifacts (Foster's Inverse).

User can choose to preprocesses their MEG data with Foster's Inverse with SSS, Foster's Inverse with mSSS, or mSSS. See instructions below in description of files!

## Description of Files
1. `fit_spheres_to_mri.py` takes a BEM model and subject-specific MRI information to dynamically fit two spherical basis sets to span the brain, returns the origins of these two optimized expansions.

2. `run_fosters_msss.py` contains the function to run different SSS-based preprocessing methods depending on two user inputs, which helper functions

## How to call the preprocessing function
1. Foster's Inverse with SSS
- `raw_fos_sss = apply_preprocessing(raw, Nt=[Ntmin,Ntmax])` with a user specific baseline period of data between [Ntmin,Ntmax] for the noise covariance. See `Example_run.py` for an example use case
2. Foster's Inverse with mSSS
- `raw_fos_msss = apply_preprocessing(raw, Nt=[Ntmin,Ntmax],origin=centers)` with a user specific baseline period of data between [Ntmin,Ntmax] for the noise covariance, and takes in the two origins calculated by `fit_spheres_to_mri.py` exactly as returned as matrix `centers`. See `Example_run.py` for an example use case
3. mSSS 
- `raw_msss = apply_preprocessing(raw, Nt=None,origin=centers)` will run mSSS using the two origins calculated by `fit_spheres_to_mri.py` exactly as returned as matrix `centers`. See `Example_run.py` for an example use case
       
This implementation of Foster's Inverse automatically used the Empirical method for estimating noise covariance included with MNE-Python. For more details and alternative noise covariance methods, see this [Repository](https://github.com/xannnimal/fosters_inverse_sss/tree/main)

3. `Example_run.py` shows how to implement the above to functions to preprocess MEG data

## Dependencies
* MNE Python

To visualize brain coverage of the optimized basis expansions:
* matplotlib
* pyvista
* pyvistaqt

Must import the required functions from this repository

 ```bash
from fit_spheres_to_mri import fit_spheres_to_mri
from run_fosters_mSSS import apply_preprocessing
