# Fosters_inverse_mSSS

These functions and Example implement the mSSS method as described in [Refined signal space separation methods for on-scalp MEG systems](https://pubmed.ncbi.nlm.nih.gov/40541227/) McPherson et al. 2025 from MATLAB into Python. See linked manuscript and this [Repository](https://github.com/xannnimal/refined_SSS_methods_onscalp_MEG) for original MATLAB implementation. 

Next, these functions run the preprocessing method Foster's Inverse using the mSSS basis as derived and verified in [Noise Optimization of Basic Signal Component Extraction for Cryogenic and On-Scalp Magnetoencephalography (MEG)](https://www.biorxiv.org/content/10.64898/2026.07.21.739883v1) McPherson et al. July 2026. 

Foster's Inverse with mSSS is a novel MEG preprocessing method modified specifically for on-scalp MEG, such as OPM-MEG, data preprocessing, where the internal brain activity is isolated using multiple overlapping and optimized vector spherical harmonic basis sets (mSSS) such that the whole brain is spanned without encroaching on the sensors, and estimated utilizing information about the sensor noise profile and artifacts (Foster's Inverse).

## Description of Files
*`fit_spheres_to_mri.py` takes a BEM model and subject-specific MRI information to dynamically fit two spherical basis sets to span the brain, returns the origins of these two optimized expansions.

*`run_fosters_msss.py` uses the two origins to preprocess and clean the raw MEG data using Foster's Inverse with mSSS. These functions can also be used to preprocess with mSSS only. This implementation of Foster's Inverse automatically used the Empirical method for estimating noise covariance included with MNE-Python. For more details, alternative implementations, and Foster's Inverse with SSS, see this [Repository](https://github.com/xannnimal/fosters_inverse_sss/tree/main)

*`Example_run.py` shows how to implement the above to functions to preprocess MEG data

## Dependencies
* MNE Python

To visualize brain coverage of the optimized basis expansions:
* matplotlib
* pyvista
* pyvistaqt

Must import the required functions from this repository

 ```bash
from fit_spheres_to_mri import fit_spheres_to_mri
from run_fosters_mSSS import apply_multi_sss
