Particle image velocimtry (PIV) is widely used in experimental fluid dynamics and high particle seeding density enables improved spatial resolution and flow characterization. PET-PIV is a moncular 3D particle image velocimetry framework that applies compact PSF engineering/coded aperture/wavefront shaping for depth encoding and deep learning for flow reconstruction. 
We contribute in two aspects:
1) extend the density tolerance of wavefront-shaping-based monocular 3D particle tracking velocimty
2) establishe, for the first time, Eulerian field reconstruction in monocular PIV that directly reconstructs velocity fields from image intensity variations.
For more technical details, please refer to https://arxiv.org/abs/2607.29163 

# env

## tracking mode 
1, run collect_psfs.ipynb and obtain a PSF collection for in-situ modeling
2, run in_situ_psf.ipynb

train_locnet.ipynb

infer_locnet.ipynb

link.ipynb

alignment.ipynb

CCM/CCM_v2_J.m

## field mode

