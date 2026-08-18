Particle image velocimetry (PIV) is widely used in experimental fluid dynamics, where high particle seeding density is critical for improving spatial resolution. PET-PIV is a monocular 3D PIV framework that combines compact PSF engineering (coded aperture/wavefront shaping) for depth encoding with deep learning-based flow reconstruction.

Our contributions are twofold:
1) We extend the density tolerance of wavefront-shaping-based monocular 3D particle tracking velocimetry
2) To further overcome the density limitation, we establish, for the first time, Eulerian field reconstruction in monocular PIV, which directly recovers 3D velocity fields from temporal image intensity variations.
For more technical details, please refer to https://arxiv.org/abs/2607.29163 

# env

# tracking mode 
1, run collect_psfs.ipynb and obtain a PSF collection for in-situ modeling
2, run in_situ_psf.ipynb

train_locnet.ipynb

infer_locnet.ipynb

link.ipynb

alignment.ipynb

CCM/CCM_v2_J.m

# field mode

