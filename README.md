Particle image velocimetry (PIV) is widely used in experimental fluid dynamics, where high particle seeding density is critical for improving spatial resolution. PET-PIV is a monocular 3D PIV framework that combines compact PSF engineering (coded aperture/wavefront shaping) for depth encoding with deep learning-based flow reconstruction.

Our contributions are twofold:
1) We extend the density tolerance of wavefront-shaping-based monocular 3D particle tracking velocimetry
2) To further overcome the density limitation, we establish, for the first time, Eulerian field reconstruction in monocular PIV, which directly recovers 3D velocity fields from temporal image intensity variations.

For more technical details, please refer to https://arxiv.org/abs/2607.29163 

# env
The code is organized as a Jupyter Notebook/JupyterLab project. To run the code, a virtual Python environment with Jupyter Notebook/Lab and a CUDA-enabled version of PyTorch is required. Other required Python packages can be conveniently installed as needed based on the error messages encountered during execution.

# tracking mode
0, switch the working directory to the tracking folder

1, run collect_psfs.ipynb which processes sparse data (images with sparse particles) and outputs a collection of PSFs

2, run in_situ_psf.ipynb

3, train_locnet.ipynb

4, infer_locnet.ipynb

5, link.ipynb

6, alignment.ipynb

7, CCM/CCM_v2_J.m

# field mode

