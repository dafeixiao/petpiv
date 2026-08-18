# PET-PIV
Particle image velocimetry (PIV) is widely used in experimental fluid dynamics, where high particle seeding density is critical for improving spatial resolution. PET-PIV is a monocular 3D PIV framework that combines compact PSF engineering (coded aperture/wavefront shaping) for depth encoding with deep learning-based flow reconstruction.

Our contributions are twofold:
1) We extend the density tolerance of wavefront-shaping-based monocular 3D particle tracking velocimetry
2) To further overcome the density limitation, we establish, for the first time, Eulerian field reconstruction in monocular PIV, which directly recovers 3D velocity fields from temporal image intensity variations.

For more technical details, please refer to https://arxiv.org/abs/2607.29163 

# environment
The code is organized as a Jupyter Notebook/JupyterLab project. 

To run the code, a virtual Python environment with Jupyter Notebook/Lab and a CUDA-enabled version of PyTorch is required. Other required Python packages can be conveniently installed as needed based on the error messages encountered during execution.

# tracking mode
1, run **collect_psfs.ipynb**. It processes sparse data (images with sparse particles) and outputs a collection of PSFs.

2, run **in_situ_psf.ipynb**. It takes the PSF collection and optical parameters, and outputs a PSF model file.

3, run **train_locnet.ipynb**. It employs the calibrated model and generates syhthetic training data (images-x, particle localizations-y).

4, run **infer_locnet.ipynb**. It applies the trained network to localize particles in experimental data. The ouput is a localizaiton list file.

5, run **alignment.ipynb**. It builds the registration/alignment between PET-PIV localizaitons and shake-the-box localizaitons (4-camera references). In trackign mode, this is optional and can be used for validation. In field mode, this is mandatory to generate high-quality training data.

6, run **link.ipynb**. It takes the (aligned) localizaiton list file and generate particle trajectories using a four-frame tracking algorithm (adapted from https://github.com/ronshnapp/MyPTV). The output includes velocities.

7, run **CCM/CCM_v2_J.m**. It takes the tracked velocities and generates volumetric velocity fields.

# field mode

1, run **training_data.ipynb**. It takes PET-PIV data, TomoPIV data, and the alignment file, and outputs paired training data (cropped PET-PIV image sequence with algined velocity components U, V, W).

2, run **train_velnet.ipynb**. This trains a VelNet for directly reconstructing velocitry fields from image intensity variations.

3, run **infer_velnet.ipynb**. The application of the trained VelNet to experimental PET-PIV data.


Please contact dafeixiao@gmail.com if you have any questions.




