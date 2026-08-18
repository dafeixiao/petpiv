
import numpy as np
import torch
import torch.nn as nn
import torch.fft as fft
import torch.nn.functional as F


class LensModelBase(nn.Module):  # basic function of a lens model: aperture, standard PSF,
    def __init__(self, param_dict):
        super().__init__()
        device = param_dict['device']
        N, N_aper = param_dict['N'], param_dict['N_aper']
        scaling_factor, rotation_factor = param_dict['scaling_factor'], param_dict['rotation_factor']
        coef_defocus = param_dict['coef_defocus']
        im_size = param_dict['im_size']
        sigma0 = param_dict['sigma0']

        self.sigma1 = param_dict['sigma1']
        self.N = N

        self.device = device
        self.idx05 = int(np.floor(N / 2))
        self.im_size = im_size
        self.r0, self.c0 = int(np.round((N - self.im_size[0]) / 2)), int(np.round((N - self.im_size[1]) / 2))

        xl = torch.linspace(-1, 1, N, device=device)
        xx, yy = torch.meshgrid(xl, xl, indexing='xy')
        rr2 = xx ** 2 + yy ** 2
        rr = torch.sqrt(rr2)
        aper = (rr < (N_aper / N)).type(torch.float32)
        self.aper = aper
        self.defocus_base = rr2 * coef_defocus

        if 'phase_mask' in param_dict:
            self.phase_mask = torch.tensor(param_dict['phase_mask'], device=device)
        else:
            self.phase_mask = torch.clone(aper)

        # gaussian blur
        g_size = 21
        g_r = int(g_size / 2)
        g_xs = torch.linspace(-g_r, g_r, g_size, device=device)
        self.g_xx, self.g_yy = torch.meshgrid([g_xs, g_xs], indexing='xy')
        self.sigma0 = sigma0
        blur_kernel0 = torch.exp(-0.5 * (self.g_xx ** 2 + self.g_yy ** 2) / sigma0 ** 2)
        self.blur_kernel0 = blur_kernel0 / blur_kernel0.sum()  # gaussian kernel before scaling

        self.scaling_factor = torch.tensor([[scaling_factor]], device=device)
        self.rotation_factor = torch.tensor([[rotation_factor]], device=device)
        self.translation_factor = N / 2
        self.tensor0 = torch.zeros((1, 1), device=device)
        self.tensor001 = torch.tensor([[0, 0, 1]], device=device)
        self.tensor_eye3 = torch.eye(3, device=device)

    def scaling(self):
        scaling_matrix = torch.cat((torch.cat((self.scaling_factor.view(1, 1), self.tensor0, self.tensor0), dim=1),
                                    torch.cat((self.tensor0, self.scaling_factor.view(1, 1), self.tensor0), dim=1),
                                    self.tensor001),
                                   dim=0)
        return scaling_matrix

    def translation(self, xys):
        translation_matrix = self.tensor_eye3.repeat(xys.shape[0], 1, 1)
        translation_matrix[:, :2, 2] = xys / self.translation_factor * (-1)
        return translation_matrix

    def rotation(self):
        rotation_matrix = torch.cat(
            (torch.cat(
                (torch.cos(self.rotation_factor.view(1, 1)), -torch.sin(self.rotation_factor.view(1, 1)), self.tensor0),
                dim=1),
             torch.cat(
                 (torch.sin(self.rotation_factor.view(1, 1)), torch.cos(self.rotation_factor.view(1, 1)), self.tensor0),
                 dim=1),
             self.tensor001), dim=0)
        return rotation_matrix

    def get_psfs(self, xyzps):  # one time of blur operation
        """
        xyzps: tensor, rank 2, N by 4
        """
        Nemitters = xyzps.shape[0]
        xys, zs, photons = xyzps[:, :2], xyzps[:, 2], xyzps[:, 3]
        # FT
        z_phase = self.defocus_base * zs.unsqueeze(1).unsqueeze(1)
        uin = self.aper * torch.exp(1j * (z_phase + self.phase_mask))  # rank 3
        uout = fft.fftshift(fft.fftn(fft.ifftshift(uin, dim=(1, 2)), dim=(1, 2)), dim=(1, 2))
        psfs = torch.abs(uout) ** 2

        # smoothing the PSF before down-sampling
        psfs = psfs.unsqueeze(1).type(torch.float32)  # rank 4
        psfs = F.conv2d(psfs, self.blur_kernel0.unsqueeze(0).unsqueeze(0), padding='same')  # apply blur before scaling, important!!

        # scaling and rotation
        at_m1 = (self.scaling() @ self.rotation()).unsqueeze(0).repeat(Nemitters, 1, 1)
        grid = F.affine_grid(at_m1[:, :2, :], psfs.size(), align_corners=False)
        psfs = F.grid_sample(psfs, grid, align_corners=False)
        # translation
        at_m2 = self.translation(xys)
        grid = F.affine_grid(at_m2[:, :2, :], psfs.size(), align_corners=False)
        psfs = F.grid_sample(psfs, grid, align_corners=False)

        psfs = psfs.squeeze(1)  # rank 3
        # photon normalization before cropping
        psfs = psfs / torch.sum(psfs, dim=(1, 2), keepdim=True) * photons.unsqueeze(1).unsqueeze(1)
        # crop
        psfs = psfs[:, self.r0:self.r0 + self.im_size[0], self.c0:self.c0 + self.im_size[1]]

        return psfs

    def forward(self, xyzps):  # two times of blur operation
        Nemitters = xyzps.shape[0]
        xys, zs, photons = xyzps[:, :2], xyzps[:, 2], xyzps[:, 3]

        z_phase = self.defocus_base * zs.unsqueeze(1).unsqueeze(1)
        uin = self.aper * torch.exp(1j * (z_phase + self.phase_mask))  # rank 3
        uout = fft.fftshift(fft.fftn(fft.ifftshift(uin, dim=(1, 2)), dim=(1, 2)), dim=(1, 2))
        psfs = torch.abs(uout) ** 2

        psfs = psfs.unsqueeze(1).type(torch.float32)  # rank 4
        psfs = F.conv2d(psfs, self.blur_kernel0.unsqueeze(0).unsqueeze(0), padding='same')  # apply blur before scaling

        # scaling, rotation, and translation
        at_m1 = (self.scaling() @ self.rotation()).unsqueeze(0).repeat(Nemitters, 1, 1)
        grid = F.affine_grid(at_m1[:, :2, :], psfs.size(), align_corners=False)
        psfs = F.grid_sample(psfs, grid, align_corners=False)
        at_m2 = self.translation(xys)
        grid = F.affine_grid(at_m2[:, :2, :], psfs.size(), align_corners=False)
        psfs = F.grid_sample(psfs, grid, align_corners=False)

        # same blur kernel for each PSF in the stack
        blur_kernel = torch.exp(-0.5 * (self.g_xx ** 2 + self.g_yy ** 2) / self.sigma1 ** 2)
        psfs = F.conv2d(psfs, blur_kernel.unsqueeze(0).unsqueeze(0), padding='same')
        psfs = psfs.squeeze(1)  # rank 3

        # photon normalization
        psfs = psfs / torch.sum(psfs, dim=(1, 2), keepdim=True) * photons.unsqueeze(1).unsqueeze(1)
        # crop
        psfs = psfs[:, self.r0:self.r0 + self.im_size[0], self.c0:self.c0 + self.im_size[1]]

        return psfs


class LensPSF(LensModelBase):
    def __init__(self, param_dict):
        super().__init__(param_dict)

    def forward(self, phase_mask, xyzs, photons, bgs):
        Nemitters = xyzs.shape[0]
        xys, zs = xyzs[:, :2], xyzs[:, 2]

        z_phase = self.defocus_base * zs.unsqueeze(1).unsqueeze(1)
        uin = self.aper * torch.exp(1j * (z_phase + phase_mask))  # rank 3
        uout = fft.fftshift(fft.fftn(fft.ifftshift(uin, dim=(1, 2)), dim=(1, 2)), dim=(1, 2))
        psfs = torch.abs(uout) ** 2

        psfs = psfs.unsqueeze(1).type(torch.float32)  # rank 4
        psfs = F.conv2d(psfs, self.blur_kernel0.unsqueeze(0).unsqueeze(0), padding='same')  # apply blur before scaling

        # scaling, rotation, and translation
        at_m1 = (self.scaling() @ self.rotation()).unsqueeze(0).repeat(Nemitters, 1, 1)
        grid = F.affine_grid(at_m1[:, :2, :], psfs.size(), align_corners=False)
        psfs = F.grid_sample(psfs, grid, align_corners=False)
        at_m2 = self.translation(xys)
        grid = F.affine_grid(at_m2[:, :2, :], psfs.size(), align_corners=False)
        psfs = F.grid_sample(psfs, grid, align_corners=False)

        # same blur kernel for each PSF in the stack
        blur_kernel = torch.exp(-0.5 * (self.g_xx ** 2 + self.g_yy ** 2) / self.sigma1 ** 2)
        psfs = F.conv2d(psfs, blur_kernel.unsqueeze(0).unsqueeze(0), padding='same')
        psfs = psfs.squeeze(1)  # rank 3

        # photon normalization
        psfs = psfs / torch.sum(psfs, dim=(1, 2), keepdim=True) * torch.abs(photons.unsqueeze(1)) + torch.abs(
            bgs.unsqueeze(1))
        # crop
        psfs = psfs[:, self.r0:self.r0 + self.im_size[0], self.c0:self.c0 + self.im_size[1]]

        return psfs


class LensModelTraining(LensModelBase):
    def __init__(self, param_dict):
        super().__init__(param_dict)  # initiate the base model
        device = param_dict['device']
        self.sigma_range = param_dict['sigma_range']
        self.baseline = param_dict['baseline']
        self.read_std = param_dict['read_std']
        # self.phase_mask = torch.tensor(param_dict['phase_mask'], device=device)
        self.bit_depth = param_dict['bit_depth']

    def blur_kernels(self, Nemitters):
        std_min, std_max = self.sigma_range
        stds = (std_min + (std_max - std_min) * torch.rand((Nemitters, 1))).to(self.device)
        gaussian_kernels = [torch.exp(-0.5 * (self.g_xx ** 2 + self.g_yy ** 2) / stds[i] ** 2) for i in
                            range(Nemitters)]
        gaussian_kernels = [kernel / kernel.sum() for kernel in gaussian_kernels]  # normalization
        gaussian_kernels = torch.stack(gaussian_kernels)
        return gaussian_kernels

    def forward(self, xyzps):
        psfs = self.get_psfs(xyzps)

        psfs = psfs.unsqueeze(0).type(torch.float32)
        blur_kernels = self.blur_kernels(xyzps.shape[0]).unsqueeze(0)
        im = F.conv2d(psfs, blur_kernels, padding='same').squeeze()  # psfs are summed to for an image

        # im = torch.poisson(im)  # shot noise  ##!!!
        im = (im + torch.poisson(im)) / 2  # * torch.rand(1, device=self.device)*1.0 # shot noise

        # im[im < 5] = 0
        # read_baseline = self.baseline[0] + torch.rand(1, device=self.device) * (self.baseline[1] - self.baseline[0])
        # read_std = self.read_std[0] + torch.rand(1, device=self.device) * (self.read_std[1] - self.read_std[0])
        # im += torch.round(read_baseline + torch.randn(im.shape, device=self.device) * read_std)  # read noise
        max_adu = 2 ** self.bit_depth - 1
        im[im > max_adu] = max_adu
        if self.bit_depth == 8:
            im = im.cpu().numpy().astype(np.uint8)
        elif self.bit_depth == 16:
            im = im.cpu().numpy().astype(np.uint16)
        return im


