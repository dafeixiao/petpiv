
from torch.utils.data import Dataset
import numpy as np
from skimage import io
import os
from torch import nn
import torch
from monai.networks.blocks import Convolution
import torch.nn.functional as F
import copy

class VolNet(nn.Module):
    def __init__(self, nc=64, nz=100):
        super().__init__()
        gn = ("GROUP", {"num_groups": 16})
        self.conv1 = Convolution(spatial_dims=2, in_channels=1, out_channels=nc, norm=gn, dilation=1, padding=1)
        self.conv2 = Convolution(spatial_dims=2, in_channels=nc+1, out_channels=nc, norm=gn, dilation=1, padding=1)
        self.conv3 = Convolution(spatial_dims=2, in_channels=nc+1, out_channels=nc, norm=gn, dilation=2, padding=2)
        self.conv4 = Convolution(spatial_dims=2, in_channels=nc+1, out_channels=nc, norm=gn, dilation=2, padding=2)
        self.conv5 = Convolution(spatial_dims=2, in_channels=nc+1, out_channels=nc, norm=gn, dilation=4, padding=4)
        self.conv6 = Convolution(spatial_dims=2, in_channels=nc+1, out_channels=nc, norm=gn, dilation=4, padding=4)
        self.conv7 = Convolution(spatial_dims=2, in_channels=nc,  out_channels=nz, kernel_size=1, padding=0, conv_only=True)

    def forward(self, x):
        y = self.conv1(x)
        y = self.conv2(torch.cat((y, x), 1)) + y
        y = self.conv3(torch.cat((y, x), 1)) + y
        y = self.conv4(torch.cat((y, x), 1)) + y
        y = self.conv5(torch.cat((y, x), 1)) + y
        y = self.conv6(torch.cat((y, x), 1)) + y
        y = F.softplus(self.conv7(y))   # positive intensity output

        return y


class UVWHead(nn.Module):
    def __init__(self, dt1_dict, dt2_dict, head_ds):
        super().__init__()
        gn = ("GROUP", {"num_groups": 8})

        self.dt1 = copy.deepcopy(dt1_dict)
        self.dt2 = copy.deepcopy(dt2_dict)

        self.use_dt1 = self.dt1["use"]
        self.use_dt2 = self.dt2["use"]

        if not (self.use_dt1 or self.use_dt2):
            raise ValueError("At least one of dt1 or dt2 must be enabled.")

        pool = nn.AvgPool3d(head_ds, head_ds) if any(s > 1 for s in head_ds) else nn.Identity()

        if self.use_dt1:
            self.compress1 = Convolution(
                3,
                in_channels=self.dt1["K"],
                out_channels=self.dt1["cv_nc"],
                kernel_size=1,
                padding=0,
                norm=gn,
            )

            self.head1 = nn.Sequential(
                Convolution(3, in_channels=self.dt1["f_nc"],
                            out_channels=self.dt1["cv_nc"], padding=1, norm=gn),
                pool,
                Convolution(3, in_channels=self.dt1["cv_nc"], out_channels=32, padding=1, norm=gn),
                nn.Dropout3d(p=0.1),
                Convolution(3, in_channels=32, out_channels=3, padding=1, conv_only=True),
            )
        else:
            self.compress1 = None
            self.head1 = None

        if self.use_dt2:
            self.compress2 = Convolution(
                3,
                in_channels=self.dt2["K"],
                out_channels=self.dt2["cv_nc"],
                kernel_size=1,
                padding=0,
                norm=gn,
            )

            self.head2 = nn.Sequential(
                Convolution(3, in_channels=self.dt2["f_nc"],
                            out_channels=self.dt2["cv_nc"], padding=1, norm=gn),
                pool,
                Convolution(3, in_channels=self.dt2["cv_nc"], out_channels=32, padding=1, norm=gn),
                nn.Dropout3d(p=0.1),
                Convolution(3, in_channels=32, out_channels=3, padding=1, conv_only=True),
            )
        else:
            self.compress2 = None
            self.head2 = None

        if self.use_dt1 and self.use_dt2:
            self.gate = nn.Sequential(
                Convolution(3, in_channels=self.dt1["f_nc"] + self.dt2["f_nc"],
                            out_channels=self.dt1["cv_nc"]+self.dt2["cv_nc"],
                            kernel_size=1,
                            padding=0,
                            norm=gn,
                            ),
                pool,
                Convolution(3, in_channels=self.dt1["cv_nc"]+self.dt2["cv_nc"], out_channels=32, padding=1, norm=gn),
                nn.Dropout3d(p=0.1),
                Convolution(3, in_channels=32, out_channels=3, kernel_size=1, padding=0, conv_only=True),
                nn.Sigmoid(),
            )
        else:
            self.gate = None


    def forward(self, f1=None, f2=None):
        if self.use_dt1 and f1 is None:
            raise ValueError("f1 is required when dt1 is enabled.")
        if self.use_dt2 and f2 is None:
            raise ValueError("f2 is required when dt2 is enabled.")

        uvw1 = self.head1(f1) if self.use_dt1 else None
        uvw2 = self.head2(f2) if self.use_dt2 else None

        if self.use_dt1 and self.use_dt2:

            alpha = self.gate(torch.cat([f1, f2], dim=1))
            uvw = alpha * uvw1 + (1.0 - alpha) * uvw2

            # uvw = (uvw1 + uvw2)/2

        elif self.use_dt1:
            uvw = uvw1
        else:
            uvw = uvw2

        return uvw


class Net0317(nn.Module):
    def __init__(self, nt, vol_dict, cv_dict, dt1_dict, dt2_dict):
        super().__init__()

        dt1_dict = copy.deepcopy(dt1_dict)
        dt2_dict = copy.deepcopy(dt2_dict)

        vol_nc = vol_dict['vol_nc']
        vol_nz = vol_dict['vol_nz']
        vol_ds = vol_dict['vol_ds']
        ds_zyx = vol_dict['ds_zyx']

        win = cv_dict['win']
        cv_ds = cv_dict['cv_ds']

        use_dt1 = dt1_dict['use']
        step1 = dt1_dict['step']
        md1_zyx = dt1_dict['md_zyx']
        cv1_nc = dt1_dict['cv_nc']

        use_dt2 = dt2_dict['use']
        step2 = dt2_dict['step']
        md2_zyx = dt2_dict['md_zyx']
        cv2_nc = dt2_dict['cv_nc']

        if use_dt1:
            assert nt > step1, "nt does not support dt2 step"
        if use_dt2:
            assert nt > step2, "nt does not support dt2 step"

        self.nt = nt
        self.vol_ds = vol_ds
        self.win = win
        self.cv_ds = cv_ds
        self.use_dt1 = use_dt1
        self.use_dt2 = use_dt2
        self.md1_zyx = md1_zyx
        self.md2_zyx = md2_zyx
        self.step1 = step1
        self.step2 = step2

        head_ds = (
            ds_zyx[0] // vol_ds[0] // cv_ds[0],
            ds_zyx[1] // vol_ds[1] // cv_ds[1],
            ds_zyx[2] // vol_ds[2] // cv_ds[2],
        )

        print(f'image volume downsampling: {vol_ds}, cv downsampling: {cv_ds}, head downsampling: {head_ds}')

        if use_dt1:
            mdz, mdy, mdx = md1_zyx
            K1 = (2 * mdz + 1) * (2 * mdy + 1) * (2 * mdx + 1)
            print(f'dt1 cost volume channels: {K1}')
            f1_nc = int(cv1_nc * (nt - step1))
        else:
            K1 = 0
            f1_nc = 0

        if use_dt2:
            mdz, mdy, mdx = md2_zyx
            K2 = (2 * mdz + 1) * (2 * mdy + 1) * (2 * mdx + 1)
            print(f'dt2 cost volume channels: {K2}')
            f2_nc = int(cv2_nc * (nt - step2))
        else:
            K2 = 0
            f2_nc = 0

        dt1_dict['K'] = K1
        dt1_dict['f_nc'] = f1_nc
        dt2_dict['K'] = K2
        dt2_dict['f_nc'] = f2_nc

        self.layer1 = VolNet(nc=vol_nc, nz=vol_nz)
        self.layer2 = UVWHead(dt1_dict, dt2_dict, head_ds)

    def cost_volume_3d(self, a, b, md_zyx, eps=1e-6):
        mdz, mdy, mdx = md_zyx
        B, Z, H, W = a.shape
        Dz, Dy, Dx = 2 * mdz + 1, 2 * mdy + 1, 2 * mdx + 1
        K = Dz * Dy * Dx
        b_pad = F.pad(b, (mdx, mdx, mdy, mdy, mdz, mdz))
        bz = b_pad.unfold(1, Z, 1)
        bzh = bz.unfold(2, H, 1)
        bzhw = bzh.unfold(3, W, 1)
        b_shifts = bzhw.contiguous().view(B, K, Z, H, W)
        cv = a.unsqueeze(1) * b_shifts
        a2 = a * a
        b2_shifts = b_shifts * b_shifts
        pad = tuple(w // 2 for w in self.win)
        cv_w = F.avg_pool3d(cv, kernel_size=self.win, stride=self.cv_ds, padding=pad)
        a2_w = F.avg_pool3d(a2.unsqueeze(1), kernel_size=self.win, stride=self.cv_ds, padding=pad)
        b2_w = F.avg_pool3d(b2_shifts, kernel_size=self.win, stride=self.cv_ds, padding=pad)
        denom = torch.sqrt(a2_w * b2_w + eps).clamp_min(eps)
        cv_norm = cv_w / denom
        return cv_norm

    def forward(self, x):
        B, T, H, W = x.shape
        assert T == self.nt, "input - T mismatch"

        x_ = x.reshape(B * T, 1, H, W)
        vx = self.layer1(x_)
        vol = F.avg_pool3d(vx.unsqueeze(1), kernel_size=self.vol_ds, stride=self.vol_ds).squeeze(1)
        vol = vol.reshape(B, T, vol.shape[1], vol.shape[2], vol.shape[3])

        f1, f2 = None, None
        if self.use_dt1:
            pairs = [(i, i + self.step1) for i in range(T - self.step1)]
            cvs = [self.layer2.compress1(self.cost_volume_3d(vol[:, i], vol[:, j], self.md1_zyx)) for i, j in pairs]
            f1 = torch.cat(cvs, dim=1)

        if self.use_dt2:
            pairs = [(i, i + self.step2) for i in range(T - self.step2)]
            cvs = [self.layer2.compress2(self.cost_volume_3d(vol[:, i], vol[:, j], self.md2_zyx)) for i, j in pairs]
            f2 = torch.cat(cvs, dim=1)

        y = self.layer2(f1, f2)

        return y


class MyDS(Dataset):
    def __init__(self, list_IDs, labels, xfolder, yfolder, nt):
        self.list_IDs = list_IDs
        self.labels = labels
        self.xfolder = xfolder
        self.yfolder = yfolder
        # take nt frames around the center which is 5 (among 11 frames)
        self.t0 = 5 - nt//2
        self.t1 = 5 + nt//2 +1


    # total number of samples in the dataset
    def __len__(self):
        return len(self.list_IDs)

    # sampling one example from the data
    def __getitem__(self, index):
        ID = self.list_IDs[index]
        y = np.load(os.path.join(self.yfolder, ID)).astype(np.float32)
        x_files = self.labels[ID]
        x_files = x_files[self.t0:self.t1]
        ims = np.stack([io.imread(os.path.join(self.xfolder, im_name)) for im_name in x_files])
        x = ims.astype(np.float32)/255.0  # image normalization
        return x, y


class MyLoss(nn.Module):
    def __init__(self, lam_w=1.0):
        super().__init__()
        self.lam_w = lam_w

    def forward(self, y, y_gt):
        Lu = F.l1_loss(y[:, 0], y_gt[:, 0])
        Lv = F.l1_loss(y[:, 1], y_gt[:, 1])
        Lw = F.l1_loss(y[:, 2], y_gt[:, 2])
        loss_data = (Lu + Lv + self.lam_w * Lw) / (2.0 + self.lam_w)  # scale-stable

        # loss_data = 100*F.mse_loss(y, y_gt)
        loss = loss_data
        return loss