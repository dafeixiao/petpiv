

import pandas as pd
import re
from scipy.interpolate import RBFInterpolator
from sklearn.metrics.pairwise import pairwise_distances
from scipy.optimize import linear_sum_assignment
import numpy as np

def read_tomo_file(file_path):
    columns = None
    data_start = None

    with open(file_path, "r", errors="ignore") as f:
        for i, line in enumerate(f):
            if line.startswith("VARIABLES"):
                columns = re.findall(r'"([^"]+)"', line)
            elif line.strip().startswith("DATAPACKING"):
                data_start = i + 1
                break
    if columns is None:
        raise ValueError("VARIABLES line not found")
    if data_start is None:
        raise ValueError("DATAPACKING line not found")
    df = pd.read_csv(
        file_path,
        sep=r"\s+",
        skiprows=data_start,
        names=columns,
        engine="python"
    )
    return df


def pair_func(xyz_set0, xyz_set1, radius):
    # xyz_set0 is ground truth

    # calculate the distance matrix for each GT to each prediction
    C = pairwise_distances(xyz_set0, xyz_set1, 'euclidean')
    # number of recovered points and GT sources
    num0 = xyz_set0.shape[0]
    num1 = xyz_set1.shape[0]

    # find the matching using the Hungarian algorithm
    ids0, ids1 = linear_sum_assignment(C)

    # number of matched points
    num_matches = len(ids0)

    # run over matched points and filter points radius away from GT
    indicatorTP = [C[ids0[i], ids1[i]] < radius for i in range(num_matches)]

    # resulting TP count
    TP = sum(indicatorTP)

    # if there's TP
    if TP:
        # pairs of TP
        ids0_TP = (ids0[indicatorTP]).tolist()
        ids1_TP = (ids1[indicatorTP]).tolist()

        return ids0_TP, ids1_TP
    else:
        return None


class XYZTransform:
    def __init__(self, param_dict):
        # param_dict--affine3x4 that maps from xyz_set1 to xyz_set0
        affine_matrix = param_dict['affine3x4']
        A = affine_matrix[:, :3]  # (3,3)
        t = affine_matrix[:, 3]  # (3,)

        xyz_set0 = param_dict['xyz_set0']  # Nx3
        xyz_set1 = param_dict['xyz_set1']  # Nx3
        residual = xyz_set0 - (xyz_set1 @ A.T + t)

        ffd = RBFInterpolator(
            xyz_set1,
            residual,
            kernel='cubic',  # 'thin_plate_spline'
            smoothing=100.0  # START HERE, then tune
        )
        self.A = A
        self.t = t
        self.ffd = ffd
        self.xyz_sign = param_dict['xyz_sign']
        self.xyz_shift = param_dict['xyz_shift']

    def direct_transform(self, xyzs):
        xyzs_aff = xyzs @ self.A.T + self.t + self.ffd(xyzs)
        return xyzs_aff

    def __call__(self, xyzs):
        # xyzs--Nx3
        xyzs[:, 0] = xyzs[:, 0] * self.xyz_sign[0] + self.xyz_shift[0]
        xyzs[:, 1] = xyzs[:, 1] * self.xyz_sign[1] + self.xyz_shift[1]
        xyzs[:, 2] = xyzs[:, 2] * self.xyz_sign[2] + self.xyz_shift[2]

        xyzs_aff = xyzs @ self.A.T + self.t + self.ffd(xyzs)
        return xyzs_aff


def calc_dice_mae(xyz_gt, xyz_rec, radius):
    # if the net didn't detect anything return None's
    if xyz_rec is None:
        print("Empty Prediction!")
        return 0.0, None, None, None
    else:
        # calculate the distance matrix for each GT to each prediction
        C = pairwise_distances(xyz_rec, xyz_gt, 'euclidean')
        # number of recovered points and GT sources
        num_rec = xyz_rec.shape[0]
        num_gt = xyz_gt.shape[0]
        # find the matching using the Hungarian algorithm
        rec_ind, gt_ind = linear_sum_assignment(C)
        # number of matched points
        num_matches = len(rec_ind)
        # run over matched points and filter points radius away from GT
        indicatorTP = [False] * num_matches
        for i in range(num_matches):
            # if the point is closer than radius then TP else it's FP
            if C[rec_ind[i], gt_ind[i]] < radius:
                indicatorTP[i] = True
        # resulting TP count
        TP = sum(indicatorTP)
        # resulting jaccard index
        # tp_ratio = TP / (num_rec + num_gt - TP)
        dice = 2*TP / (num_rec + num_gt)
        # tp_ratio = TP / num_rec
        # tp_ratio = TP / num_gt
        # if there's TP
        if TP:
            # pairs of TP
            rec_ind_TP = (rec_ind[indicatorTP]).tolist()
            gt_ind_TP = (gt_ind[indicatorTP]).tolist()
            xyz_rec_TP = xyz_rec[rec_ind_TP, :]
            xyz_gt_TP = xyz_gt[gt_ind_TP, :]
            # calculate mean RMSE in xy, z, and xyz
            mae_xy = np.mean(np.sqrt(np.sum((xyz_rec_TP[:, :2] - xyz_gt_TP[:, :2]) ** 2, 1)))
            mae_z = np.mean(np.sqrt(np.sum((xyz_rec_TP[:, 2:] - xyz_gt_TP[:, 2:]) ** 2, 1)))
            mae_xyz = np.mean(np.sqrt(np.sum((xyz_rec_TP - xyz_gt_TP) ** 2, 1)))
            return dice, mae_xy, mae_z, mae_xyz
        else:
            return dice, None, None, None


def image_cropping(img, num_r, num_c):
    hh, ww = img.shape
    h = int(np.ceil(hh / num_r))
    w = int(np.ceil(ww / num_c))

    sub_imgs = []
    for i in range(num_r):
        c_list = []
        for j in range(num_c):
            sub_img = img[i * h:min((i + 1) * h, hh), j * w:min((j + 1) * w, ww)]
            c_list.append(sub_img)
        sub_imgs.append(c_list)

    return sub_imgs


def combine_tr_sm(tr, sm):
    # tr, ndarray, Nx6 [track id, x, y, z, intensity, frame index]
    # sm, ndarray, Mx11, [track id, x, y, z, u, v, w, ax, ay, az, frame index], M<N- the head and taile of each track is filtered out

    tr = np.hstack([tr[:, :4], np.zeros((tr.shape[0], 6), dtype=tr.dtype), tr[:, 5:]])  # add uvw and axyz columns
    lookup = {(pid, t): row[4:10] for row, pid, t in zip(sm, sm[:, 0], sm[:, -1])}
    for i in range(tr.shape[0]):
        key = (tr[i, 0], tr[i, -1])
        if key in lookup:
            tr[i, 4:10] = lookup[key]

    return tr


from pathlib import Path


def export_psfe_to_dat(ps_list, outdir, base_name="B", eps=1e-6):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    variables = [
        "x[mm]", "y[mm]", "z[mm]", "I",
        "Vx[m/s]", "Vy[m/s]", "Vz[m/s]", "|V|[m/s]",
        "trackID",
        "Ax[m/s²]", "Ay[m/s²]", "Az[m/s²]", "|A|[m/s²]",
        "Absolute intensity[]"
    ]

    for idx, ps in enumerate(ps_list):
        trackID = ps[:, 0]
        x, y, z = ps[:, 1], ps[:, 2], ps[:, 3]
        u, v, w = ps[:, 4], ps[:, 5], ps[:, 6]
        ax, ay, az = ps[:, 7], ps[:, 8], ps[:, 9]

        Vmag = np.sqrt(u ** 2 + v ** 2 + w ** 2)
        Amag = np.sqrt(ax ** 2 + ay ** 2 + az ** 2)

        filler = np.full_like(trackID, eps)

        out = np.column_stack([
            x, y, z,
            filler,  # I
            u, v, w,
            Vmag,
            trackID,
            ax, ay, az,
            Amag,
            filler  # Absolute intensity
        ])

        idx = idx + 1

        fname = outdir / f"{base_name}{(idx):06d}.dat"

        with open(fname, "w") as f:
            f.write(f'TITLE = "{base_name}{idx:06d}"\n')
            f.write(
                'VARIABLES = ' +
                ' '.join(f'"{v}"' for v in variables) +
                '\n'
            )
            f.write(f'ZONE T="Snapshot {idx:04d}"\n')
            f.write(f'STRANDID=1, SOLUTIONTIME=None\n')
            f.write(
                f'I={out.shape[0]}, J=1, K=1, ZONETYPE=Ordered\n'
            )
            f.write('DATAPACKING=POINT\n')

            np.savetxt(f, out, fmt="%.6g")

        # print(f"written: {fname}")
    print('done')



def read_IJK_from_tecplot_dat(fname):
    """
    Read I, J, K dimensions from the first ZONE line
    of a Tecplot ASCII .dat file.

    Parameters
    ----------
    fname : str
        Path to the .dat file

    Returns
    -------
    I, J, K : int
        Grid dimensions
    """
    I = J = K = None

    with open(fname, "r") as f:
        for line in f:
            if "ZONE" in line:
                mI = re.search(r"I\s*=\s*(\d+)", line)
                mJ = re.search(r"J\s*=\s*(\d+)", line)
                mK = re.search(r"K\s*=\s*(\d+)", line)

                if mI:
                    I = int(mI.group(1))
                if mJ:
                    J = int(mJ.group(1))
                if mK:
                    K = int(mK.group(1))

                break

    if I is None or J is None or K is None:
        raise ValueError("Could not find I, J, K in ZONE header")

    return I, J, K