
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import dctn
import plotly.graph_objects as go
import matplotlib as mpl

# show three cross-sections of volumetric data
def plot_quadrants(ax, data, fixed_coord, maxv, minv, cmap='viridis_r'):
    array = data
    D, H, W = array.shape
    xx, yy, zz = np.meshgrid(np.arange(0, W), np.arange(0, H), np.arange(0, D), indexing='xy')
    xxx, yyy, zzz = xx.transpose(2, 0, 1), yy.transpose(2, 0, 1), zz.transpose(2, 0, 1)

    nz, ny, nx = array.shape
    index = {
        'x': (slice(None), slice(None), nx // 2),
        'y': (slice(None), ny//2, slice(None)),
        'z': (nz//2, slice(None), slice(None)),
    }[fixed_coord]

    plane_data = array[index]
    plane_x = xxx[index]
    plane_y = yyy[index]
    plane_z = zzz[index]

    n0, n1 = plane_data.shape
    row_split = n0 // 2
    column_split = n1 // 2
    quadrants = [
        plane_data[:row_split+1, :column_split+1],
        plane_data[:row_split+1, column_split:],
        plane_data[row_split:, :column_split+1],
        plane_data[row_split:, column_split:]
    ]

    quadrants_x = [
        plane_x[:row_split+1, :column_split+1],
        plane_x[:row_split+1, column_split:],
        plane_x[row_split:, :column_split+1],
        plane_x[row_split:, column_split:]
    ]

    quadrants_y = [
        plane_y[:row_split+1, :column_split+1],
        plane_y[:row_split+1, column_split:],
        plane_y[row_split:, :column_split+1],
        plane_y[row_split:, column_split:]
    ]
    quadrants_z = [
        plane_z[:row_split+1, :column_split+1],
        plane_z[:row_split+1, column_split:],
        plane_z[row_split:, :column_split+1],
        plane_z[row_split:, column_split:]
    ]
    min_val = minv
    max_val = maxv

    cmap = plt.get_cmap(cmap)

    for i, quadrant in enumerate(quadrants):
        facecolors = cmap((quadrant - min_val) / (max_val - min_val))

        if fixed_coord == 'x':
            X, Y, Z = quadrants_x[i], quadrants_y[i], quadrants_z[i]
            ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=facecolors, shade=False)

        elif fixed_coord == 'y':
            X, Y, Z = quadrants_x[i], quadrants_y[i], quadrants_z[i]
            ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=facecolors, shade=False)
        elif fixed_coord == 'z':
            X, Y, Z = quadrants_x[i], quadrants_y[i], quadrants_z[i]
            ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=facecolors, shade=False)






def plot_cross_sections(dhw_data, cmap='viridis_r'):
    D, H, W = dhw_data.shape
    xx, yy, zz = np.meshgrid(np.arange(0, W), np.arange(0, H), np.arange(0, D), indexing='xy')
    xx, yy, zz = xx.transpose(2, 0, 1), yy.transpose(2, 0, 1), zz.transpose(2, 0, 1)

    # xxx, yyy, zzz, array = data
    nz, ny, nx = xx.shape
    plane_dict = {
        'x': (slice(None), slice(None), nx // 2),
        'y': (slice(None), ny//2, slice(None)),
        'z': (nz//2, slice(None), slice(None)),
    }

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for fixed_coord in list(plane_dict.keys()):
        plane_idx = plane_dict[fixed_coord]
        plane_data = dhw_data[plane_idx]
        plane_x = xx[plane_idx]
        plane_y = yy[plane_idx]
        plane_z = zz[plane_idx]

        n0, n1 = plane_data.shape
        row_split = n0 // 2
        column_split = n1 // 2

        quadrants = [
            plane_data[:row_split+1, :column_split+1],
            plane_data[:row_split+1, column_split:],
            plane_data[row_split:, :column_split+1],
            plane_data[row_split:, column_split:]
        ]

        quadrants_x = [
            plane_x[:row_split+1, :column_split+1],
            plane_x[:row_split+1, column_split:],
            plane_x[row_split:, :column_split+1],
            plane_x[row_split:, column_split:]
        ]

        quadrants_y = [
            plane_y[:row_split+1, :column_split+1],
            plane_y[:row_split+1, column_split:],
            plane_y[row_split:, :column_split+1],
            plane_y[row_split:, column_split:]
        ]
        quadrants_z = [
            plane_z[:row_split+1, :column_split+1],
            plane_z[:row_split+1, column_split:],
            plane_z[row_split:, :column_split+1],
            plane_z[row_split:, column_split:]
        ]

        min_val, max_val = dhw_data.min(), dhw_data.max()

        cmap = plt.get_cmap(cmap)
        for i, quadrant in enumerate(quadrants):
            facecolors = cmap((quadrant - min_val) / (max_val - min_val))

            if fixed_coord == 'x':
                X, Y, Z = quadrants_x[i], quadrants_y[i], quadrants_z[i]
                ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=facecolors, shade=False)

            elif fixed_coord == 'y':
                X, Y, Z = quadrants_x[i], quadrants_y[i], quadrants_z[i]
                ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=facecolors, shade=False)
            elif fixed_coord == 'z':
                X, Y, Z = quadrants_x[i], quadrants_y[i], quadrants_z[i]
                ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=facecolors, shade=False)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()


def plot_uvw(uvw1, uvw2):
    fig = plt.figure()


    # uu = np.stack((uvw1[0], uvw2[0]))
    # maxv, minv = uu.max(), uu.min()
    maxv, minv = uvw2[0].max(), uvw2[0].min()
    ax = fig.add_subplot(321, projection='3d')
    plot_quadrants(ax, uvw1[0], 'x', maxv, minv, cmap='viridis_r')
    plot_quadrants(ax, uvw1[0], 'y', maxv, minv, cmap='viridis_r')
    plot_quadrants(ax, uvw1[0], 'z', maxv, minv, cmap='viridis_r')
    ax = fig.add_subplot(322, projection='3d')
    plot_quadrants(ax, uvw2[0], 'x', maxv, minv, cmap='viridis_r')
    plot_quadrants(ax, uvw2[0], 'y', maxv, minv, cmap='viridis_r')
    plot_quadrants(ax, uvw2[0], 'z', maxv, minv, cmap='viridis_r')

    # vv = np.stack((uvw1[1], uvw2[1]))
    maxv, minv = uvw2[1].max(), uvw2[1].min()
    ax = fig.add_subplot(323, projection='3d')
    plot_quadrants(ax, uvw1[1], 'x', maxv, minv, cmap='viridis_r')
    plot_quadrants(ax, uvw1[1], 'y', maxv, minv, cmap='viridis_r')
    plot_quadrants(ax, uvw1[1], 'z', maxv, minv, cmap='viridis_r')
    ax = fig.add_subplot(324, projection='3d')
    plot_quadrants(ax, uvw2[1], 'x', maxv, minv, cmap='viridis_r')
    plot_quadrants(ax, uvw2[1], 'y', maxv, minv, cmap='viridis_r')
    plot_quadrants(ax, uvw2[1], 'z', maxv, minv, cmap='viridis_r')

    # ww = np.stack((uvw1[2], uvw2[2]))
    maxv, minv = uvw2[2].max(), uvw2[2].min()
    ax = fig.add_subplot(325, projection='3d')
    plot_quadrants(ax, uvw1[2], 'x', maxv, minv, cmap='viridis_r')
    plot_quadrants(ax, uvw1[2], 'y', maxv, minv, cmap='viridis_r')
    plot_quadrants(ax, uvw1[2], 'z', maxv, minv, cmap='viridis_r')
    ax = fig.add_subplot(326, projection='3d')
    plot_quadrants(ax, uvw2[2], 'x', maxv, minv, cmap='viridis_r')
    plot_quadrants(ax, uvw2[2], 'y', maxv, minv, cmap='viridis_r')
    plot_quadrants(ax, uvw2[2], 'z', maxv, minv, cmap='viridis_r')

    plt.show()


def plot_vol2(vol1, vol0):  # vol0 is the reference
    fig = plt.figure(figsize=(10, 5))

    maxv, minv = vol0.max(), vol0.min()
    ax1 = fig.add_subplot(121, projection='3d')
    plot_quadrants(ax1, vol1, 'x', maxv, minv, cmap='turbo')
    plot_quadrants(ax1, vol1, 'y', maxv, minv, cmap='turbo')
    plot_quadrants(ax1, vol1, 'z', maxv, minv, cmap='turbo')
    ax2 = fig.add_subplot(122, projection='3d')
    plot_quadrants(ax2, vol0, 'x', maxv, minv, cmap='turbo')
    plot_quadrants(ax2, vol0, 'y', maxv, minv, cmap='turbo')
    plot_quadrants(ax2, vol0, 'z', maxv, minv, cmap='turbo')

    cmap = plt.get_cmap('turbo')
    norm = mpl.colors.Normalize(vmin=minv, vmax=maxv)
    mappable = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array([])

    fig.colorbar(mappable, ax=[ax1, ax2], shrink=0.5, pad=0.05)


    plt.show()




def go_plot(data):
    D, H, W = data.shape
    xx, yy, zz = np.meshgrid(np.arange(0, W), np.arange(0, H), np.arange(0, D), indexing='xy')
    X, Y, Z = xx.transpose(2, 0, 1), yy.transpose(2, 0, 1), zz.transpose(2, 0, 1)
    values = data
    # minv, maxv = data.min(), data.max()
    fig = go.Figure(data=go.Volume(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=values.flatten(),
        # isomin=minv,
        # isomax=maxv,
        opacity=0.1,  # needs to be small to see through all surfaces
        surface_count=21,  # needs to be a large number for good volume rendering
    ))
    fig.show()



# D, W, H = 16, 32, 32
# nkd, nkh, nkw = 2, 2, 2
# coeffs = (np.random.rand(nkd, nkh, nkw) * 2 - 1)*0.2
# u = dctn(np.pad(coeffs, ((0, D - nkd), (0, H - nkh), (0, W - nkw))), type=3)
# coeffs = (np.random.rand(nkd, nkh, nkw) * 2 - 1)*0.2
# v = dctn(np.pad(coeffs, ((0, D - nkd), (0, H - nkh), (0, W - nkw))), type=3)
# coeffs = (np.random.rand(nkd, nkh, nkw) * 2 - 1)*0.2
# w = dctn(np.pad(coeffs, ((0, D - nkd), (0, H - nkh), (0, W - nkw))), type=3)
# uvw = np.stack((u, v, w))

# # demo 1
# plot_cross_sections(u)
#
#
# # demo 2
#
# xx, yy, zz = np.meshgrid(np.arange(0, W), np.arange(0, H), np.arange(0, D), indexing='xy')
# xx, yy, zz = xx.transpose(2, 0, 1), yy.transpose(2, 0, 1), zz.transpose(2, 0, 1)
#
# array = u
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# minv, maxv = array.min(), array.max()
# # ax.set_box_aspect((32, 32, 16))
# plot_quadrants(ax, [xx, yy, zz, array], 'x', minv, maxv, cmap='viridis_r')
# plot_quadrants(ax, [xx, yy, zz, array], 'y', minv, maxv, cmap='viridis_r')
# plot_quadrants(ax, [xx, yy, zz, array], 'z', minv, maxv, cmap='viridis_r')
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')
# plt.show()

# plot_uvw(uvw, uvw)


#
# # demo 3
# go_plot(u)








