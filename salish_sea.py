import os
from matplotlib import patches
from matplotlib.colors import ListedColormap
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import matplotlib.patches as patches
from scipy.io import loadmat
from datetime import datetime, timedelta

class Salish_Sea:
    def __init__(self, ss_dir, bty_file_path=None):
        self.ss_dir = ss_dir
        self.files = np.sort(os.listdir(ss_dir))

        # Get the NetCDF Variables
        self.var_names = [
            'x', 'y', 'lon', 'lat',
            'xc', 'yc', 'lonc', 'latc',
            'siglay', 'siglev',
            'siglay_center', 'siglev_center',
            'h_center', 'h', 'time',
            'zeta', 'u', 'v', 'tauc',
            'uard_obcn', 'xflux_obc',
            'omega', 'ww',
            'ua', 'va',
            'temp', 'salinity',
            'viscofm',
            'km', 'kh', 'kq',
            'q2', 'q2l', 'l',
        ]

        data = Dataset(os.path.join(ss_dir, self.files[0]), 'r')
        for var in self.var_names:
            arr = data.variables[var][:]
            setattr(
                self,
                var,
                np.zeros((*arr.shape, len(self.files)), dtype=arr.dtype)
            )

        for i, fname in enumerate(self.files):
            ds = Dataset(os.path.join(ss_dir, fname), 'r')
            for var in self.var_names:
                getattr(self, var)[(..., i)] = ds.variables[var][:]

            ds.close()
        
        # Declare Datetimes
        self.dt = np.empty(self.time.shape, dtype='datetime64[ns]')
        for i in range(self.time.shape[0]):
            for j in range(self.time.shape[1]):
                mjd = self.time[i,j]
                dt = datetime(1858, 11, 17) + timedelta(days=float(mjd))
                self.dt[i, j] = dt

        # Unchanging Variables
        self.x = self.x[:,0]
        self.y = self.y[:,0]
        self.lon = self.lon[:,0]
        self.lat = self.lat[:,0]
        self.xc = self.xc[:,0]
        self.yc = self.yc[:,0]
        self.lonc = self.lonc[:,0]
        self.latc = self.latc[:,0]
        self.siglay = self.siglay[:, :, 0]
        self.siglev = self.siglev[:, :, 0]
        self.siglay_center = self.siglay_center[:, :, 0]
        self.siglev_center = self.siglev_center[:, :, 0]
        self.h = self.h[:, 0]
        self.h_center = self.h_center[:, 0]

        t_size = self.dt.shape[0] * self.dt.shape[1]
        time = np.empty(t_size, dtype='datetime64[ns]')
        zeta = np.zeros((self.zeta.shape[1], t_size))
        u = np.zeros((self.u.shape[1], self.u.shape[2], t_size))
        v = np.zeros((self.v.shape[1], self.v.shape[2], t_size))
        tauc = np.zeros((self.tauc.shape[1], t_size))
        uard_obcn = np.zeros((self.uard_obcn.shape[1], t_size))
        xflux_obc = np.zeros((self.xflux_obc.shape[1], self.xflux_obc.shape[2], t_size))
        omega = np.zeros((self.omega.shape[1], self.omega.shape[2], t_size))
        ww = np.zeros((self.ww.shape[1], self.ww.shape[2], t_size))
        ua = np.zeros((self.ua.shape[1], t_size))
        va = np.zeros((self.va.shape[1], t_size))
        temp = np.zeros((self.temp.shape[1], self.temp.shape[2], t_size))
        salinity = np.zeros((self.salinity.shape[1], self.salinity.shape[2], t_size))
        viscofm = np.zeros((self.viscofm.shape[1], self.viscofm.shape[2], t_size))
        km = np.zeros((self.km.shape[1], self.km.shape[2], t_size))
        kh = np.zeros((self.kh.shape[1], self.kh.shape[2], t_size))
        kq = np.zeros((self.kq.shape[1], self.kq.shape[2], t_size))
        q2 = np.zeros((self.q2.shape[1], self.q2.shape[2], t_size))
        q2l = np.zeros((self.q2l.shape[1], self.q2l.shape[2], t_size))
        l = np.zeros((self.l.shape[1], self.l.shape[2], t_size))

        # Sort for time columns
        t_idx = 0
        for d in range(self.dt.shape[1]):
            for h in range(self.dt.shape[0]):
                time[t_idx] = self.dt[h, d]
                zeta[:, t_idx] = self.zeta[h,:,d]
                u[:, :, t_idx] = self.u[h, :, :, d]
                v[:, :, t_idx] = self.v[h, :, :, d]
                tauc[:, t_idx] = self.tauc[h, :, d]
                uard_obcn[:, t_idx] = self.uard_obcn[h, :, d]
                xflux_obc[:, :, t_idx] = self.xflux_obc[h, :, :, d]
                omega[:, :, t_idx] = self.omega[h, :, :, d]
                ww[:, :, t_idx] = self.ww[h, :, :, d]
                ua[:, t_idx] = self.ua[h, :, d]
                va[:, t_idx] = self.va[h, :, d]
                temp[:, :, t_idx] = self.temp[h, :, :, d]
                salinity[:, :, t_idx] = self.salinity[h, :, :, d]
                viscofm[:, :, t_idx] = self.viscofm[h, :, :, d]
                km[:, :, t_idx] = self.km[h, :, :, d]
                kh[:, :, t_idx] = self.kh[h, :, :, d]
                kq[:, :, t_idx] = self.kq[h, :, :, d]
                q2[:, :, t_idx] = self.q2[h, :, :, d]
                q2l[:, :, t_idx] = self.q2l[h, :,:, d]
                l[:, :, t_idx] = self.l[h, :, :, d]

                # Increment the Time Index for the 
                t_idx += 1

        # Resave the sorted data
        self.time = time
        self.zeta = zeta
        self.u = u
        self.v = v
        self.tauc = tauc
        self.uard_obcn = uard_obcn
        self.xflux_obc = xflux_obc
        self.omega = omega
        self.ww = ww
        self.ua = ua
        self.va = va
        self.temp = temp
        self.salinity = salinity
        self.viscofm = viscofm
        self.km = km
        self.kh = kh
        self.kq = kq
        self.q2 = q2
        self.q2l = q2l
        self.l = l

        self.z = np.zeros((self.siglay.shape[0], self.siglay.shape[1], self.u.shape[2]))

        # Compute the vertical coordinate (z) based on the bathymetry and water level
        for d in range(self.z.shape[0]):
            for p in range(self.z.shape[1]):
                for t in range(self.z.shape[2]):
                    self.z[d,p,t] = -self.zeta[p,t] - self.siglay[d,p] * (self.h[p] + self.zeta[p,t])

        if bty_file_path is not None:
            self.bty = np.flipud(loadmat(bty_file_path)['bath_map']) * -1
            self.lon_range = loadmat(bty_file_path)['lon_range']
            self.lat_range = loadmat(bty_file_path)['lat_range']
            self.lon_min = np.min(self.lon_range)
            self.lon_max = np.max(self.lon_range)
            self.lat_min = np.min(self.lat_range)
            self.lat_max = np.max(self.lat_range)
            mask = (self.lon >= self.lon_min) & (self.lon <= self.lon_max) & (self.lat >= self.lat_min) & (self.lat <= self.lat_max)
            mask_c = (self.lonc >= self.lon_min) & (self.lonc <= self.lon_max) & (self.latc >= self.lat_min) & (self.latc <= self.lat_max)
            

    def plot_bty(self, trackline=None, save_path=None):
        # Plot Bathymetry
        plt.figure(figsize=(7, 10))
        plt.imshow(
            self.bty,
            extent=[
                self.lon_range.min(),
                self.lon_range.max(),
                self.lat_range.min(),
                self.lat_range.max()
            ],
            origin='lower',
            aspect='auto',
            cmap='plasma',
            vmin=-50,
            vmax=200
        )
        cbar =plt.colorbar()
        cbar.ax.tick_params(labelsize=12)
        cbar.set_label('Depth (m)', fontsize=12)

        # Add Trackline
        if trackline is not None:
            lon = trackline[:, 0]
            lat = trackline[:, 1]
            plt.plot(lon, lat, 'k--', linewidth=2)
            # Source (green dot)
            plt.plot(
                trackline[0,0],
                trackline[0,1],
                'go',   # g = green, o = circle
                markersize=8,
                label='Source'
            )

            # Receiver (red star)
            plt.plot(
                trackline[1,0],
                trackline[1,1],
                'r*',   # r = red, * = star
                markersize=12,
                label='Receiver'
            )

        plt.xlabel("Longitude (˚)", fontsize=12)
        plt.ylabel("Latitude (˚)", fontsize=12)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.xticks(rotation=45)
        
        legend = plt.legend()
        for text in legend.get_texts():
            text.set_fontsize(12)
            text.set_fontname('Arial')
            text.set_color('black')
        plt.tight_layout()

        if save_path is not None:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            plt.close()
        else:
            plt.show()


    def plot_current(self, start_t, end_t=None, save_path=None):

        start_t = np.datetime64(start_t)
        t0 = np.argmin(np.abs(self.time - start_t))

        # Static Plot for single time
        if end_t is None:

            fig, axes = plt.subplots(2, 5, figsize=(24, 10), constrained_layout=True)
            axes = axes.ravel()

            # Precompute land mask once
            land = np.ma.masked_where(self.bty >= 0, self.bty)

            for k, ax in enumerate(axes):

                u = self.u[k, :, t0]
                v = self.v[k, :, t0]
                z_layer = self.z[k, :, t0]

                sc = ax.scatter(
                    self.lon,
                    self.lat,
                    c=z_layer,
                    s=5,
                    cmap='viridis'
                )

                plt.colorbar(sc, ax=ax, label='Layer Depth (m)')

                ax.imshow(
                    land,
                    extent=[
                        self.lon_range.min(), self.lon_range.max(),
                        self.lat_range.min(), self.lat_range.max()
                    ],
                    origin='lower',
                    aspect='auto',
                    cmap='Blues_r',
                    vmin=-50,
                    vmax=0
                )

                q = ax.quiver(
                    self.lonc,
                    self.latc,
                    u,
                    v,
                    color='k',
                    scale=4,
                    width=0.003
                )

                # Quiver key only once
                if k == 0:
                    box = patches.FancyBboxPatch(
                        (0.025, 0.945),
                        0.31,
                        0.02,
                        boxstyle="round,pad=0.01",
                        transform=ax.transAxes,
                        facecolor='white',
                        edgecolor='none',
                        alpha=0.5,
                        zorder=2
                    )
                    ax.add_patch(box)

                    qk = ax.quiverkey(
                        q,
                        X=0.15,
                        Y=0.95,
                        U=0.5,
                        label='0.5 m/s',
                        labelpos='E',
                        coordinates='axes'
                    )
                    qk.set_zorder(10)
                    qk.text.set_zorder(11)

                # Titles
                if k == 0:
                    ax.set_title("Surface Currents")
                elif k == 9:
                    ax.set_title("Bottom Currents")
                else:
                    ax.set_title(f"Level {k+1}")

                ax.set_xlim(self.lon_range.min(), self.lon_range.max())
                ax.set_ylim(self.lat_range.min(), self.lat_range.max())

            fig.suptitle(f"Dabob Bay Currents - {str(self.time[t0])}", fontsize=14)

            if save_path:
                plt.savefig(save_path, dpi=300, bbox_inches="tight")

            plt.show()

        # Animate over time
        else:
            end_t = np.datetime64(end_t)
            t1 = np.argmin(np.abs(self.time - end_t))

            fig, axes = plt.subplots(2, 5, figsize=(24, 10), constrained_layout=True)
            axes = axes.ravel()

            land = np.ma.masked_where(self.bty >= 0, self.bty)

            # draw static bathymetry
            for ax in axes:
                ax.imshow(
                    land,
                    extent=[
                        self.lon_range.min(), self.lon_range.max(),
                        self.lat_range.min(), self.lat_range.max()
                    ],
                    origin='lower',
                    aspect='auto',
                    cmap='Blues_r',
                    vmin=-50,
                    vmax=0
                )
                ax.set_xlim(self.lon_range.min(), self.lon_range.max())
                ax.set_ylim(self.lat_range.min(), self.lat_range.max())

            # initial quivers
            quivers = []
            for k in range(10):
                q = axes[k].quiver(
                    self.lonc,
                    self.latc,
                    self.u[k, :, t0],
                    self.v[k, :, t0],
                    color='k',
                    scale=4,
                    width=0.003
                )
                quivers.append(q)

            def update(frame):

                t_idx = t0 + frame
                current_time = self.time[t_idx]

                for k in range(10):
                    u = self.u[k, :, t_idx]
                    v = self.v[k, :, t_idx]
                    quivers[k].set_UVC(u, v)

                    if k == 0:
                        axes[k].set_title("Surface Currents")
                    elif k == 9:
                        axes[k].set_title("Bottom Currents")
                    else:
                        axes[k].set_title(f"Level {k+1}")

                fig.suptitle(f"Dabob Bay Currents - {str(current_time)}", fontsize=14)

                return quivers

            anim = animation.FuncAnimation(
                fig,
                update,
                frames=(t1 - t0),
                interval=100,
                blit=False
            )

            if save_path:
                anim.save(save_path, dpi=200, fps=8)

            plt.show()


# Run the Salish Sea Model
if __name__ == "__main__":
    # Dabob Bathymetry and Source/Receiver Coordinates
    ss_dir = os.path.join(os.getcwd(), 'data', 'SS_Data')
    bty_file_path =  os.path.join(os.getcwd(), 'data', 'bty.mat')
    source_coords = np.array([[-122.802983 , 47.77295]])
    receiver_coords = np.array([[-122.84, 47.73]])
    trackline = np.array([source_coords[0], receiver_coords[0]])
    bty_save = os.path.join(os.getcwd(), 'save', 'plots', 'bty.png')

    # Time Range
    start_t = datetime(2020, 1, 21, 0, 0, 0)
    end_t = datetime(2020, 1, 30, 0, 0, 0)

    # Salish Sea Model
    salish_sea = Salish_Sea(ss_dir, bty_file_path=bty_file_path)

    # Plot Bathymetry
    salish_sea.plot_bty(trackline=trackline, save_path=bty_save)

    # Plot Current
    salish_sea.plot_current(start_t=start_t, save_path=os.path.join(os.getcwd(), 'save', 'plots', '28_surface_curr.png'))