import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from matplotlib.animation import FuncAnimation
import numpy as np
import h5py
import scipy.io as io

class LiveOceanData:
    def __init__(self, filepath):
        self.filepath = filepath
        self.file = None
        self.lat = None
        self.lon = None
        self.h = None
        self.dabob_group = None
        self.temp_refs = None
        self.salt_refs = None
        self.zr_refs = None
        self.time_refs = None
        self.load_data()

    def load_data(self):
        self.file =h5py.File(self.filepath, 'r')
        self.dabob_group = self.file['DABOB']
        self.temp_refs = self.dabob_group["temp"][:]
        self.salt_refs = self.dabob_group["salt"][:]
        self.zr_refs = self.dabob_group["zr"][:]
        self.time = np.array([datetime(2020, 1, 1, 0, 0, 0) + i * timedelta(hours=1) for i in range(int((datetime(2020, 1, 31, 23, 0, 0) - datetime(2020, 1, 1, 0, 0, 0)).total_seconds() / 3600) + 1)])
        self.lat = np.squeeze(self.file['lat'][()])
        self.lon = np.squeeze(self.file['lon'][()])
        self.h = self.file['h'][()]

    def get_temp(self, index):
        ref = self.temp_refs[index]
        if isinstance(ref, np.ndarray):
            ref = ref.item() 
        temp_data = self.file[ref][()]
        temp = np.transpose(temp_data, (1, 2, 0))
        return np.array(temp)
    
    def get_salt(self, index):
        ref = self.salt_refs[index]
        if isinstance(ref, np.ndarray):
            ref = ref.item() 
        salt_data = self.file[ref][()]
        salt = np.transpose(salt_data, (1, 2, 0))
        return np.array(salt)
    
    def get_zr(self, index):
        ref = self.zr_refs[index]
        if isinstance(ref, np.ndarray):
            ref = ref.item() 
        zr_data = self.file[ref][()]
        zr = np.transpose(zr_data, (1, 2, 0))
        return np.array(zr)
    
    def ssp(self, depth, temp, salt):
        """Mackenzie sound speed formula."""
        return (
            1448.96
            + 4.591 * temp
            - 5.304e-2 * temp**2
            + 2.374e-4 * temp**3
            + 1.340 * (salt - 35)
            + 1.630e-2 * depth
            + 1.675e-7 * depth**2
            - 1.025e-2 * temp * (salt - 35)
            - 7.139e-13 * temp * depth**3
        )

    
    def animate_surface(self, save_path=None, fps=30):
        # Define spatial limits
        lon_min, lon_max = -123.05, -122.70
        lat_min, lat_max = 47.56, 47.85

        # Find index bounds within the domain
        lon_mask = (self.lon >= lon_min) & (self.lon <= lon_max)
        lat_mask = (self.lat >= lat_min) & (self.lat <= lat_max)
        lon_inds = np.where(lon_mask)[0]
        lat_inds = np.where(lat_mask)[0]

        lon_slice = slice(lon_inds.min(), lon_inds.max() + 1)
        lat_slice = slice(lat_inds.min(), lat_inds.max() + 1)

        # Setup figure (1 row × 4 columns)
        fig, axes = plt.subplots(1, 4, figsize=(18, 6 * 1.3), constrained_layout=True)
        ax_temp_surf, ax_temp_bot, ax_salt_surf, ax_salt_bot = axes
        fig.patch.set_facecolor("white")

        # Initial data
        temp0 = self.get_temp(0)[lat_slice, lon_slice, :]
        salt0 = self.get_salt(0)[lat_slice, lon_slice, :]

        surface_temp = temp0[:, :, 0]
        bottom_temp = temp0[:, :, -1]
        surface_salt = salt0[:, :, 0]
        bottom_salt = salt0[:, :, -1]

        # Shared vmin/vmax (non-normalized)
        tmin, tmax = 6.5, 10.5
        smin, smax = 26, 31.5

        # Plot initial frames
        im_temp_surf = ax_temp_surf.imshow(
            surface_temp, cmap='inferno', origin='lower',
            extent=[self.lon[lon_slice].min(), self.lon[lon_slice].max(),
                    self.lat[lat_slice].min(), self.lat[lat_slice].max()],
            aspect=1.3, vmin=tmin, vmax=tmax
        )
        im_temp_bot = ax_temp_bot.imshow(
            bottom_temp, cmap='inferno', origin='lower',
            extent=[self.lon[lon_slice].min(), self.lon[lon_slice].max(),
                    self.lat[lat_slice].min(), self.lat[lat_slice].max()],
            aspect=1.3, vmin=tmin, vmax=tmax
        )
        im_salt_surf = ax_salt_surf.imshow(
            surface_salt, cmap='viridis', origin='lower',
            extent=[self.lon[lon_slice].min(), self.lon[lon_slice].max(),
                    self.lat[lat_slice].min(), self.lat[lat_slice].max()],
            aspect=1.3, vmin=smin, vmax=smax
        )
        im_salt_bot = ax_salt_bot.imshow(
            bottom_salt, cmap='viridis', origin='lower',
            extent=[self.lon[lon_slice].min(), self.lon[lon_slice].max(),
                    self.lat[lat_slice].min(), self.lat[lat_slice].max()],
            aspect=1.3, vmin=smin, vmax=smax
        )

        # Titles and labels
        ax_temp_surf.set_title("Surface Temperature (°C)")
        ax_temp_bot.set_title("Bottom Temperature (°C)")
        ax_salt_surf.set_title("Surface Salinity (psu)")
        ax_salt_bot.set_title("Bottom Salinity (psu)")

        for ax in axes:
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")
            ax.set_xlim(lon_min, lon_max)
            ax.set_ylim(lat_min, lat_max)
            ax.tick_params(axis='x', rotation=45)
            ax.set_aspect(1.3)

        # Shared colorbars
        cbar_temp = fig.colorbar(im_temp_surf, ax=[ax_temp_surf, ax_temp_bot], shrink=0.8, location='bottom', pad=0.1)
        cbar_salt = fig.colorbar(im_salt_surf, ax=[ax_salt_surf, ax_salt_bot], shrink=0.8, location='bottom', pad=0.1)
        cbar_temp.set_label("Temperature (°C)")
        cbar_salt.set_label("Salinity (psu)")

        # Main title
        suptitle = fig.suptitle(f"LiveOcean Surface & Bottom Fields — {self.time[0].strftime('%Y-%m-%d %H:%M')}")

        # Update function for animation
        def update(frame):
            print(frame)
            temp = self.get_temp(frame)[lat_slice, lon_slice, :]
            salt = self.get_salt(frame)[lat_slice, lon_slice, :]

            im_temp_surf.set_data(temp[:, :, 0])
            im_temp_bot.set_data(temp[:, :, -1])
            im_salt_surf.set_data(salt[:, :, 0])
            im_salt_bot.set_data(salt[:, :, -1])

            suptitle.set_text(f"LiveOcean Surface & Bottom Fields — {self.time[frame].strftime('%Y-%m-%d %H:%M')}")
            return [im_temp_surf, im_temp_bot, im_salt_surf, im_salt_bot]

        # Animate
        anim = FuncAnimation(fig, update, frames=len(self.time), interval=1000/fps, blit=False)

        if save_path:
            print(f"Saving animation to {save_path}...")
            anim.save(save_path, writer='ffmpeg', fps=fps, dpi=200)
            print("✅ Animation saved successfully.")
        else:
            plt.show()

