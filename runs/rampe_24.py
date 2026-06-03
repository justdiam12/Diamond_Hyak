# #############################################################################
# This program runs RAM PE using input files LiveOcean SSP and Dabob Bay BTY.
#
# Originally Written by: David Dall'Osto (MATLAB)
# Modified by: Justin Diamond (Python)
#
# #############################################################################

########################
### Import Libraries ###
########################
from datetime import datetime
import time
import os
from concurrent.futures import ProcessPoolExecutor
import subprocess
import hdf5storage
import numpy as np
import scipy.io as io
import pandas as pd
from live_ocean import LiveOceanData
from helper import *
from scipy.interpolate import RegularGridInterpolator, interp1d
from pyproj import Geod
import matplotlib.pyplot as plt

if __name__ == "__main__":
    ########################
    ### File Directories ###
    ########################
    model_dir = os.path.join(os.getcwd(), "model")
    rampe_exe = os.path.join(model_dir, "ramrun")
    data_dir = os.path.join(os.getcwd(), "data")
    bty_file = os.path.join(data_dir, "bty.mat")
    liveocean_file = os.path.join(data_dir, "liveOcean2020_hourly.mat")
    save_file_dir = os.path.join(os.getcwd(), "save", "RAMPE")


    ##################
    ### Bathymetry ### 
    ##################
    # Source/Recvier Coordinates
    lon_start = -122.802983; 
    lon_end = -122.84;
    lat_start = 47.77295; 
    lat_end = 47.73;
    bty_data = io.loadmat(bty_file)
    lon_range = bty_data['lon_range'].flatten()
    lat_range = bty_data['lat_range'].flatten()
    bath_map = bty_data['bath_map']
    lon_track = np.linspace(lon_start, lon_end, 5000)
    lat_track = np.linspace(lat_start, lat_end, 5000)

    # Interpolate bathymetry along track
    interp_func = RegularGridInterpolator(
        (lat_range, lon_range),
        bath_map,
        bounds_error=False,
        fill_value=np.nan
    )

    track_points = np.column_stack((lat_track, lon_track))

    # Depth array (m)
    dl = -interp_func(track_points)

    # Compute cumulative geodesic range
    geod = Geod(ellps="WGS84")
    rl = np.zeros(5000)
    for i in range(1, 5000):
        _, _, dist = geod.inv(
            lon_track[i - 1],
            lat_track[i - 1],
            lon_track[i],
            lat_track[i]
        )
        rl[i] = rl[i - 1] + dist

    depth_interp = interp1d(
        rl,
        dl,
        bounds_error=False,
        fill_value="extrapolate"
    )

    # Final Range and Depth Arrays for Bathymetry
    range_m = np.arange(0, rl[-1], 1)
    depth = depth_interp(range_m)
    lat_line = np.interp(range_m, rl, lat_track)
    lon_line = np.interp(range_m, rl, lon_track)


    #######################
    ### Frequency Array ###
    #######################
    freqs = np.arange(3450, 3551) # Frequency Array (Hz)


    #################################
    ### Sound Speed Profile (SSP) ###
    #################################
    # Load LiveOcean Data
    live_ocean = LiveOceanData(liveocean_file)

    # Time Range for ARMS (January 24-30, 2020)
    arms_start = datetime(2020, 1, 24, 0, 0, 0)
    arms_end = datetime(2020, 1, 24, 23, 0, 0)
    arms_time = pd.date_range(start=arms_start, end=arms_end, freq='1h') # ARMS Time Range

    # Master SSP Array
    depth_size = 201
    zw_array = np.zeros((depth_size, len(arms_time)))             # SSP Depths (m)
    cw_array = np.zeros((len(range_m), depth_size, len(arms_time))) # SSP (m/s)

    print("----------------------------------------------------------\n")
    print("Collecting LiveOcean Sound Speed Profiles\n")
    print("----------------------------------------------------------\n")
    ramrun_dir = os.path.join(os.getcwd(), "models", "ramrun")

    for t in arms_time:
        # Time index for LiveOcean Data
        arms_t_idx = np.where(live_ocean.time == t)[0][0]

        # Extract Water Column Properties for Time Index
        temp = live_ocean.get_temp(arms_t_idx)
        salt = live_ocean.get_salt(arms_t_idx)
        zr = live_ocean.get_zr(arms_t_idx)

        depths = np.arange(temp.shape[2]) # Depth Array for Interpolation (m)

        zr_interp = RegularGridInterpolator((live_ocean.lat, live_ocean.lon, depths), zr)
        temp_interp = RegularGridInterpolator((live_ocean.lat, live_ocean.lon, depths), temp)
        salt_interp = RegularGridInterpolator((live_ocean.lat, live_ocean.lon, depths), salt)

        temp_line = np.zeros((len(range_m), len(depths)))
        salt_line = np.zeros((len(range_m), len(depths)))
        zr_line = np.zeros((len(range_m), len(depths)))

        for k, z in enumerate(depths):
            pts = np.column_stack((lat_line, lon_line, np.full(len(range_m), z)))
            temp_line[:, k] = temp_interp(pts)
            salt_line[:, k] = salt_interp(pts)
            zr_line[:, k] = zr_interp(pts)

        zr_positive = np.abs(zr_line)
        zr_positive = np.clip(zr_positive, 0, None)
        ssp_line = live_ocean.ssp(zr_positive, temp_line, salt_line)

        # Append to SSP Array
        depth_grid = np.linspace(np.nanmin(zr_positive), np.nanmax(zr_positive), depth_size)
        zw_array[:, arms_time.get_loc(t)] = depth_grid
   
        ssp_interp = np.zeros((len(range_m), depth_size))
        for r_idx in range(len(range_m)):

            z_profile = zr_positive[r_idx, :]
            c_profile = ssp_line[r_idx, :] 

            # Ensure monotonic depth ordering
            sort_idx = np.argsort(z_profile)

            z_sorted = z_profile[sort_idx]
            c_sorted = c_profile[sort_idx]

            # Remove duplicate depths if necessary
            z_unique, unique_idx = np.unique(z_sorted, return_index=True)
            c_unique = c_sorted[unique_idx]

            ssp_interp[r_idx, :] = np.interp(
                depth_grid,
                z_unique,
                c_unique,
                left=c_unique[0],
                right=c_unique[-1]
            )
        cw_array[:, :, arms_time.get_loc(t)] = ssp_interp

    # Load Rough Surface
    x = range_m                                   # Rough Surface Range (m)
    surface = x * 0 + np.sin(x*2*np.pi/1000) * 0  # Rough Surface Profile (m)


    ######################################
    ### Time Loop for RAMPE Simulation ###
    ######################################
    for t_idx in range(len(arms_time)):
        print("----------------------------------------------------------\n")
        print(f"Running at {arms_time[t_idx]:%Y-%m-%d %H:%M:%S}\n")
        print("----------------------------------------------------------\n")

        t_curr = arms_time[t_idx] # Current Time for RAMPE Simulation
        c0 = 1500             # Reference Sound Speed (m/s)
        SD = depth[0] - 4     # Source Depth (4m above surface)
        RD = 55               # Receiver Depth for PEoutput.line
        RR = range_m[-1]      # Receiver Range (m)
        dz =  0.05            # Depth Grid Spacing
        dr =  0.05            # Range Grid Spacing 
        ndr = 20              # Number of Range Outputs
        ndz = 20              # Number of Depth Outputs

        # Range and Bathymetry Properties
        rw = range_m               # Range Array for SSP (m)
        rhob = np.array([1, 1, 1, 1, 1]) * 1.50          # Bottom Density (g/cm^3)
        attnb = 0.2 * np.array([1, 1, 2, 5, 100])        # Bottom Attenuation (dB/lambda)
        zmplt = 600                                      # Maximum Grid Output Depth

        # Surficial Air Layer
        za =   np.array([150, 20, 15, 10, dz]) * -1      # Air Depth Layers (m)                           
        ca =    np.array([340, 340, 340, 340, 340])      # Air Sound Speed (m/s)
        rhoa =  np.array([1, 1, 1, 1, 1]) * 0.00128      # Air Density (g/cm^3)
        attna = np.array([100, 10, 1, 0, 0])             # Air Attenuation (dB/lambda)
        
        # Bathymetry for model input
        rb = range_m                      # Bathymetry range (m)
        z_rb = depth - za[0] + dz             # Bathymetry depth (m)

        nproc = os.cpu_count()  # or set manually (e.g., 16)

        tasks = []

        for f_idx in range(len(freqs)):
            tasks.append((
                t_idx,
                f_idx,
                freqs[f_idx],
                SD, RD, RR,
                dr, dz, ndr, ndz,
                rb, z_rb,
                za, ca, rhoa, attna,
                rhob, attnb,
                zmplt, c0,
                surface, x,
                ramrun_dir
            ))

        Psv = None

        with ProcessPoolExecutor(
            max_workers=nproc,
            initializer=init_worker,
            initargs=(cw_array, zw_array, arms_time)
        ) as ex:
            results = list(ex.map(run_ram_case, tasks))

        # sort results back into frequency order
        results.sort(key=lambda x: x[0])

        # initialize once
        Range, Depth, _ = rampe_output_reader(t_curr, freqs[0])
        Psv = np.zeros((len(Depth), len(Range), len(freqs)), dtype=np.complex64)

        for f_idx, freq, result, dt in results:
            Range, Depth, P = rampe_output_reader(t_curr, freq)
            Psv[:, :, f_idx] = P

            print(f"freq {freq} Hz finished in {dt:.2f} s")
            print(result.returncode)

        ### Save Output Files ###
        os.makedirs(save_file_dir, exist_ok=True)
        save_file = os.path.join(save_file_dir, f"rampe_time_{arms_time[t_idx]:%Y%m%d_%H%M%S}.mat")

        # Save HDF5 / MATLAB v7.3-compatible file
        hdf5storage.savemat(
            save_file,
            {
                "Psv": Psv,
                "Range": Range,
                "Depth": Depth
            },
            format='7.3'
        )

    # Save frequency/time file
    freq_time_file = os.path.join(save_file_dir, "rampe_freq_time_24.mat")
    hdf5storage.savemat(
        freq_time_file,
        {
            "freqs": freqs,
            "arms_time": np.array(arms_time.astype(str)),
            "cw_array": cw_array,
            "zw_array": zw_array
        },
        format='7.3'
    )