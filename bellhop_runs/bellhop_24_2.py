# #############################################################################
# This program runs BELLHOP using input files LiveOcean SSP and Dabob Bay BTY.
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
    bellhop_exe = os.path.join(model_dir, "bellhop")
    data_dir = os.path.join(os.getcwd(), "data")
    bty_file = os.path.join(data_dir, "bty.mat")
    liveocean_file = os.path.join(data_dir, "liveOcean2020_hourly.mat")
    arms_depth_dir = os.path.join(data_dir, "array_depths")
    save_file_dir = os.path.join(os.getcwd(), "save", "BELLHOP")


    ############################
    ### Extract Array Depths ###
    ############################
    arms_24 = io.loadmat(os.path.join(arms_depth_dir, "depths_20200124.mat"))
    arms_24_depth = arms_24['unique_depths'].flatten()
    arms_28 = io.loadmat(os.path.join(arms_depth_dir, "depths_20200128.mat"))
    arms_28_depth = arms_28['unique_depths'].flatten()
    arms_29 = io.loadmat(os.path.join(arms_depth_dir, "depths_20200129.mat"))
    arms_29_depth = arms_29['unique_depths'].flatten()
    arms_depths = np.sort(np.concatenate([arms_24_depth, arms_28_depth, arms_29_depth]))


    ##################
    ### Bathymetry ### 
    ##################
    # Source/Recvier Coordinates
    lon_start = -122.802983 
    lon_end = -122.84
    lat_start = 47.77295
    lat_end = 47.73
    receiver_range = 5500
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
    arms_start = datetime(2020, 1, 24, 1, 0, 0)
    arms_end = datetime(2020, 1, 24, 2, 0, 0)
    arms_time = pd.date_range(start=arms_start, end=arms_end, freq='1h') # ARMS Time Range

    # Master SSP Array
    depth_size = 201
    zw_array = np.zeros((depth_size, len(arms_time)))             # SSP Depths (m)
    cw_array = np.zeros((len(range_m), depth_size, len(arms_time))) # SSP (m/s)

    print("----------------------------------------------------------\n")
    print("Collecting LiveOcean Sound Speed Profiles\n")
    print("----------------------------------------------------------\n")
    bellhop_dir = os.path.join(os.getcwd(), "models", "bellhop.exe")

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

        zr_positive = zr_line * -1
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

    # Interpolate SSP and Depths for Minute Resolution
    arms_time_min = pd.date_range(start=arms_start, end=arms_end, freq='1min')
    t_hour = (arms_time - arms_time[0]).total_seconds().to_numpy()
    t_min = (arms_time_min - arms_time_min[0]).total_seconds().to_numpy()

    # Interpolation Functions for SSP and Depths
    cw_interp = interp1d(
        t_hour,
        cw_array,
        axis=2,
        kind='linear',
        bounds_error=False,
        fill_value='extrapolate'
    )
    zw_interp = interp1d(
        t_hour,
        zw_array,
        axis=1,
        kind='linear',
        bounds_error=False,
        fill_value='extrapolate'
    )

    # Minute-Resolution SSP, Depth Arrays and Time Array for ARMS
    cw_array = cw_interp(t_min)
    zw_array = zw_interp(t_min)
    arms_time = arms_time_min

    # Load Rough Surface
    r_ati = range_m                                   # Rough Surface Range (m)
    z_ati = r_ati * 0 + np.sin(r_ati*2*np.pi/1000) * 0  # Rough Surface Profile (m)


    ########################################
    ### Time Loop for BELLHOP Simulation ###
    ########################################
    for t_idx in range(len(arms_time)):
        print("----------------------------------------------------------\n")
        print(f"Running at {arms_time[t_idx]:%Y-%m-%d %H:%M:%S}\n")
        print("----------------------------------------------------------\n")

        t_curr = arms_time[t_idx] # Current Time for BELLHOP Simulation
        rw = range_m               # Range Array for SSP (m)

        # SSP Options
        nmedia = 1 # Arbitrary for Bellhop
        sspopt_1 = 'Q' # Cubic Spline Interpolation for SSP
        sspopt_2 = 'V' # Acoustic Half-Space at Surface
        sspopt_3 = 'F' # Thorpe Attenuation for the water
        sspopt_4 = ' ' # Default Parameter
        sspopt_5 = '*' # Includes .ati File for altimetry
        sspopt = sspopt_1 + sspopt_2 + sspopt_3 + sspopt_4 + sspopt_5

        # Bottom Type
        bot_type_1 = 'A'
        bot_type_2 = '*'
        bot_type = bot_type_1 + bot_type_2
        rough = 0.0

        # Bottom Properties
        cb = 1700         # Bottom Sound Speed (m/s)
        ssb = 0.0         # Bottom Shear Speed (m/s)
        rhob = 1.50       # Bottom Density (g/cm^3)
        attnb = 0.0       # Bottom Attenuation (dB/lambda)
        

        # Surficial Air Layer                         
        ca = 340          # Air Sound Speed (m/s)
        rhoa = 0.00128    # Air Density (g/cm^3)
        attna = 100       # Air Attenuation (dB/lambda)
        
        # Bathymetry for model input
        rb = range_m      # Bathymetry range (m)
        z_rb = depth      # Bathymetry depth (m)

        # Source Depths
        NSD = 1             # Number of Source Depths
        SD = z_rb[1] - 4    # Source Depth (m)

        # Receiver Depths and Ranges
        NRD = len(arms_depths)  # Number of Receiver Depths
        RD = arms_depths        # Receiver Depths (m)
        NRR = 1                 # Number of Receiver Ranges
        RR = [5500]               # Receiver Ranges (m)

        # Run Specifications
        run_type = 'A'
        angles = np.arange(-89, 89.1, 0.1)
        num_beam = len(angles)
        min_ang = np.min(angles)
        max_ang = np.max(angles)
        step = 0.1

        nproc = 8  # or set manually (e.g., 16)
        # nproc = os.cpu_count()  # or set manually (e.g., 16)

        tasks = []

        for f_idx in range(len(freqs)):
            tasks.append((
                t_idx,
                f_idx,
                freqs[f_idx],
                nmedia,
                sspopt,
                range_m,
                rb,
                z_rb,
                r_ati,
                z_ati,
                bot_type,
                rough,
                cb,
                ssb,
                rhob,
                attnb,
                NSD,
                SD,
                NRD,
                RD,
                NRR,
                RR,
                run_type,
                num_beam,
                min_ang,
                max_ang,
                step,
                bellhop_dir
            ))

        Psv = None

        with ProcessPoolExecutor(
            max_workers=nproc,
            initializer=init_worker,
            initargs=(cw_array, zw_array, arms_time)
        ) as ex:
            results = list(ex.map(run_bellhop_case, tasks))


        # sort results back into frequency order
        results.sort(key=lambda x: x[0])

        # initialize once
        arrivals = bellhop_arr_output_reader(t_curr, freqs[0])
        arr_save = np.empty((len(freqs), len(arms_depths)), dtype=object) 

        for f_idx, freq, result, dt in results:
            arrivals = bellhop_arr_output_reader(t_curr, freq)
            
            arrivals = arrivals['arrivals']

            for a in range(len(arrivals)):
                amp = np.array([])
                amp_time = np.array([])
                phase_deg = np.array([])
                curr = arrivals[a]                
                for b in range(len(curr)):
                    amp = np.append(amp, curr[b]['amp'])
                    amp_time = np.append(amp_time, curr[b]['time_s'])
                    phase_deg = np.append(phase_deg, curr[b]['phase_deg'])
                
                arr_save[f_idx, a] = {'amp': amp, 'time': amp_time, 'phase': phase_deg}
            
            print(f"freq {freq} Hz finished in {dt:.2f} s")
            print(result.returncode)

        ### Save Output Files ###
        os.makedirs(save_file_dir, exist_ok=True)
        save_file = os.path.join(save_file_dir, f"bellhop_time_{arms_time[t_idx]:%Y%m%d_%H%M%S}.mat")

        # Save HDF5 / MATLAB v7.3-compatible file
        hdf5storage.savemat(
            save_file,
            {
                "arrivals": arr_save,
                "freqs": freqs,
                "arms_depths": arms_depths,
                "arms_time": arms_time[t_idx]
            },
            format='7.3'
        )

        # Remove BELLHOP input/output files to save space
        bellhop_input_dir = os.path.join(
            os.getcwd(),
            "input_output",
            "BELLHOP",
            t_curr.strftime("%Y-%m-%d_%H-%M-%S"),
        )

        clean_dir = "rm -rf " + bellhop_input_dir
        subprocess.run(clean_dir, shell=True)