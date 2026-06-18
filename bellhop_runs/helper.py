##########################################
### RAMPE and BELLHOP Helper Functions ###
##########################################
import numpy as np
import os
import subprocess
import time

CW_ARRAY = None
ZW_ARRAY = None
ARMS_TIME = None

def init_worker(cw_array_, zw_array_, arms_time_):
    global CW_ARRAY, ZW_ARRAY, ARMS_TIME

    CW_ARRAY = cw_array_
    ZW_ARRAY = zw_array_
    ARMS_TIME = arms_time_


def run_ram_case(args):
    # Input Arguements
    (
        t_idx,f_idx,freq,SD,RD,RR,dr,dz,ndr,ndz,
        rb,z_rb,za,ca,rhoa,attna,
        rhob,attnb,zmplt,c0,surface,x,ramrun_dir
    ) = args

    # Display Current Frequency
    print(f"Running Frequency: {freq:.1f} Hz\n")

    t_curr = ARMS_TIME[t_idx]
    cw = CW_ARRAY[:, :, t_idx] # SSP (m/s)
    zw = ZW_ARRAY[:, t_idx];   # Depth Array for SSP (m)
    attnw = cw * 0        # Attenuation (dB/lambda)
    rhow = cw * 0 + 1     # Density (g/cm^3)
    zb = zw[-1] + np.array([0, 400, 500, 600, 700])  # Bottom Depth (Layers (m))
    cb =1700 * np.ones_like(zb)                      # Bottom Sound Speed (Layer Sound Speeds)

    # Write RAMPE Input File for Current Time/Frequency Step
    rampe_input_writer(
        t_curr,
        freq,
        SD, RD, RR,
        dr, dz, ndr, ndz,
        rb, z_rb,
        zw, cw,
        rhow, attnw,
        za, ca, rhoa, attna,
        zb, cb, rhob, attnb,
        zmplt, c0,
        surface, x
    )

    # Current RAMPE Input/Output Directory
    ram_input_dir = os.path.join(
        os.getcwd(),
        "input_output",
        "RAMPE",
        t_curr.strftime("%Y-%m-%d_%H-%M-%S"),
        str(freq)
    )

    tic = time.time()

    result = subprocess.run(
        [ramrun_dir],
        cwd=ram_input_dir,
        capture_output=True,
        text=True
    )

    toc = time.time()

    return f_idx, freq, result, toc - tic


def run_bellhop_case(args):
    # Input Arguements
    (
        t_idx,f_idx,freq,
        nmedia,sspopt,rw,
        rb,z_rb,r_ati,z_ati,
        bot_type,rough,cb,ssb,
        rhob,attnb,NSD,SD,
        NRD,RD,NRR,RR,
        run_type,num_beam,min_ang,
        max_ang,step,bellhop_dir
    ) = args

    # Display Current Frequency
    print(f"Running Frequency: {freq:.1f} Hz\n")

    t_curr = ARMS_TIME[t_idx]
    cw = np.transpose(CW_ARRAY[:, :, t_idx]) # SSP (m/s)
    zw = ZW_ARRAY[:, t_idx];   # Depth Array for SSP (m)
    zw[-1] = 200.0
    # Write RAMPE Input File for Current Time/Frequency Step
    bellhop_input_writer(
        t_curr,freq,nmedia,
        sspopt,cw,zw,rw,
        rb,z_rb,r_ati,
        z_ati,bot_type,
        rough,cb,ssb,rhob,
        attnb,NSD,SD,NRD,RD,
        NRR,RR,run_type,
        num_beam,min_ang,
        max_ang,step
    )

    # Current RAMPE Input/Output Directory
    bellhop_input_dir = os.path.join(
        os.getcwd(),
        "input_output",
        "BELLHOP",
        t_curr.strftime("%Y-%m-%d_%H-%M-%S"),
        str(freq)
    )

    tic = time.time()

    result = subprocess.run(
        [bellhop_dir, "bellhop"],
        cwd=bellhop_input_dir,
        capture_output=True,
        text=True
    )

    print(f"\n===== {freq} Hz =====")
    print("STDOUT:")
    print(result.stdout)

    print("STDERR:")
    print(result.stderr)

    print("RETURN CODE:")
    print(result.returncode)

    toc = time.time()

    return f_idx, freq, result, toc - tic


def rampe_input_writer(t_curr,freq,SD,RD,RR,dr,dz,ndr,ndz,rb,z_rb,zw,cw,rhow,attnw,za,ca,rhoa,attna,zb,cb,rhob,attnb,zmplt,c0,surface,x):
    ###########################################
    # Writing Waveguide Properties input file #
    ###########################################
    ramin_dir = os.path.join(os.getcwd(), "input_output", "RAMPE", t_curr.strftime("%Y-%m-%d_%H-%M-%S"), str(freq))
    os.makedirs(ramin_dir, exist_ok=True)
    filename = 'ram.in'
    
    # Save Variable Names
    fs = freq
    zs = SD - za[0]
    zr = zs
    rmax = RR
    zmax = zb[-1] - za[0]
    npade = 8
    ns = 1
    rs = 0 
    k0 = 2 * np.pi * fs / c0
    
    isrc_type = 0  # Source Type (0 = Monopole, 1 = Beam)
    r0  = 0.0 
    r1  = 0.0  
    theta0 = 5              
    beam_width = 15


    with open(os.path.join(ramin_dir, filename), "w") as fid:
        ### Header ###
        fid.write("Test input file for range-dependent RAM_PE\n")
        fid.write(f"{isrc_type:d} {r0:6.6f} {r1:6.6f} {theta0:6.6f} {beam_width:6.6f}\n")
        fid.write(f"{fs:6.6f} {zs:6.6f} {zr:6.6f}\n")
        fid.write(f"{rmax:6.6f} {dr:6.6f} {ndr:6.0f}\n")
        fid.write(f"{zmax:6.6f} {dz:6.6f} {ndz:6.0f} {zmplt:6.6f}\n")
        fid.write(f"{c0:6.6f} {npade:6.0f} {ns:6.0f} {rs:6.0f}\n")
        
        ### Bathemetry ###
        RINDP = np.vstack((rb, z_rb))  # same as MATLAB [rb; z_rb]
        for r, z in zip(RINDP[0], RINDP[1]):
            fid.write(f"{r:6.6f} {z:6.6f}\n")
        fid.write(f"{-1:6.0f} {-1:6.0f}\n")
        
        ### Water Column Properties ###
        zasave = za
        casave = ca
        rhoasave = rhoa
        attnasave = attna
        zwsave = zw
        cwsave = np.transpose(cw)
        rhowsave = np.transpose(rhow)
        attnwsave = np.transpose(attnw)
        
        ### Define Water Column at each Surface Height Range ###
        for hdex in range(len(x)):

            # --- Reset to original profiles ---
            za = zasave.copy()
            ca = casave.copy()
            rhoa = rhoasave.copy()
            attna = attnasave.copy()
            zw = zwsave.copy()

            cw = cwsave[:, hdex].copy()
            rhow = rhowsave[:, hdex].copy()
            attnw = attnwsave[:, hdex].copy()

            # --- Surface adjustment ---
            h = -surface[hdex]
            za[-1] = za[-1] - h
            za = np.append(za, za[-1] + dz)   # avoid overlap at interface

            ca = np.append(ca, cw[0])
            rhoa = np.append(rhoa, rhow[0])
            attna = np.append(attna, attnw[0])

            # --- Remove overlapping water column points ---
            zdex = np.where(zw <= za[-1])[0]

            zw = np.delete(zw, zdex)
            cw = np.delete(cw, zdex)
            rhow = np.delete(rhow, zdex)
            attnw = np.delete(attnw, zdex)

            # =========================
            # Range-dependent block tag
            # =========================
            if hdex > 0:
                fid.write(f"{x[hdex]:7.2f}             rp\n")

            # =========================
            # AIR / WATER SSP
            # =========================
            z_aw = np.concatenate((za, zw)) - za[0]
            c_aw = np.concatenate((ca, cw))

            SSPFL_AW = np.vstack((z_aw, c_aw))
            for i in range(SSPFL_AW.shape[1]):
                fid.write(f"{SSPFL_AW[0, i]:6.6f} {SSPFL_AW[1, i]:6.6f}\n")
            fid.write("-1 -1\n")

            # =========================
            # AIR / WATER DENSITY
            # =========================
            rho_aw = np.concatenate((rhoa, rhow))

            DENPFL_AW = np.vstack((z_aw, rho_aw))
            for i in range(DENPFL_AW.shape[1]):
                fid.write(f"{DENPFL_AW[0, i]:6.6f} {DENPFL_AW[1, i]:6.6f}\n")
            fid.write("-1 -1\n")

            # =========================
            # AIR / WATER ATTENUATION
            # =========================
            attn_aw = np.concatenate((attna, attnw))

            ATTNPFL_AW = np.vstack((z_aw, attn_aw))
            for i in range(ATTNPFL_AW.shape[1]):
                fid.write(f"{ATTNPFL_AW[0, i]:6.6f} {ATTNPFL_AW[1, i]:6.6f}\n")
            fid.write("-1 -1\n")

            # =========================
            # SEABED SSP
            # =========================
            z_b = zb - za[0]

            SSPFL_B = np.vstack((z_b, cb))
            for i in range(SSPFL_B.shape[1]):
                fid.write(f"{SSPFL_B[0, i]:6.6f} {SSPFL_B[1, i]:6.6f}\n")
            fid.write("-1 -1\n")

            # =========================
            # SEABED DENSITY
            # =========================
            DENPFL_B = np.vstack((z_b, rhob))
            for i in range(DENPFL_B.shape[1]):
                fid.write(f"{DENPFL_B[0, i]:6.6f} {DENPFL_B[1, i]:6.6f}\n")
            fid.write("-1 -1\n")

            # =========================
            # SEABED ATTENUATION
            # =========================
            ATTNPFL_B = np.vstack((z_b, attnb))
            for i in range(ATTNPFL_B.shape[1]):
                fid.write(f"{ATTNPFL_B[0, i]:6.6f} {ATTNPFL_B[1, i]:6.6f}\n")
            fid.write("-1 -1\n")
        

def rampe_output_reader(t_curr, freq):
    ###########################################
    # Reading Waveguide Properties input file #
    ###########################################
    ramout_dir = os.path.join(os.getcwd(), "input_output", "RAMPE", t_curr.strftime("%Y-%m-%d_%H-%M-%S"), str(freq))
    filename = os.path.join(ramout_dir, "ram.in")

    with open(filename, "r") as fid_input:

        # Skip header lines
        tline = fid_input.readline()
        tline = fid_input.readline()

        # Parse First Line
        tline = fid_input.readline()
        val = np.fromstring(tline, sep=' ')

        freq = val[0]   # Source Frequency (Hz)
        zs   = val[1]
        zr   = val[2]

        # Parse Second Line
        tline = fid_input.readline()
        val = np.fromstring(tline, sep=' ')

        rmax = val[0]      # Maximum Range (m)
        dr   = val[1]      # Range Grid Spacing
        ndr  = int(val[2]) # Number of Receiver Outputs

        # Parse Third Line
        tline = fid_input.readline()
        val = np.fromstring(tline, sep=' ')

        zmax = val[0]      # Maximum Depth (m)
        dz   = val[1]      # Depth Grid Spacing
        ndz  = int(val[2]) # Number of Depth Outputs
        zmplt = val[3]

        # Parse Fourth Line
        tline = fid_input.readline()
        val = np.fromstring(tline, sep=' ')

        c0 = val[0]        # Reference Sound Speed (m/s)
        np_val = int(val[1])
        ns = int(val[2])
        rs = val[3]

        # Parse Fifth Line
        tline = fid_input.readline()
        val = np.fromstring(tline, sep=' ')

        rbs = val[0]
        zbs = val[1]

    #########################################
    # Convert to Operational Variables
    #########################################
    rge = rmax
    dep = zmplt
    k0 = 2 * np.pi * freq / c0
    rhow = 1000.0  # Assumption of Constant Impedance

    #########################################
    # Read tl.grid Binary File
    #########################################
    dat_filename = os.path.join(ramout_dir, "tl.grid")

    with open(dat_filename, "rb") as fid:
        PEdata1 = np.fromfile(fid, dtype=np.float32)
    PEdata = PEdata1.copy()

    #########################################
    # Remove Wrapped Values
    #########################################
    wrapdex = np.where(PEdata == PEdata[3])[0]
    PEdata = np.delete(PEdata, wrapdex)
    PEdata = PEdata[3:]

    #########################################
    # Determine Output Dimensions
    #########################################
    xsize = (wrapdex[1] - wrapdex[0]) - 1
    ysize = int(len(PEdata) / xsize / 2)

    #########################################
    # Reshape Complex Pressure Field
    #########################################
    OUT = PEdata.reshape((xsize, 2, ysize), order='F')
    OUT = OUT[:, :, 1:]
    OUT = OUT[:, 0, :] + 1j * OUT[:, 1, :]

    #########################################
    # Create Depth and Range Vectors
    #########################################
    Depth = np.arange(dz, xsize * ndz * dz, ndz * dz)
    Range = np.arange(ndr * dr,
                    dr * ndr * (ysize - 1) + ndr * dr,
                    ndr * dr)

    #########################################
    # Reallocate range exp(ik0r) phase dependence
    #########################################
    phase = np.exp(1j * k0 * Range - 1j * np.pi / 4)
    P = OUT * phase[np.newaxis, :]

    #########################################
    # Compute Pressure Gradients
    #########################################
    dPz = np.gradient(P, dz, axis=0)
    dPr = np.gradient(P, dr, axis=1)

    #########################################
    # Particle Velocity
    #########################################
    omega = k0 * c0
    Vz = -dPz / (-1j * omega * rhow)
    Vr = -dPr / (-1j * omega * rhow)

    #########################################
    # Acoustic Intensity
    #########################################
    Iz = np.real(P * np.conj(Vz))
    Ir = np.real(P * np.conj(Vr))

    # Return Range, Depth, and Pressure Arrays
    return Range, Depth, P


def bellhop_input_writer(t_curr,freq,nmedia,sspopt,cw,zw,rw,rb,z_rb,r_ati,z_ati,bot_type,rough,cb,ssb,rhob,attnb,NSD,SD,NRD,RD,NRR,RR,run_type,num_beam,min_ang,max_ang,step):

    ###################
    # Create ENV File #
    ###################
    bellhop_dir = os.path.join(os.getcwd(), "input_output", "BELLHOP", t_curr.strftime("%Y-%m-%d_%H-%M-%S"), str(freq))
    os.makedirs(bellhop_dir, exist_ok=True)
    filename = os.path.join(bellhop_dir, "bellhop.env")
    with open(filename, "w") as fid:

        ### Header ###
        fid.write('bellhop.env\t\t\t! TITLE\n')
        fid.write(f'{freq:.1f}\t\t\t! FREQ (Hz)\n')
        fid.write(f'{nmedia:.0f}\t\t\t! NMEDIA\n')
        fid.write(f'{sspopt}\t\t\t! SSPOPT\n')

        ### SSP Specification ###
        fid.write(f'{len(zw):.0f} {np.min(zw):.2f} {np.max(zw):.2f}\t\t\t! DEPTH of bottom (m)\n')
        SSP_DATA = np.vstack((zw, cw[:,1]))
        for i in range(SSP_DATA.shape[1]):
            fid.write(f'{SSP_DATA[0,i]:.2f} {SSP_DATA[1,i]:.2f} /\n')
        fid.write('\n')

        # Bottom Type, Roughness, and Bottom Properties
        fid.write(f'\'{bot_type}\' {rough:.1f} \t\t\t! BOTTOM Type, roughness \n')
        fid.write(f'{np.max(z_rb):.1f} {cb:.2f} {ssb:.1f} {rhob:.1f} {attnb:.1f} /\t\t\t! Bottom depth, compressional speed, shear speed, density, and attenuation \n')
        fid.write('\n')
        
        # Source Depths
        fid.write(f'{NSD:.0f}\t\t\t! NSD: Number of source depths \n')
        fid.write(f'{SD:.1f} /\t\t\t! Source Depth (m) \n')
        fid.write('\n')

        # Receiver Depths
        fid.write(f'{NRD:.0f}\t\t\t! NRD: Number of receiver depths \n')
        fid.write(' '.join(f'{z:.1f}' for z in RD) + ' /\t\t\t! Receiver depths (m) \n')
        fid.write('\n')

        # Receiver Ranges
        fid.write(f'{NRR:.0f}\t\t\t! NRR: Number of receiver ranges \n')
        fid.write(' '.join(f'{r/1000:.1f}' for r in RR) + ' /\t\t\t! Receiver ranges (km) \n')
        fid.write('\n')

        # Run Type and Angle Specs
        fid.write(f'\'{run_type}\'\t\t\t ! R: Ray Tracing, C: Coherent TL, I: Incoherent TL, S: Arrivals\n')
        fid.write(f'{num_beam:.0f} \t\t\t! Number of beams \n')
        fid.write(f'{min_ang:.1f} {max_ang:.1f} /\t\t\t! Launch Angles (degrees) \n')
        fid.write('\n')

        # Step Size, Max Depth, Max Range
        fid.write(f'{step:.1f} {np.max(zw):.2f} {np.max(rw/1000)-0.001:.3f} \t\t\t! Step Size (m), Max Depth (m), Max Range (km)\n')


    ###################
    # Create SSP File #
    ###################
    filename = os.path.join(bellhop_dir, "bellhop.ssp")
    with open(filename, "w") as fid:

        # Define SSP
        fid.write(f'{len(rw):.0f} \n');

        # Specify the Ranges
        for c in range(len(rw)):
            if c == len(rw) - 1:
                fid.write(f'{rw[c]/1000:.3f} \n')
            else:
                fid.write(f'{rw[c]/1000:.3f} ')
        # SSP at Range
        for z in range(len(zw)):
            for r in range(len(rw)):
                if r == len(rw) - 1:
                    fid.write(f'{cw[z, r]:.2f} \n')
                else:
                    fid.write(f'{cw[z, r]:.2f} ')

        fid.write('\n')


    ###################
    # Create BTY File #
    ###################
    filename = os.path.join(bellhop_dir, "bellhop.bty")
    with open(filename, "w") as fid:

        # Define BTY
        fid.write('\'L\'\n')
        fid.write(f'{len(rb):.0f},\n')
        BTY_DATA = np.vstack((rb/1000, z_rb))
        for i in range(BTY_DATA.shape[1]):
            fid.write(f'{BTY_DATA[0,i]:.3f} {BTY_DATA[1,i]:.3f} /\n')
        fid.write('\n')

                
    ###################
    # Create ATI File #
    ###################
    filename = os.path.join(bellhop_dir, "bellhop.ati")
    with open(filename, "w") as fid:

        # Define ATI
        fid.write('\'L\'\n')
        fid.write(f'{len(r_ati):.0f},\n')
        ATI_DATA = np.vstack((r_ati/1000, z_ati))
        for i in range(ATI_DATA.shape[1]):
            fid.write(f'{ATI_DATA[0,i]:.3f} {ATI_DATA[1,i]:.3f} /\n')
        fid.write('\n')


def bellhop_arr_output_reader(t_curr, freq):
    # Output File Name
    bellhop_dir = os.path.join(os.getcwd(), "input_output", "BELLHOP", t_curr.strftime("%Y-%m-%d_%H-%M-%S"), str(freq))
    filename = os.path.join(bellhop_dir, "bellhop.arr")

    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

        idx = 0

        dimension = lines[idx].replace("'", "")
        idx += 1

        freq = float(lines[idx])
        idx += 1

        # Source depths
        vals = lines[idx].split()
        nsd = int(vals[0])
        sd = np.array([float(v) for v in vals[1:]])
        idx += 1

        # Receiver depths
        vals = lines[idx].split()
        nrd = int(vals[0])
        rd = np.array([float(v) for v in vals[1:]])
        idx += 1

        # Receiver ranges
        vals = lines[idx].split()
        nrr = int(vals[0])
        rr = np.array([float(v) for v in vals[1:]])
        idx += 1

        # Max Number of Arrivals
        vals = lines[idx].split()
        max_arrivals = int(vals[0])
        idx += 1

        arrivals = []

        for ir in range(nrr):

            for idepth in range(nrd):

                narr = int(lines[idx])
                idx += 1

                depth_arrivals = []

                for _ in range(narr):

                    vals = lines[idx].split()
                    idx += 1

                    depth_arrivals.append({
                        "amp": float(vals[0]),
                        "phase_deg": float(vals[1]),
                        "time_s": float(vals[2]),
                        "imag_time": float(vals[3]),
                        "src_angle_deg": float(vals[4]),
                        "rcv_angle_deg": float(vals[5]),
                        "surface_bounces": int(vals[6]),
                        "bottom_bounces": int(vals[7])
                    })

                arrivals.append(depth_arrivals)

    return {
        "dimension": dimension,
        "frequency": freq,
        "source_depths": sd,
        "receiver_depths": rd,
        "receiver_ranges": rr,
        "arrivals": arrivals
    }


def bellhop_ctl_output_reader(t_curr, freq):
    # Output File Name
    bellhop_dir = os.path.join(os.getcwd(), "input_output", "BELLHOP", t_curr.strftime("%Y-%m-%d_%H-%M-%S"), str(freq))
    filename = os.path.join(bellhop_dir, "bellhop.shd")

    xs = np.nan
    ys = np.nan
        
    with open(filename, "rb") as fid:
        # =========================
        # Record 1
        # =========================
        recl = np.fromfile(fid, dtype=np.int32, count=1)[0]

        PlotTitle = fid.read(80).decode("ascii", errors="ignore").strip()

        fid.seek(4 * recl, 0)

        PlotType = fid.read(10).decode("ascii", errors="ignore").strip()

        # =========================
        # Record 2
        # =========================
        fid.seek(2 * 4 * recl, 0)

        Nfreq  = np.fromfile(fid, dtype=np.int32, count=1)[0]
        Ntheta = np.fromfile(fid, dtype=np.int32, count=1)[0]

        Nsx = np.fromfile(fid, dtype=np.int32, count=1)[0]
        Nsy = np.fromfile(fid, dtype=np.int32, count=1)[0]
        Nsz = np.fromfile(fid, dtype=np.int32, count=1)[0]

        Nrz = np.fromfile(fid, dtype=np.int32, count=1)[0]
        Nrr = np.fromfile(fid, dtype=np.int32, count=1)[0]

        freq0 = np.fromfile(fid, dtype=np.float64, count=1)[0]
        atten = np.fromfile(fid, dtype=np.float64, count=1)[0]

        # =========================
        # Record 3
        # =========================
        fid.seek(3 * 4 * recl, 0)

        freqVec = np.fromfile(fid, dtype=np.float64, count=Nfreq)

        # =========================
        # Record 4
        # =========================
        fid.seek(4 * 4 * recl, 0)

        theta = np.fromfile(fid, dtype=np.float64, count=Ntheta)

        # =========================
        # Record 5
        # =========================
        fid.seek(5 * 4 * recl, 0)

        sx = np.fromfile(fid, dtype=np.float64, count=2)
        sx = np.linspace(sx[0], sx[-1], Nsx)

        # =========================
        # Record 6
        # =========================
        fid.seek(6 * 4 * recl, 0)

        sy = np.fromfile(fid, dtype=np.float64, count=2)
        sy = np.linspace(sy[0], sy[-1], Nsy)

        # =========================
        # Record 7
        # =========================
        fid.seek(7 * 4 * recl, 0)

        sz = np.fromfile(fid, dtype=np.float32, count=Nsz)

        # =========================
        # Record 8
        # =========================
        fid.seek(8 * 4 * recl, 0)

        rz = np.fromfile(fid, dtype=np.float32, count=Nrz)

        # =========================
        # Record 9
        # =========================
        fid.seek(9 * 4 * recl, 0)

        rr = np.fromfile(fid, dtype=np.float64, count=Nrr)

        # =========================
        # Allocate pressure array
        # theta x sd x rd x rr
        # =========================
        pressure = np.zeros(
            (Ntheta, Nsz, Nrz, Nrr),
            dtype=np.complex64
        )

        # =========================
        # Frequency index
        # =========================
        freq = freqVec[0]

        freqdiff = np.abs(freqVec - freq)
        ifreq = np.argmin(freqdiff)

        # =========================
        # Read pressure records
        # =========================
        for itheta in range(Ntheta):

            for isz in range(Nsz):

                for irz in range(Nrz):

                    recnum = (
                        10
                        + ifreq * Ntheta * Nsz * Nrz
                        + itheta * Nsz * Nrz
                        + isz * Nrz
                        + irz
                    )

                    fid.seek(recnum * 4 * recl, 0)

                    temp = np.fromfile(
                        fid,
                        dtype=np.float32,
                        count=2 * Nrr
                    )

                    pressure[itheta, isz, irz, :] = (
                        temp[0::2] + 1j * temp[1::2]
                    )

    # =========================
    # Final outputs
    # =========================
    P = np.squeeze(pressure[0, 0, :, :])


    # Output Arrays (Same Output as RAMPE Output Reader)
    P = np.squeeze(pressure[itheta, isz, :, :])
    Depth = rz
    Range = rr


    return Range, Depth, P