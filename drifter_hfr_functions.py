"""
drifter_hfr_functions.py
Drifter-HFR Comparison Functions Module
Contains all processing and plotting functions for comparing drifter data 
with HF Radar radial and total velocities.
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime, timedelta
import pyproj
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
from windrose import WindroseAxes
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cool_maps.plot as cplt
from oceans.ocfis import uv2spdir, spdir2uv
from hfradarpy.radials import Radial
import cmocean.cm as cmo
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # Convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    r = 6371000  # Radius of earth in meters
    return c * r

def calculate_bearing_distance(lat1, lon1, lat2, lon2):
    """
    Calculate bearing and distance between two points using pyproj
    
    Args:
        lat1, lon1: Latitude and longitude of point 1 (degrees)
        lat2, lon2: Latitude and longitude of point 2 (degrees)
        
    Returns:
        bearing: Forward azimuth in degrees (0° at North, clockwise positive)
        distance: Distance in meters
    """
    geodesic = pyproj.Geod(ellps='WGS84')
    fwd_azimuth, back_azimuth, distance = geodesic.inv(lon1, lat1, lon2, lat2)
    
    # Normalize to 0-360 range
    bearing = fwd_azimuth % 360
    
    return bearing, distance

def get_site_bearing_distance(site_lat, site_lon, drifter_lat, drifter_lon):
    """
    Calculate bearing and distance from a site to a drifter position
    """
    geodesic = pyproj.Geod(ellps='WGS84')
    fwd_azimuth, back_azimuth, distance = geodesic.inv(
        site_lon, site_lat,
        drifter_lon, drifter_lat
    )
    
    # Normalize to 0-360 range
    bearing = fwd_azimuth % 360
    
    return bearing, distance

def get_time_subset(drifter_ds, days=3):
    """Get a time subset from the middle of the drifter time series"""
    drifter_start = pd.to_datetime(drifter_ds.time.values[0])
    drifter_end = pd.to_datetime(drifter_ds.time.values[-1])
    t_middle = drifter_start + 0.5 * (drifter_end - drifter_start)
    half = timedelta(days=days/2)
    return t_middle - half, t_middle + half

# ============================================================================
# DRIFTER DATA PROCESSING
# ============================================================================

def import_drifter(file):
    """Import drifter data from a file."""
    ds = xr.open_dataset(file)
    return ds

def add_drifter_velocity(drifter_ds):
    """Calculate velocity and direction from drifter positions"""
    # Extract coordinates and times
    lats = drifter_ds.lat.values
    lons = drifter_ds.lon.values
    times = drifter_ds.time.values
    
    # Calculate time differences in seconds
    time_diffs = np.diff(times).astype('timedelta64[s]').astype(float)
    
    # Calculate bearings and distances between consecutive points
    bearings = []
    distances = []
    for i in range(len(lats)-1):
        bearing, distance = calculate_bearing_distance(
            lats[i], lons[i], 
            lats[i+1], lons[i+1]
        )
        bearings.append(bearing)
        distances.append(distance)
    
    # Convert to arrays
    bearings = np.array(bearings)
    distances = np.array(distances)
    
    # Calculate velocities (m/s)
    velocities = distances / time_diffs
    
    # Calculate u and v components
    u = velocities * np.sin(np.radians(bearings))  # Eastward velocity
    v = velocities * np.cos(np.radians(bearings))  # Northward velocity
    
    # Add calculated values to dataset with padding for the last point
    drifter_ds['velocity'] = (('time'), np.pad(velocities, (0, 1), 'edge'))
    drifter_ds['bearing'] = (('time'), np.pad(bearings, (0, 1), 'edge'))
    drifter_ds['u'] = (('time'), np.pad(u, (0, 1), 'edge'))
    drifter_ds['v'] = (('time'), np.pad(v, (0, 1), 'edge'))
    
    return drifter_ds

# ============================================================================
# QUALITY CONTROL FUNCTIONS
# ============================================================================

def apply_velocity_qc_roc(drifter_ds, threshold=0.5):
    """
    Apply QC by removing velocities with excessive rate of change.
    
    Arguments:
    - drifter_ds: xarray Dataset with calculated velocities.
    - threshold: float, maximum allowed rate of change in velocity (m/s).
    
    Returns:
    - cleaned_drifter_ds: xarray Dataset with anomalous velocities set to NaN.
    """
    # Create a copy of the dataset to avoid modifying the original
    cleaned_drifter_ds = drifter_ds.copy()
    
    # Calculate the rate of change in velocity (difference between consecutive values)
    velocity_diff = np.abs(np.diff(cleaned_drifter_ds['velocity'].values))
    
    # Identify points where the rate of change exceeds the threshold
    bad_points = np.insert(velocity_diff > threshold, 0, False)
    
    # Apply NaN to anomalous values in velocity, u, and v components
    cleaned_drifter_ds['velocity'].values[bad_points] = np.nan
    cleaned_drifter_ds['u'].values[bad_points] = np.nan
    cleaned_drifter_ds['v'].values[bad_points] = np.nan
    
    return cleaned_drifter_ds

def apply_velocity_qc_std(drifter_ds, std_levels=3):
    """
    Apply QC by removing velocities with excessive deviations from the mean.
    
    Arguments:
    - drifter_ds: xarray Dataset with calculated velocities.
    - std_levels: int, number of standard deviations to use as threshold.
    
    Returns:
    - cleaned_drifter_ds: xarray Dataset with anomalous velocities set to NaN.
    """
    # Create a copy of the dataset to avoid modifying the original
    cleaned_drifter_ds = drifter_ds.copy()
    
    # Calculate the mean and standard deviation of the velocity
    velocity_mean = np.nanmean(cleaned_drifter_ds['velocity'].values)
    velocity_std = np.nanstd(cleaned_drifter_ds['velocity'].values)
    
    # Ensure velocity_std is not zero to avoid division errors
    if velocity_std == 0:
        print("Warning: Zero standard deviation detected. No QC applied.")
        return cleaned_drifter_ds
    
    # Identify points where the velocity exceeds the threshold
    bad_points = np.abs(cleaned_drifter_ds['velocity'].values - velocity_mean) > (std_levels * velocity_std)
    
    # Apply NaN to anomalous values in velocity, u, and v components
    cleaned_drifter_ds['velocity'] = cleaned_drifter_ds['velocity'].where(~bad_points, np.nan)
    cleaned_drifter_ds['u'] = cleaned_drifter_ds['u'].where(~bad_points, np.nan)
    cleaned_drifter_ds['v'] = cleaned_drifter_ds['v'].where(~bad_points, np.nan)
    
    return cleaned_drifter_ds

def apply_distance_qc(drifter_ds, dist_threshold=10):
    """
    Apply QC by removing data points with distances between consecutive positions exceeding a threshold.
    
    Arguments:
    - drifter_ds: xarray Dataset with calculated velocities.
    - dist_threshold: float, maximum allowed distance between consecutive positions (meters).
    
    Returns:
    - cleaned_drifter_ds: xarray Dataset with anomalous data set to NaN.
    """
    # Create a copy of the dataset to avoid modifying the original
    cleaned_drifter_ds = drifter_ds.copy()
    
    # Extract coordinates
    lats = cleaned_drifter_ds.lat.values
    lons = cleaned_drifter_ds.lon.values
    
    # Calculate distances between consecutive points
    distances = []
    for i in range(len(lats) - 1):
        _, distance = calculate_bearing_distance(lats[i], lons[i], lats[i+1], lons[i+1])
        distances.append(distance)
    distances = np.array(distances)
    
    # Identify points where the distance exceeds the threshold
    bad_points = distances > dist_threshold
    # Pad bad_points to match the length of the dataset
    bad_points = np.insert(bad_points, -1, False)  # Add False at the end to match length
    
    # Apply NaN to anomalous values in velocity, u, and v components
    cleaned_drifter_ds['velocity'] = cleaned_drifter_ds['velocity'].where(~bad_points, np.nan)
    cleaned_drifter_ds['u'] = cleaned_drifter_ds['u'].where(~bad_points, np.nan)
    cleaned_drifter_ds['v'] = cleaned_drifter_ds['v'].where(~bad_points, np.nan)
    
    return cleaned_drifter_ds

def load_drifter(file, qc_params):
    """
    Load drifter data and apply quality control.
    
    Arguments:
    - file: str, path to the drifter data file.
    - qc_params: dict, QC parameters with keys:
        - 'roc_threshold': rate of change threshold (m/s)
        - 'std_levels': standard deviation levels
        - 'dist_threshold': distance threshold (m)
        - 'apply_roc': bool, whether to apply ROC QC
        - 'apply_std': bool, whether to apply STD QC
        - 'apply_dist': bool, whether to apply distance QC
    
    Returns:
    - cleaned_drifter_ds: xarray Dataset with applied QC.
    """
    ds = import_drifter(file)
    ds = add_drifter_velocity(ds)

    if qc_params.get('apply_roc', True):
        ds = apply_velocity_qc_roc(ds, threshold=qc_params.get('roc_threshold', 0.15))
    if qc_params.get('apply_std', True):
        ds = apply_velocity_qc_std(ds, std_levels=qc_params.get('std_levels', 3))
    if qc_params.get('apply_dist', True):
        ds = apply_distance_qc(ds, dist_threshold=qc_params.get('dist_threshold', 10000))

    return ds

# ============================================================================
# RADIAL DATA PROCESSING
# ============================================================================

def get_freq_resolution_for_site(site, site_frequency_dict, freq_resolution):
    """Get frequency resolution for a specific site"""
    for freq, sites in site_frequency_dict.items():
        if site in sites:
            return freq_resolution[freq]
    return None

def radial_velocity_to_site(site_lat, site_lon, drifter_lat_arr, drifter_lon_arr, u_arr, v_arr):
    """
    Calculate radial velocity from drifter to site
    
    Parameters
    ----------
    site_lat, site_lon : float
        Radar site coordinates (deg).
    drifter_lat_arr, drifter_lon_arr : 1-D arrays
        Drifter positions (deg) – same length as u_arr / v_arr.
    u_arr, v_arr : 1-D arrays
        Drifter east- and north-ward velocities (m s-1).

    Returns
    -------
    vr : 1-D np.ndarray
        Radial velocity (m s-1) positive AWAY from the radar.
    """
    geod = pyproj.Geod(ellps="WGS84")
    # Forward azimuth from radar to drifter (deg, CW from north)
    fwd_az, back_az, _ = geod.inv(
        np.full_like(drifter_lon_arr, site_lon),
        np.full_like(drifter_lat_arr, site_lat),
        drifter_lon_arr,
        drifter_lat_arr,
    )
    theta = np.deg2rad(back_az)  # to radians

    # LOS unit vector components
    ux = np.sin(theta)   # east component of LOS
    uy = np.cos(theta)   # north component of LOS

    # Dot product: (u,v) · (ux,uy)
    vr = u_arr * ux + v_arr * uy
    return vr

def get_nearest_radial(drifter_ds, site_frequency_dict, freq_resolution, 
                       site_loc_df, radial_src_dir, sites_to_process='all', 
                       proximity=None):
    """
    Get nearest radial data for each site near the drifter track
    
    Arguments:
    - drifter_ds: xarray Dataset with drifter data
    - site_frequency_dict: dict mapping frequency to site lists
    - freq_resolution: dict mapping frequency to distance resolution
    - site_loc_df: DataFrame with site locations
    - radial_src_dir: str, base directory for radial files
    - sites_to_process: str or list, 'all' or specific sites to process
    - proximity: float, distance in km to search for radial data (None = use default)
    
    Returns:
    - point_selected_site_radial_datasets: dict of xarray Datasets with radial data
    - full_site_radial_lists: dict of file paths used
    """
    geodesic = pyproj.Geod(ellps='WGS84')

    lats = drifter_ds['lat'].values
    lons = drifter_ds['lon'].values
    times = pd.to_datetime(drifter_ds['time'].values)

    filename_prefix = 'RDL'
    type_of_radial = {'ideal': 'i', 'measured': 'm'}
    suffix_rad = '.ruv'

    # Determine which sites to process
    if sites_to_process == 'all':
        sites = []
        for freq_sites in site_frequency_dict.values():
            if isinstance(freq_sites, list):
                sites.extend(freq_sites)
    elif isinstance(sites_to_process, str) and sites_to_process in site_frequency_dict:
        sites = site_frequency_dict[sites_to_process]
    elif isinstance(sites_to_process, list):
        sites = sites_to_process
    else:
        raise ValueError(f"Invalid sites_to_process: {sites_to_process}")

    # Remove duplicates
    sites = list(set(sites))

    freq_max_distance_to_site = {
        'long': 175,   # km
        'medium': 110,  # km
        'short': 70    # km
    }
    
    point_selected_site_radial_datasets = {}
    full_site_radial_lists = {}
    
    for site in sites:
        # Get site frequency for resolution
        site_freq = None
        for freq, freq_sites in site_frequency_dict.items():
            if site in freq_sites:
                site_freq = freq
                break
        
        if site_freq is None:
            print(f"WARNING: Could not determine frequency for site {site}")
            continue
            
        distance_to_vector_limit = freq_resolution[site_freq]
        
        if proximity is None:
            proximity_km = freq_max_distance_to_site[site_freq]
        else:
            proximity_km = proximity
            
        print(f"\n=== Processing site: {site} (freq: {site_freq}) ===")
        print(f"Proximity to site requirement: {proximity_km} km")
        print(f"Distance to vector limit: {distance_to_vector_limit} km")
        
        site_info = site_loc_df.loc[site_loc_df['site'] == site]
        if site_info.empty:
            print(f"WARNING: Site {site} not found in site_loc_df")
            continue

        site_lat = site_info['lat'].values[0]
        site_lon = site_info['lon'].values[0]
        print(f"Site coordinates: lat {site_lat}, lon {site_lon}")

        # Calculate distances from drifter positions to site
        fwd_az, back_az, distances = geodesic.inv(
            np.full_like(lons, site_lon),
            np.full_like(lats, site_lat),
            lons, lats
        )
        distances_km = distances / 1000.0
        within_range_times = times[distances_km <= proximity_km]
        print(f"Found {len(within_range_times)} times within {proximity_km} km of {site}")
        if len(within_range_times) == 0:
            print(f"No times within {proximity_km}km for {site}")
            continue

        time_list = []
        velocity_list = []
        u_list = []
        v_list = []
        heading_list = []
        radial_files = []
        lon_list = []
        lat_list = []
        distance_to_vector_list = []
        
        for dtime in within_range_times:
            dpos = drifter_ds.sel(time=dtime)
            drifter_lat = float(dpos.lat.values)
            drifter_lon = float(dpos.lon.values)

            # Build filename
            time_fmt = dtime.strftime('%Y_%m_%d_%H%M')
            year_month = dtime.strftime('%Y_%m')
            radial_file = f'{radial_src_dir}/{site}/{year_month}/{filename_prefix}{type_of_radial["measured"]}_{site}_{time_fmt}{suffix_rad}'
            # print(f"Checking for radial file: {radial_file}")
            radial_files.append(radial_file)
            
            if not os.path.exists(radial_file):
                continue

            try:
                # print(f"Reading radial file: {radial_file}")
                r = Radial(radial_file)
                if r.data.empty:
                    continue
                radial_data = r.data

                # Apply QC filtering
                if 'VFLG' in radial_data.columns:
                    radial_data = radial_data.query('VFLG == 0')
                if 'VELO' in radial_data.columns:
                    radial_data = radial_data.query('VELO < 32767')

                if radial_data.empty:
                    continue

                vec_lons = radial_data['LOND'].values
                vec_lats = radial_data['LATD'].values

                # Find nearest vector
                distances_vec = [haversine(drifter_lat, drifter_lon, vlat, vlon)
                                 for vlat, vlon in zip(vec_lats, vec_lons)]
                distances_vec = np.array(distances_vec)

                if np.all(distances_vec > distance_to_vector_limit*1000):
                    print(f"Max distances_vec: {np.nanmax(distances_vec)}, Min distances_vec: {np.nanmin(distances_vec)}")
                    # Save distances vector limit max and min info for debugging
                    print(f"Distance to vector limit: {distance_to_vector_limit*1000} m")
                    print(f"No vectors within {distance_to_vector_limit} km at time {dtime} for site {site}")
                    continue
                
                if len(distances_vec) == 0:
                    print(f"No vectors found at time {dtime} for site {site}")
                    continue

                idx_min = np.nanargmin(distances_vec)
                if np.isnan(distances_vec[idx_min]):
                    print(f"All distances are NaN at time {dtime} for site {site}")
                    continue
                
                if distances_vec[idx_min] > distance_to_vector_limit*1000:
                    print(f"Nearest vector {distances_vec[idx_min]/1000:.2f} km exceeds limit at time {dtime} for site {site}")
                    continue

                chosen = radial_data.iloc[idx_min]
                velocity_cm_s = chosen['VELO']
                velocity_m_s = velocity_cm_s / 100.0

                heading_raw = chosen['HEAD']
                heading_degrees = heading_raw * 0.1
                heading_rad = np.radians(heading_degrees)
                u = velocity_m_s * np.sin(heading_rad)
                v = velocity_m_s * np.cos(heading_rad)
                lon = chosen['LOND']
                lat = chosen['LATD']
                distance_to_vector = distances_vec[idx_min]

                time_list.append(dtime)
                velocity_list.append(velocity_m_s)
                u_list.append(u)
                v_list.append(v)
                heading_list.append(heading_degrees)
                lon_list.append(lon)
                lat_list.append(lat)
                distance_to_vector_list.append(distance_to_vector)

            except Exception as e:
                print(f"Error reading {radial_file}: {e}")
                continue

        if len(time_list) > 0:
            ds_out = xr.Dataset(
                {
                    'velocity': (['time'], velocity_list),
                    'u': (['time'], u_list),
                    'v': (['time'], v_list),
                    'heading': (['time'], heading_list),
                    'lon': (['time'], lon_list),
                    'lat': (['time'], lat_list),
                    'distance_to_vector': (['time'], distance_to_vector_list)
                },
                coords={'time': time_list}
            )
            point_selected_site_radial_datasets[site] = ds_out
            print(f"Final dataset created for {site} with {len(time_list)} points")
        else:
            print(f"No valid data points for {site}")
            
        full_site_radial_lists[site] = radial_files

    return point_selected_site_radial_datasets, full_site_radial_lists

def calculate_drifter_radial_velocities(drifter_ds, radial_dss, site_loc_df):
    """Calculate drifter radial velocities for comparison with site radials"""
    drifter_radial_vel_ds = {}
    
    for site, radial_ds in radial_dss.items():
        times = radial_ds.time.values
        site_info = site_loc_df.loc[site_loc_df['site'] == site]
        if site_info.empty:
            continue
            
        site_lat = site_info['lat'].values[0]
        site_lon = site_info['lon'].values[0]
        
        d_sel = drifter_ds.sel(time=times)
        d_u = d_sel['u'].values
        d_v = d_sel['v'].values
        d_lat = d_sel['lat'].values
        d_lon = d_sel['lon'].values
        
        # Calculate radial velocity
        drifter_radial_velocities = radial_velocity_to_site(
            site_lat, site_lon, d_lat, d_lon, d_u, d_v
        )
        
        drifter_radial_vel_ds[site] = xr.Dataset(
            data_vars={
                'drifter_radial_velocity': (('time'), drifter_radial_velocities),
                'u': (('time'), d_u),
                'v': (('time'), d_v),
            },
            coords={'time': times},
            attrs={
                'site': site,
                'description': "Drifter radial velocity (m s-1) positive toward radar",
                'site_lat': site_lat,
                'site_lon': site_lon,
            }
        )
    
    return drifter_radial_vel_ds

# ============================================================================
# TOTALS DATA PROCESSING
# ============================================================================

def get_totals_times_dict(drifter_ds, totals_ds, freq_resolution, frequency='long'):
    """
    Returns a dict mapping each drifter timestamp to the nearest-total-point Dataset.
    """
    max_km = freq_resolution[frequency]
    totals_times_dict = {}

    # crop totals to drifter time window once
    t0, t1 = drifter_ds.time.values[[0, -1]]
    tot = totals_ds[['u','v']].sel(time=slice(t0, t1))

    for t in drifter_ds.time.values:
        # 1) select the closest time‐slice
        ds_time = tot.sel(time=t, method='nearest')

        # 2) get drifter lat/lon at that time
        drlat = float(drifter_ds.sel(time=t).lat)
        drlon = float(drifter_ds.sel(time=t).lon)

        # 3) find the nearest grid‐cell in that slice
        ds_pt = ds_time.sel(
            lat = drlat,
            lon = drlon,
            method = 'nearest'
        )

        # 4) compute true haversine distance (km)
        toplat = float(ds_pt.lat)
        toplon = float(ds_pt.lon)
        dkm = haversine(drlat, drlon, toplat, toplon) / 1000.0

        # 5) only keep it if within your threshold
        if dkm <= max_km:
            ds_pt = ds_pt.assign_coords(distance = dkm, time = t)
            totals_times_dict[t] = ds_pt

    return totals_times_dict

def get_totals_times_dict_interp(drifter_ds, totals_ds, freq_resolution, 
                                frequency='long', max_time_offset='15min'):
    """
    Map each drifter timestamp to a bilinearly-interpolated HF-radar vector.
    """
    max_km = freq_resolution[frequency]
    out = {}

    # restrict radar time range once
    t0, t1 = drifter_ds.time.values[[0, -1]]
    tot = totals_ds[['u', 'v']].sel(time=slice(t0, t1))

    for t in drifter_ds.time.values:
        # pick the radar time-slice closest to this fix
        slice_time = tot.sel(time=t, method='nearest')

        # option: enforce a maximum time mismatch
        if max_time_offset is not None:
            dt = abs((np.datetime64(t) - slice_time.time.values).astype('timedelta64[s]'))
            if dt > pd.Timedelta(max_time_offset):
                continue

        # drifter position
        drlat = float(drifter_ds.sel(time=t).lat)
        drlon = float(drifter_ds.sel(time=t).lon)

        # distance to centre of nearest grid cell (for QC)
        nearest_pt = slice_time.sel(lat=drlat, lon=drlon, method='nearest')
        d_km = haversine(drlat, drlon,
                         float(nearest_pt.lat), float(nearest_pt.lon)) / 1000.0
        if d_km > max_km:
            continue

        # bilinear interpolation to the exact drifter coord
        interp_pt = slice_time.interp(lat=drlat, lon=drlon, method='linear')

        # make sure lon/lat coords are exactly the drifter pos
        interp_pt = interp_pt.assign_coords(
            lon=drlon, lat=drlat, time=t, distance=d_km
        )

        out[t] = interp_pt

    return out

# ============================================================================
# RADIAL PLOTTING FUNCTIONS
# ============================================================================

def plot_radial_drifter_track(radial_dss, drifter_ds, site_name, extent=None, save_dir=None):
    """Plot drifter track with radial velocities for a specific site"""
    if site_name not in radial_dss:
        print(f"Site {site_name} not found in radial data")
        return
        
    radial_ds = radial_dss[site_name]
    
    if extent is None:
        extent = [-75.5, -69, 37, 42]
        
    fig, ax = cplt.create(extent=extent, proj=ccrs.Mercator())
    ax.set_title(f"Radial Velocities and Drifter Path for {site_name}")
    
    # Plot radial velocities
    sc = plt.scatter(radial_ds.lon, radial_ds.lat, c=radial_ds.velocity, 
                    cmap='viridis', s=10, transform=ccrs.PlateCarree(), 
                    label='Radial Velocities')
    
    # Plot drifter path
    plt.plot(drifter_ds.lon, drifter_ds.lat, 'r-', transform=ccrs.PlateCarree(), 
            label='Drifter Path', alpha=0.4)
    
    plt.colorbar(sc, ax=ax, label='Velocity (m/s)')
    plt.legend(loc='lower right')
    
    if save_dir is not None:
        plt.savefig(os.path.join(save_dir, f'{site_name}_track_map.png'), dpi=300, bbox_inches='tight')
    plt.close()

def plot_radial_velocity_timeseries(drifter_radial_vel_ds, radial_dss, site_name, 
                                   save_dir=None, pad_time=True):
    """Plot time series comparison of drifter vs site radial velocities"""
    if site_name not in drifter_radial_vel_ds:
        return
        
    drifter_ds = drifter_radial_vel_ds[site_name]
    radial_ds = radial_dss[site_name]
    
    if pad_time and len(radial_ds.time) > 1:
        # Pad to hourly intervals
        t0 = radial_ds.time.values[0]
        t1 = radial_ds.time.values[-1]
        time_range = pd.date_range(t0, t1, freq='1H')
        drifter_ds = drifter_ds.reindex(time=time_range, method='nearest', 
                                       tolerance=pd.Timedelta('1H'))
        radial_ds = radial_ds.reindex(time=time_range, method='nearest', 
                                     tolerance=pd.Timedelta('1H'))
    
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1)
    
    drifter_ds['drifter_radial_velocity'].plot(ax=ax, label="Drifter Calculated Radial Velocity")
    radial_ds['velocity'].plot(ax=ax, label=f"Gathered Site Radial Velocity")
    
    plt.ylim(-0.5, 0.5)
    plt.yticks(np.arange(-0.5, 0.6, 0.1))
    plt.legend()
    plt.title(f"Drifter Radial Velocity vs. Site Radial Velocity\n{site_name}")
    plt.xlabel("Time")
    plt.ylabel("Velocity (m/s)")
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=7))
    plt.gca().xaxis.set_minor_locator(mdates.DayLocator(interval=1))
    plt.xticks(rotation=45)
    plt.grid(which='both', lw=0.2)
    
    if save_dir is not None:
        plt.savefig(os.path.join(save_dir, f'{site_name}_velocity_timeseries.png'), 
                   dpi=300, bbox_inches='tight')
    plt.close()

def plot_radial_velocity_difference(drifter_radial_vel_ds, radial_dss, site_name, 
                                   save_dir=None, pad_time=True):
    """Plot time series of velocity differences"""
    if site_name not in drifter_radial_vel_ds:
        return
        
    drifter_ds = drifter_radial_vel_ds[site_name]
    radial_ds = radial_dss[site_name]
    
    if pad_time and len(radial_ds.time) > 1:
        # Pad to hourly intervals
        t0 = radial_ds.time.values[0]
        t1 = radial_ds.time.values[-1]
        time_range = pd.date_range(t0, t1, freq='1H')
        drifter_ds = drifter_ds.reindex(time=time_range, method='nearest', 
                                       tolerance=pd.Timedelta('1H'))
        radial_ds = radial_ds.reindex(time=time_range, method='nearest', 
                                     tolerance=pd.Timedelta('1H'))
    
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1)
    
    diff = drifter_ds['drifter_radial_velocity'] - radial_ds['velocity']
    diff.plot(ax=ax, label="Difference in Radial Velocity")
    
    plt.ylim(-0.5, 0.5)
    plt.yticks(np.arange(-0.5, 0.6, 0.1))
    plt.legend()
    plt.title(f"Difference in Drifter Radial Velocity vs. Site Radial Velocity\n{site_name}")
    plt.xlabel("Time")
    plt.ylabel("Velocity (m/s)")
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=7))
    plt.gca().xaxis.set_minor_locator(mdates.DayLocator(interval=1))
    plt.xticks(rotation=45)
    plt.grid(which='both', lw=0.2)
    
    if save_dir is not None:
        plt.savefig(os.path.join(save_dir, f'{site_name}_difference_timeseries.png'), 
                   dpi=300, bbox_inches='tight')
    plt.close()

def plot_radial_velocity_with_distance(drifter_radial_vel_ds, radial_dss, site_name,
                                      site_frequency_dict, freq_resolution, save_dir=None):
    """Plot radial velocity timeseries with distance to vector overlay"""
    if site_name not in drifter_radial_vel_ds:
        return
        
    ds = drifter_radial_vel_ds[site_name]
    radial_ds = radial_dss[site_name]
    
    if len(radial_ds.time) > 1:
        # pad to hourly intervals
        t0 = radial_ds.time.values[0]
        t1 = radial_ds.time.values[-1]
        time_range = pd.date_range(t0, t1, freq='1H')
        ds = ds.reindex(time=time_range, method='nearest',
                       tolerance=pd.Timedelta('1H'))
        radial_ds = radial_ds.reindex(time=time_range, method='nearest',
                                     tolerance=pd.Timedelta('1H'))
    
    # Get time, distance, velocities
    time = ds['drifter_radial_velocity'].time
    dist = radial_ds['distance_to_vector'] / 1000  # in km
    dr_v = ds['drifter_radial_velocity'].values * 100  # m/s → cm/s
    hr_v = radial_ds['velocity'].values * 100

    fig, ax_vel = plt.subplots(figsize=(12, 4))

    # Distance envelope on the right‐axis
    ax_dist = ax_vel.twinx()
    ax_dist.fill_between(time, dist, -dist,
                         color='red', alpha=0.3,
                         label='Distance to Closest Radial Measurement')
    ax_dist.set_ylabel('Distance (km)', color='red')
    ylims = get_freq_resolution_for_site(site_name, site_frequency_dict, freq_resolution)
    if ylims:
        ax_dist.set_ylim(-ylims, ylims)

    # Radial velocities on the left‐axis
    ax_vel.plot(time, dr_v, color='green', label='Drifter Radial Velocity')
    ax_vel.plot(time, hr_v, color='blue', label=f'{site_name} Radial Velocity')
    ax_vel.set_ylabel('Radial Velocity (cm/s)')
    ax_vel.set_ylim(-100, 100)

    # Formatting
    ax_vel.xaxis.set_major_locator(mdates.DayLocator(interval=7))
    ax_vel.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax_vel.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
    ax_vel.xaxis.set_tick_params(rotation=45)
    ax_vel.set_xlabel('Time')
    ax_vel.set_title(f'Drifter vs {site_name} Radial Velocity\nwith Distance to Nearest Measurement')
    ax_vel.grid(True, which='both', lw=0.2)

    # Legend
    ax_vel.legend(loc='upper left')
    
    plt.tight_layout()
    
    if save_dir is not None:
        plt.savefig(os.path.join(save_dir, f'{site_name}_velocity_with_distance.png'), 
                   dpi=300, bbox_inches='tight')
    plt.close()

def plot_radial_correlation_scatter(drifter_radial_vel_ds, radial_dss, site_name, save_dir=None):
    """Create scatter plot comparing drifter and HFR radial velocities"""
    if site_name not in drifter_radial_vel_ds:
        return
        
    drifter_site_ds = drifter_radial_vel_ds[site_name]
    hfr_site_ds = radial_dss[site_name]
    
    # Filter NaN values
    valid_idx = ~np.isnan(drifter_site_ds.drifter_radial_velocity) & ~np.isnan(hfr_site_ds.velocity)
    drifter_valid = drifter_site_ds.drifter_radial_velocity[valid_idx]
    hfr_valid = hfr_site_ds.velocity[valid_idx]
    time_valid = pd.to_datetime(drifter_site_ds.time.values)[valid_idx]
    
    if len(drifter_valid) < 2:
        print(f"Not enough valid points for scatter plot for {site_name}")
        return
    
    # Create scatter plot
    fig, ax = plt.subplots(figsize=(10, 10))
    sc = ax.scatter(drifter_valid, hfr_valid, c=mdates.date2num(time_valid), 
                   cmap='jet', edgecolor='k', s=50)
    
    # Add line of best fit
    z = np.polyfit(drifter_valid, hfr_valid, 1)
    p = np.poly1d(z)
    x_vals = np.linspace(drifter_valid.min().values.item(), 
                        drifter_valid.max().values.item(), 100)
    slope = z[0]
    ax.plot(x_vals, p(x_vals), "r--", label=f'Best-Fit Line - Slope: {slope:.2f}')
    
    # Add 1:1 correlation line
    axlim = [min(drifter_valid.min(), hfr_valid.min()) - 0.1, 
             max(drifter_valid.max(), hfr_valid.max()) + 0.1]
    ax.plot(axlim, axlim, 'k--', label='1:1 Line')
    ax.set_xlim(axlim)
    ax.set_ylim(axlim)
    
    # Add colorbar
    cb = plt.colorbar(sc, label='Time')
    cb.ax.yaxis.set_major_formatter(mdates.DateFormatter('%m/%d %H:%M'))
    cb.ax.yaxis.set_major_locator(mdates.AutoDateLocator())
    # Make colorbar text larger and bold
    cb.set_label('Time', fontsize=16, fontweight='bold')
    cb.ax.tick_params(labelsize=14)
    for label in cb.ax.get_yticklabels():
        label.set_fontweight('bold')
    
    # Customize plot
    ax.grid(True)
    ax.set_title(f'Drifter vs HFR Radial Velocities for {site_name}', fontsize=18, fontweight='bold')
    ax.set_xlabel('Drifter Radial Velocity (m/s)', fontsize=16, fontweight='bold')
    ax.set_ylabel('HFR Radial Velocity (m/s)', fontsize=16, fontweight='bold')
    # Make axis tick labels larger and bold
    ax.tick_params(axis='both', which='major', labelsize=14)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')
    ax.legend(loc='upper left', fontsize=14, prop={'weight': 'bold'})
    
    plt.tight_layout()
    
    if save_dir is not None:
        plt.savefig(os.path.join(save_dir, f'{site_name}_correlation_scatter.png'), 
                   dpi=300, bbox_inches='tight')
    plt.close()

# ============================================================================
# TOTALS PLOTTING FUNCTIONS
# ============================================================================

def plot_totals_drifter_track(drifter_ds, save_dir=None):
    """Plot the drifter track with points colored by time"""
    drifter_lon_min = drifter_ds.lon.min().values
    drifter_lon_max = drifter_ds.lon.max().values
    drifter_lat_min = drifter_ds.lat.min().values
    drifter_lat_max = drifter_ds.lat.max().values
    dx = 0.5
    dy = 0.5
    extent = [drifter_lon_min-dx, drifter_lon_max+dx, drifter_lat_min-dy, drifter_lat_max+dy]
    
    fig, ax = cplt.create(extent=extent, proj=ccrs.Mercator(), figsize=(10, 8))
    
    times = mdates.date2num(pd.to_datetime(drifter_ds.time.values)) 
    
    # Plot drifter track
    sc = ax.scatter(drifter_ds.lon, drifter_ds.lat, c=times, 
                   cmap='jet', s=20, transform=ccrs.PlateCarree())
    ax.plot(drifter_ds.lon, drifter_ds.lat, 'k-', linewidth=1, transform=ccrs.PlateCarree())
    
    # Add colorbar
    cb = plt.colorbar(sc, label='Time')
    cb.ax.yaxis.set_major_formatter(mdates.DateFormatter('%m/%d %H:%M'))
    cb.ax.yaxis.set_major_locator(mdates.AutoDateLocator())
    
    ax.grid(True)
    plt.title(f'Drifter Track for Drifter: {int(drifter_ds.drifter_id)}', fontsize=16, fontweight='bold')
    
    if save_dir is not None:
        plt.savefig(os.path.join(save_dir, 'track_map.png'), dpi=300, bbox_inches='tight')
    plt.close()

def plot_current_roses(drifter_ds, totals_times_ds, save_dir=None):
    """Plot current roses comparing drifter and HFR data"""
    fig = plt.figure(figsize=(15, 5))
    
    # Convert drifter velocities to direction and speed
    drifter_dir, drifter_spd = uv2spdir(drifter_ds.u.values, drifter_ds.v.values)
    
    # Get HFR velocities
    hfr_u = totals_times_ds.u.values
    hfr_v = totals_times_ds.v.values
    hfr_dir, hfr_spd = uv2spdir(hfr_u, hfr_v)
    
    # Calculate difference
    err_u = hfr_u - drifter_ds.u.values
    err_v = hfr_v - drifter_ds.v.values
    err_dir, err_spd = uv2spdir(err_u, err_v)
    
    # Plot drifter rose
    ax1 = fig.add_subplot(131, projection='windrose')
    ax1.bar(drifter_dir, drifter_spd, bins=np.arange(0, 1.1, 0.1), 
            nsector=32, normed=True)
    ax1.set_title('Drifter')
    
    # Plot HFR rose
    ax2 = fig.add_subplot(132, projection='windrose')
    ax2.bar(hfr_dir, hfr_spd, bins=np.arange(0, 1.1, 0.1), 
            nsector=32, normed=True)
    ax2.set_title('HFR Totals')
    
    # Plot difference rose
    ax3 = fig.add_subplot(133, projection='windrose')
    ax3.bar(err_dir, err_spd, bins=np.arange(0, 0.8, 0.1),
            nsector=32, normed=True)
    ax3.set_title('HFR - Drifter')
    
    for ax in (ax1, ax2, ax3):
        ax.set_ylim(0, 6)
        ticks = np.arange(0, 8, 2)
        ax.set_yticks(ticks)
        ax.set_yticklabels([f'{int(t)}%' for t in ticks])
    
    ax3.legend(loc='lower left', bbox_to_anchor=(1, 0.61))
    
    plt.tight_layout()
    
    if save_dir is not None:
        plt.savefig(os.path.join(save_dir, 'current_roses.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    return fig

def plot_totals_uv_timeseries(drifter_ds, totals_times_ds, save_dir=None):
    """Plot u and v component time series comparison"""
    fig = plt.figure(figsize=(15, 10))
    
    hfr_u = totals_times_ds.u.values
    hfr_v = totals_times_ds.v.values
    hfr_times = totals_times_ds.time.values
    
    # U component
    ax1 = fig.add_subplot(311)
    ax1.plot(drifter_ds.time, drifter_ds.u, label='Drifter u')
    ax1.plot(hfr_times, hfr_u, label='HFR u')
    ax1.xaxis.set_minor_locator(mdates.DayLocator())
    ax1.set_title('u Component')
    ax1.legend(loc='upper right')
    ax1.grid(True)
    ax1.set_ylabel('Velocity (m/s)')
    
    plt.subplots_adjust(hspace=0.3)
    
    # V component
    ax2 = fig.add_subplot(312)
    ax2.plot(drifter_ds.time, drifter_ds.v, label='Drifter v')
    ax2.plot(hfr_times, hfr_v, label='HFR v')
    ax2.xaxis.set_minor_locator(mdates.DayLocator())
    ax2.set_title('v Component')
    ax2.legend(loc='upper right')
    ax2.grid(True)
    ax2.set_ylabel('Velocity (m/s)')
    
    # Difference
    ax3 = fig.add_subplot(313)
    u_diff = hfr_u - drifter_ds.sel(time=hfr_times).u
    v_diff = hfr_v - drifter_ds.sel(time=hfr_times).v
    ax3.plot(hfr_times, u_diff, label='HFR - Drifter u')
    ax3.plot(hfr_times, v_diff, label='HFR - Drifter v')
    ax3.xaxis.set_minor_locator(mdates.DayLocator())
    ax3.set_title('Difference in u and v Components')
    ax3.legend(loc='upper right')
    ax3.grid(True)
    ax3.set_ylabel('Velocity Difference (m/s)')
    ax3.set_xlabel('Time')
    
    plt.tight_layout()
    
    if save_dir is not None:
        plt.savefig(os.path.join(save_dir, 'uv_timeseries.png'), dpi=300, bbox_inches='tight')
    plt.close()

def create_animation_frames(drifter_ds, totals_times_ds, totals_ds, save_dir=None, 
                          days=3, start_time=None, end_time=None):
    """Create frames for animation of drifter moving through totals field"""
    if start_time is None or end_time is None:
        start_time, end_time = get_time_subset(drifter_ds, days=days)
    else:
        start_time = pd.to_datetime(start_time)
        end_time = pd.to_datetime(end_time)
    
    drifter_ds_sel = drifter_ds.sel(time=slice(start_time, end_time))
    totals_ds_sel = totals_times_ds.sel(time=slice(start_time, end_time))
    
    extent = [drifter_ds_sel.lon.min().values - 0.05, 
              drifter_ds_sel.lon.max().values + 0.05, 
              drifter_ds_sel.lat.min().values - 0.05, 
              drifter_ds_sel.lat.max().values + 0.05]
    
    frames_dir = os.path.join(save_dir, 'animation_frames') if save_dir else 'animation_frames'
    os.makedirs(frames_dir, exist_ok=True)
    
    n_time_steps = drifter_ds_sel.time.size
    
    for i in range(1, n_time_steps + 1):
        drifter_lim = drifter_ds_sel.isel(time=slice(0, i))
        totals_lim = totals_ds_sel.isel(time=slice(0, i))
        
        frame_name = f'frame_{i:03d}.png'
        
        # Create the plot
        fig, ax = cplt.create(extent=extent, proj=ccrs.Mercator(), figsize=(10, 8))
        
        # Plot background totals grid
        latest_time = pd.to_datetime(drifter_lim.time.values[-1])
        totals_last = totals_ds.sel(time=latest_time, method='nearest')
        totals_extent = [extent[0]-1, extent[1]+1, extent[2]-1, extent[3]+1]
        totals_last = totals_last.sel(lon=slice(totals_extent[0], totals_extent[1]),
                                     lat=slice(totals_extent[2], totals_extent[3]))
        
        # Plot drifter track
        times_num = mdates.date2num(pd.to_datetime(drifter_lim.time.values))
        sc = ax.scatter(drifter_lim.lon, drifter_lim.lat, c=times_num, cmap='turbo', 
                       s=30, transform=ccrs.PlateCarree())
        ax.plot(drifter_lim.lon, drifter_lim.lat, 'grey', linewidth=0.5, 
               transform=ccrs.PlateCarree())
        
        # Plot velocity vectors
        ax.quiver(drifter_lim.lon, drifter_lim.lat, drifter_lim.u, drifter_lim.v,
                 color='blue', scale=10, width=0.004, label='Drifter Velocity', 
                 transform=ccrs.PlateCarree())
        ax.quiver(drifter_lim.lon, drifter_lim.lat, totals_lim.u, totals_lim.v,
                 color='red', scale=10, width=0.004, label='HFR Velocity', 
                 transform=ccrs.PlateCarree())
        
        ax.set_title(f'Drifter Track Velocity Comparison\nFrame {i}/{n_time_steps}')
        ax.legend()
        
        plt.savefig(os.path.join(frames_dir, frame_name), dpi=150, bbox_inches='tight')
        plt.close()
    
    print(f"Created {n_time_steps} animation frames in {frames_dir}")
    return frames_dir

def plot_totals_correlation_scatter(drifter_ds, totals_times_ds, save_dir=None):
    """Create scatter plots for U and V velocity correlations"""
    dr_u = drifter_ds.u.values
    dr_v = drifter_ds.v.values
    dr_time = drifter_ds.time.values
    tot_u = totals_times_ds.u.values
    tot_v = totals_times_ds.v.values
    
    # Filter NaN values
    valid_idx = ~np.isnan(dr_u) & ~np.isnan(tot_u)
    dr_u_valid = dr_u[valid_idx]
    dr_v_valid = dr_v[valid_idx]
    tot_u_valid = tot_u[valid_idx]
    tot_v_valid = tot_v[valid_idx]
    dr_time_valid = pd.to_datetime(dr_time)[valid_idx]
    
    # Create scatter plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    
    # Plot U
    sc1 = ax1.scatter(dr_u_valid, tot_u_valid, c=mdates.date2num(dr_time_valid), 
                     cmap='turbo', edgecolor='k', s=50)
    
    # Add line of best fit
    z1 = np.polyfit(dr_u_valid, tot_u_valid, 1)
    p1 = np.poly1d(z1)
    x_vals1 = np.linspace(dr_u_valid.min(), dr_u_valid.max(), 100)
    slope1 = z1[0]
    ax1.plot(x_vals1, p1(x_vals1), "r--", label=f'Best-Fit Line - Slope: {slope1:.2f}')
    
    # Add 1:1 line
    axlim1 = [min(dr_u_valid.min(), tot_u_valid.min()) - 0.1, 
              max(dr_u_valid.max(), tot_u_valid.max()) + 0.1]
    ax1.plot(axlim1, axlim1, 'k--', label='1:1 Line')
    ax1.set_xlim(axlim1)
    ax1.set_ylim(axlim1)
    
    # Customize
    ax1.grid(True)
    ax1.set_title('Drifter vs HFR U Velocities', fontweight='bold', fontsize=16)
    ax1.set_xlabel('Drifter U Velocity (m/s)')
    ax1.set_ylabel('Totals U Velocity (m/s)')
    ax1.legend(loc='upper left')
    
    # Plot V
    sc2 = ax2.scatter(dr_v_valid, tot_v_valid, c=mdates.date2num(dr_time_valid), 
                     cmap='turbo', edgecolor='k', s=50)
    
    # Add line of best fit
    z2 = np.polyfit(dr_v_valid, tot_v_valid, 1)
    p2 = np.poly1d(z2)
    x_vals2 = np.linspace(dr_v_valid.min(), dr_v_valid.max(), 100)
    slope2 = z2[0]
    ax2.plot(x_vals2, p2(x_vals2), "r--", label=f'Best-Fit Line - Slope: {slope2:.2f}')
    
    # Add 1:1 line
    axlim2 = [min(dr_v_valid.min(), tot_v_valid.min()) - 0.1, 
              max(dr_v_valid.max(), tot_v_valid.max()) + 0.1]
    ax2.plot(axlim2, axlim2, 'k--', label='1:1 Line')
    ax2.set_xlim(axlim2)
    ax2.set_ylim(axlim2)
    
    # Customize
    ax2.grid(True)
    ax2.set_title('Drifter vs HFR V Velocities', fontweight='bold', fontsize=16)
    ax2.set_xlabel('Drifter V Velocity (m/s)')
    ax2.set_ylabel('Totals V Velocity (m/s)')
    ax2.legend(loc='upper left')
    
    # Add colorbars
    cb1 = plt.colorbar(sc1, ax=ax1, label='Time')
    cb1.ax.yaxis.set_major_formatter(mdates.DateFormatter('%m/%d %H:%M'))
    cb2 = plt.colorbar(sc2, ax=ax2, label='Time')
    cb2.ax.yaxis.set_major_formatter(mdates.DateFormatter('%m/%d %H:%M'))
    
    plt.tight_layout()
    
    if save_dir is not None:
        plt.savefig(os.path.join(save_dir, 'correlation_scatter.png'), dpi=300, bbox_inches='tight')
    plt.close()

def plot_improved_uv_timeseries(drifter_ds, totals_times_ds, save_dir=None):
    """Create improved U/V time series plots with better formatting"""
    dr_u = drifter_ds.u.values
    dr_v = drifter_ds.v.values
    dr_time = drifter_ds.time.values
    tot_u = totals_times_ds.u.values
    tot_v = totals_times_ds.v.values
    tot_time = totals_times_ds.time.values
    
    # Create time series plot
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10))
    
    # Plot U
    ax1.plot(dr_time, dr_u, label='Drifter U Velocity', color='blue')
    ax1.plot(tot_time, tot_u, label='Totals U Velocity', color='red')
    ax1.set_ylim(-0.6, 0.6)
    ax1.set_title('Drifter vs HFR U Velocities', fontweight='bold', fontsize=16)
    ax1.set_ylabel('Velocity (m/s)')
    ax1.legend(loc='upper right')
    ax1.grid(True, which='both')
    ax1.xaxis.set_major_locator(mdates.DayLocator(interval=7))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax1.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
    ax1.xaxis.set_tick_params(rotation=45, labelbottom=False)
    
    # Plot V
    ax2.plot(dr_time, dr_v, label='Drifter V Velocity', color='blue')
    ax2.plot(tot_time, tot_v, label='Totals V Velocity', color='red')
    ax2.set_ylim(-0.6, 0.6)
    ax2.set_title('Drifter vs HFR V Velocities', fontweight='bold', fontsize=16)
    ax2.set_ylabel('Velocity (m/s)')
    ax2.legend(loc='upper right')
    ax2.grid(True, which='both')
    ax2.xaxis.set_major_locator(mdates.DayLocator(interval=7))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax2.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
    ax2.xaxis.set_tick_params(rotation=45, labelbottom=False)
    
    # Plot difference
    diff_u = dr_u - tot_u
    diff_v = dr_v - tot_v
    ax3.plot(dr_time, diff_u, label='Drifter U - Totals U', color='green')
    ax3.plot(dr_time, diff_v, label='Drifter V - Totals V', color='orange')
    
    ax3.set_ylim(-0.6, 0.6)
    ax3.set_title('Drifter vs HFR U and V Velocity Differences', fontweight='bold', fontsize=16)
    ax3.set_xlabel('Time')
    ax3.set_ylabel('Velocity Difference (m/s)')
    ax3.legend(loc='upper right')
    ax3.grid(True, which='both')
    ax3.xaxis.set_major_locator(mdates.DayLocator(interval=7))
    ax3.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax3.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
    ax3.xaxis.set_tick_params(rotation=45, labelbottom=True)

    plt.tight_layout()
    
    if save_dir is not None:
        plt.savefig(os.path.join(save_dir, 'improved_uv_timeseries.png'), dpi=300, bbox_inches='tight')
    plt.close()

# ============================================================================
# MASTER PROCESSING FUNCTIONS
# ============================================================================

def process_radial_comparison(drifter_file, qc_params, radial_params, figure_dict, 
                            output_dir, drifter_id):
    """
    Master function to process radial comparisons
    
    Parameters:
    -----------
    drifter_file : str
        Path to drifter file
    qc_params : dict
        QC parameters for drifter data
    radial_params : dict
        Parameters for radial processing including:
        - 'site_file': path to site locations Excel file
        - 'radial_src_dir': base directory for radial files
        - 'site_frequency_dict': dict mapping frequencies to sites
        - 'freq_resolution': dict of frequency resolutions
        - 'sites_to_process': 'all' or list of specific sites
        - 'proximity': distance threshold in km (optional)
    figure_dict : dict
        Dictionary of which figures to create
    output_dir : str
        Base output directory
    drifter_id : str
        Drifter ID for output naming
    
    Returns:
    --------
    dict : Results and paths to saved figures
    """
    # Create output directory
    timestamp = datetime.now().strftime('%Y_%m_%dT%H')
    comparison_dir = os.path.join(output_dir, f'{drifter_id}_{timestamp}_hfr_drifter_comparison')
    radial_dir = os.path.join(comparison_dir, 'Radials')
    os.makedirs(radial_dir, exist_ok=True)
    
    # Load drifter data with QC
    print("Loading drifter data with QC...")
    drifter_ds = load_drifter(drifter_file, qc_params)
    
    # Load site locations
    site_loc_df = pd.read_excel(radial_params['site_file'])
    
    # Get nearest radial data
    print("Getting nearest radial data...")
    radial_dss, radial_paths = get_nearest_radial(
        drifter_ds, 
        radial_params['site_frequency_dict'],
        radial_params['freq_resolution'],
        site_loc_df,
        radial_params['radial_src_dir'],
        radial_params.get('sites_to_process', 'all'),
        radial_params.get('proximity', None)
    )
    
    if not radial_dss:
        print("No radial data found for any sites")
        return {'status': 'no_data', 'output_dir': comparison_dir}
    
    total_times = drifter_ds.time.size
    for site, ds in radial_dss.items():
        found_times = ds.time.size
        print("\n")
        print("="*30)
        print(f"Site {site}: Found radial data for {found_times} out of {total_times} drifter time points")
        print("="*30)
        print("\n")
    # Calculate drifter radial velocities
    print("Calculating drifter radial velocities...")
    drifter_radial_vel_ds = calculate_drifter_radial_velocities(
        drifter_ds, radial_dss, site_loc_df
    )
    
    # Generate figures for each site
    results = {'sites_processed': [], 'figures_created': {}}
    
    for site in radial_dss.keys():
        print(f"\nGenerating figures for site {site}...")
        site_dir = os.path.join(radial_dir, site)
        os.makedirs(site_dir, exist_ok=True)
        
        site_figures = []
        
        if figure_dict.get('track_map', True):
            plot_radial_drifter_track(radial_dss, drifter_ds, site, save_dir=site_dir)
            site_figures.append('track_map')
            
        if figure_dict.get('velocity_timeseries', True):
            plot_radial_velocity_timeseries(drifter_radial_vel_ds, radial_dss, site, 
                                          save_dir=site_dir)
            site_figures.append('velocity_timeseries')
            
        if figure_dict.get('difference_timeseries', True):
            plot_radial_velocity_difference(drifter_radial_vel_ds, radial_dss, site, 
                                          save_dir=site_dir)
            site_figures.append('difference_timeseries')
            
        if figure_dict.get('velocity_with_distance', True):
            plot_radial_velocity_with_distance(drifter_radial_vel_ds, radial_dss, site,
                                             radial_params['site_frequency_dict'],
                                             radial_params['freq_resolution'],
                                             save_dir=site_dir)
            site_figures.append('velocity_with_distance')
            
        if figure_dict.get('correlation_scatter', True):
            plot_radial_correlation_scatter(drifter_radial_vel_ds, radial_dss, site, 
                                          save_dir=site_dir)
            site_figures.append('correlation_scatter')
        
        results['sites_processed'].append(site)
        results['figures_created'][site] = site_figures
    
    results['status'] = 'success'
    results['output_dir'] = comparison_dir
    results['radial_dss'] = radial_dss
    results['drifter_radial_vel_ds'] = drifter_radial_vel_ds
    
    return results

def process_totals_comparison(drifter_file, qc_params, totals_params, figure_dict, 
                            output_dir, drifter_id):
    """
    Master function to process totals comparisons
    
    Parameters:
    -----------
    drifter_file : str
        Path to drifter file
    qc_params : dict
        QC parameters for drifter data
    totals_params : dict
        Parameters for totals processing including:
        - 'totals_url': URL to totals dataset
        - 'frequency': 'long', 'medium', or 'short'
        - 'freq_resolution': dict of frequency resolutions
        - 'use_interpolation': bool, whether to use interpolation
    figure_dict : dict
        Dictionary of which figures to create
    output_dir : str
        Base output directory
    drifter_id : str
        Drifter ID for output naming
    
    Returns:
    --------
    dict : Results and paths to saved figures
    """
    # Create output directory
    timestamp = datetime.now().strftime('%Y_%m_%dT%H')
    comparison_dir = os.path.join(output_dir, f'{drifter_id}_{timestamp}_hfr_drifter_comparison')
    totals_dir = os.path.join(comparison_dir, 'Totals')
    os.makedirs(totals_dir, exist_ok=True)
    
    # Load drifter data with QC
    print("Loading drifter data with QC...")
    drifter_ds = load_drifter(drifter_file, qc_params)
    
    # Load totals data
    print("Loading totals data...")
    totals_ds = xr.open_dataset(totals_params['totals_url'], chunks={'time': 'auto'})
    totals_ds = totals_ds[['u', 'v']].squeeze()
    
    # Subset to drifter time range
    drifter_start = drifter_ds.time.values[0]
    drifter_end = drifter_ds.time.values[-1]
    totals_ds = totals_ds.sel(time=slice(drifter_start, drifter_end))
    
    # Get totals at drifter positions
    print("Extracting totals at drifter positions...")
    if 'use_interpolation' in totals_params and totals_params['use_interpolation']:
        print("WARNING: Using interpolation for totals extraction, this will increase run time.")
    if totals_params.get('use_interpolation', False):
        totals_times_dict = get_totals_times_dict_interp(
            drifter_ds, totals_ds, totals_params['freq_resolution'],
            totals_params['frequency']
        )
    else:
        totals_times_dict = get_totals_times_dict(
            drifter_ds, totals_ds, totals_params['freq_resolution'],
            totals_params['frequency']
        )
    
    if not totals_times_dict:
        print("No totals data found at drifter positions")
        return {'status': 'no_data', 'output_dir': comparison_dir}
    
    totals_times_ds = xr.concat(list(totals_times_dict.values()), dim='time')
    totals_times_ds.load()
    
    # Generate figures
    results = {'figures_created': []}
    
    if figure_dict.get('track_map', True):
        plot_totals_drifter_track(drifter_ds, save_dir=totals_dir)
        results['figures_created'].append('track_map')
        
    if figure_dict.get('current_roses', True):
        plot_current_roses(drifter_ds, totals_times_ds, save_dir=totals_dir)
        results['figures_created'].append('current_roses')
        
    if figure_dict.get('uv_timeseries', True):
        plot_totals_uv_timeseries(drifter_ds, totals_times_ds, save_dir=totals_dir)
        results['figures_created'].append('uv_timeseries')
        
    if figure_dict.get('animation_frames', True):
        frames_dir = create_animation_frames(drifter_ds, totals_times_ds, totals_ds, 
                                           save_dir=totals_dir)
        results['figures_created'].append('animation_frames')
        results['frames_dir'] = frames_dir
        
    if figure_dict.get('correlation_scatter', True):
        plot_totals_correlation_scatter(drifter_ds, totals_times_ds, save_dir=totals_dir)
        results['figures_created'].append('correlation_scatter')
        
    if figure_dict.get('improved_uv_timeseries', True):
        plot_improved_uv_timeseries(drifter_ds, totals_times_ds, save_dir=totals_dir)
        results['figures_created'].append('improved_uv_timeseries')
    
    results['status'] = 'success'
    results['output_dir'] = comparison_dir
    results['totals_times_ds'] = totals_times_ds
    
    return results

