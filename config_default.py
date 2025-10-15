"""
config-default.py
"""
# ============================================================================
# FILE PATHS - UPDATE THESE
# ============================================================================

# Drifter data
DRIFTER_FILE = '/path/to/your/drifter_file.nc'
DRIFTER_ID = 'your_drifter_id'

# Output directory
OUTPUT_DIR = '/path/to/your/output'

# HFR site locations Excel file
SITE_FILE = '/path/to/HFR_sites.xlsx' # Excel file with HFR site metadata downloaded from repository or custom made

# Base directory for radial files
RADIAL_SRC_DIR = '/path/to/radials' # Directory containing subdirectories for each HFR site with radial files

# Totals data URLs
TOTALS_URLS = {
    'long': 'https://link-to-your-thredds-file-link-for-06km-data.ncd',
    'medium': 'https://links-to-your-thredds-file-link-for-02km-data.ncd',
    'short': 'https://link-to-your-thredds-file-link-for-01km-data.ncd'
}

# ============================================================================
# COMPARISON SELECTION
# ============================================================================

# Choose which comparisons to run
RUN_RADIAL_COMPARISON = True  # Set to False to skip radial comparison
RUN_TOTALS_COMPARISON = True  # Set to False to skip totals comparison

# ============================================================================
# QUALITY CONTROL PARAMETERS
# ============================================================================

QC_PARAMS = {
    # Rate of Change Test
    'apply_roc': True,      # Apply rate of change test
    'roc_threshold': 0.15,  # Maximum rate of change in velocity (m/s)
    
    # Standard Deviation Test
    'apply_std': True,      # Apply standard deviation test
    'std_levels': 3,        # Number of standard deviations from mean
    
    # Distance/Spike Test
    'apply_dist': True,     # Apply distance threshold test
    'dist_threshold': 10000 # Maximum distance between points (meters)
}

# ============================================================================
# RADIAL COMPARISON PARAMETERS
# ============================================================================

# Sites organized by frequency
SITE_FREQUENCY_DICT = {
    'long': ['NANT', 'BLCK', 'AMAG', 'MRCH', 'HEMP', 'HOOK', 'LOVE', 'BRIG', 'WILD'],
    'medium': ['SEAB', 'BRAD', 'SPRK', 'HLGT', 'BRMR', 'RATH', 'WOOD', 'CMPT'],
    'short': ['SILD', 'OLDB', 'PORT', 'HLPN', 'LEWE', 'CAPE']
}

# Maximum distance to radial vector for each frequency (km)
FREQ_RESOLUTION = {
    'long': 5.825,
    'medium': 3.020,
    'short': 1.001
}

RADIAL_PARAMS = {
    'site_file': SITE_FILE,
    'radial_src_dir': RADIAL_SRC_DIR,
    'site_frequency_dict': SITE_FREQUENCY_DICT,
    'freq_resolution': FREQ_RESOLUTION,
    'sites_to_process': ['SITE'],  # 'all' or list of specific site(s) in SITE_FREQUENCY_DICT i.e. long_range, medium_range, short_range or ['NANT', 'BLCK'] or ['HEMP']
    'proximity': None  # Set to specific km value to search for nearest radial site within that distance, or None for default set per frequency
}

# ============================================================================
# TOTALS COMPARISON PARAMETERS
# ============================================================================

# Select which totals dataset to use
TOTALS_FREQUENCY = 'long'  # Options: 'long', 'medium', 'short'

TOTALS_URLS = {
    'long': 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/codar/totals/oi/06km_realtime_qartod/realtime_maracoos_06km_totals_qartod_best.ncd',
    'medium': 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/codar/totals/oi/06km_realtime_qartod/realtime_maracoos_02km_totals_qartod_best.ncd',
    'short': 'https://tds.marine.rutgers.edu/thredds/dodsC/cool/codar/totals/oi/06km_realtime_qartod/realtime_maracoos_01km_totals_qartod_best.ncd'
}

TOTALS_PARAMS = {
    'totals_url': TOTALS_URLS['long'],  # Change to 'medium' or 'short' as needed
    'frequency': 'long',  # Must match the URL selected above
    'freq_resolution': FREQ_RESOLUTION,
    'use_interpolation': True  # Set to True for bilinear interpolation
}

# ============================================================================
# FIGURE SELECTION
# ============================================================================

# Radial comparison figures
RADIAL_FIGURES = {
    'track_map': True,              # Drifter track map with radial velocities
    'velocity_timeseries': True,    # Time padded radial velocity time series
    'difference_timeseries': True,  # Padded time series of velocity differences  
    'velocity_with_distance': True, # Velocities with red distance to vector
    'correlation_scatter': True     # Scatter plot with correlation line
}

# Totals comparison figures
TOTALS_FIGURES = {
    'track_map': True,              # Drifter track map
    'current_roses': True,          # Current roses (drifter, totals, difference)
    'uv_timeseries': True,          # Basic U and V time series
    'animation_frames': True,       # Animation frames of drifter through totals
    'correlation_scatter': True,    # U and V correlation scatter plots
    'improved_uv_timeseries': True  # Enhanced U and V time series plots
}
