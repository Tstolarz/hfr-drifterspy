# Drifter-HFR Comparison Tool

A Python tool for comparing drifter trajectory data with High-Frequency Radar (HFR) measurements, including both radial velocities and total current fields.

## Overview

This tool processes drifter GPS trajectory data and compares it against HF Radar observations to validate radar measurements and analyze ocean current patterns. It supports:

- **Radial Velocity Comparison**: Compare drifter-derived radial velocities with HFR radial measurements from individual radar sites
- **Total Velocity Comparison**: Compare drifter velocities with gridded HFR total current fields
- **Quality Control**: Apply multiple QC filters to drifter data (rate of change, standard deviation, distance thresholds)
- **Visualization**: Generate comprehensive figures including time series, scatter plots, current roses, and animations

## Requirements

### Python Dependencies

```bash
numpy
pandas
xarray
matplotlib
cartopy
pyproj
windrose
cmocean
cool_maps
hfradarpy
oceans
```

### Installation

Install required packages using pip:

```bash
pip install numpy pandas xarray matplotlib cartopy pyproj windrose cmocean
pip install cool_maps hfradarpy oceans
```

### Data Requirements

1. **Drifter Data**: NetCDF file containing drifter GPS positions with:
   - `time`: timestamp dimension
   - `lat`: latitude coordinates
   - `lon`: longitude coordinates
   - `drifter_id`: drifter identifier

2. **HFR Radial Data**: Directory structure with radial files organized as:
   ```
   radials/
   ├── SITE1/
   │   ├── YYYY_MM/
   │   │   └── RDLm_SITE1_YYYY_MM_DD_HHMM.ruv
   ├── SITE2/
   │   └── ...
   ```

3. **HFR Site Locations**: Excel file (`HFR_sites.xlsx`) with site metadata (see format below)

4. **HFR Totals Data**: THREDDS URL pointing to gridded total current dataset with `u` and `v` variables

## Quick Start

### 1. Set Up Your Configuration

Copy and edit the default configuration file:

```bash
cp config_default.py config_myproject.py
```

Edit `config_myproject.py` to specify:
- Path to your drifter NetCDF file
- Output directory for results
- Path to HFR site locations file
- Directory containing radial files
- THREDDS URLs for totals data (if using totals comparison)

### 2. Update the Run Script

Edit `run_drifter_comparison.py` line 28 to point to your config:

```python
import config_myproject as config
```

### 3. Run the Analysis

```bash
python run_drifter_comparison.py
```

The script will automatically:
- Load and QC your drifter data
- Find matching HFR observations
- Generate comparison figures
- Save results to your specified output directory

## Configuration Guide

### File Paths

```python
DRIFTER_FILE = '/path/to/your/drifter_file.nc'
DRIFTER_ID = 'your_drifter_id'
OUTPUT_DIR = '/path/to/your/output'
SITE_FILE = '/path/to/HFR_sites.xlsx'
RADIAL_SRC_DIR = '/path/to/radials'
```

### Comparison Selection

```python
RUN_RADIAL_COMPARISON = True   # Compare with radial velocities
RUN_TOTALS_COMPARISON = True   # Compare with total currents
```

### Quality Control Parameters

```python
QC_PARAMS = {
    'apply_roc': True,          # Apply rate of change test
    'roc_threshold': 0.15,      # Maximum velocity change (m/s)
    'apply_std': True,          # Apply standard deviation test
    'std_levels': 3,            # Number of std deviations
    'apply_dist': True,         # Apply distance threshold test
    'dist_threshold': 10000     # Maximum distance between points (m)
}
```

### Radial Comparison Parameters

```python
RADIAL_PARAMS = {
    'site_file': SITE_FILE,
    'radial_src_dir': RADIAL_SRC_DIR,
    'site_frequency_dict': SITE_FREQUENCY_DICT,
    'freq_resolution': FREQ_RESOLUTION,
    'sites_to_process': 'all',   # or ['SITE1', 'SITE2'] or 'long'/'medium'/'short'
    'proximity': None            # km, or None for default
}
```

### Totals Comparison Parameters

```python
TOTALS_PARAMS = {
    'totals_url': TOTALS_URLS['long'],
    'frequency': 'long',               # 'long', 'medium', or 'short'
    'freq_resolution': FREQ_RESOLUTION,
    'use_interpolation': True          # Bilinear interpolation
}
```

### Figure Selection

Control which figures are generated:

```python
# Radial comparison figures
RADIAL_FIGURES = {
    'track_map': True,               # Drifter track with radial velocities
    'velocity_timeseries': True,     # Time series of velocities
    'difference_timeseries': True,   # Time series of differences
    'velocity_with_distance': True,  # Velocities with distance overlay
    'correlation_scatter': True      # Scatter plot with correlation
}

# Totals comparison figures
TOTALS_FIGURES = {
    'track_map': True,               # Drifter track map
    'current_roses': True,           # Current rose diagrams
    'uv_timeseries': True,           # U/V component time series
    'animation_frames': True,        # Animation frame sequence
    'correlation_scatter': True,     # U/V correlation scatter
    'improved_uv_timeseries': True   # Enhanced time series plots
}
```

## HFR Site Configuration

### Adding New Sites

To add a new HFR site, edit the `HFR_sites.xlsx` file with the following columns:

| Column | Type    | Description                           | Example      |
|--------|---------|---------------------------------------|--------------|
| site   | string  | 4-letter site code (uppercase)        | HEMP         |
| range  | string  | Frequency range: 'long'/'medium'/'short' | long      |
| lat    | float   | Latitude in decimal degrees           | 40.969333    |
| lon    | float   | Longitude in decimal degrees          | -72.123700   |

**Example entry:**
```
site,range,lat,lon
HEMP,long,40.969333,-72.123700
SEAB,medium,39.753317,-73.991567
PORT,short,38.329800,-75.089367
```

### Site Frequency Categories

The default configuration includes three frequency ranges:

- **Long Range (5-6 km resolution)**: NANT, BLCK, AMAG, MRCH, HEMP, HOOK, LOVE, BRIG, WILD
- **Medium Range (2-3 km resolution)**: SEAB, BRAD, SPRK, HLGT, BRMR, RATH, WOOD, CMPT
- **Short Range (1 km resolution)**: SILD, OLDB, PORT, HLPN, LEWE, CAPE

When adding a new site, include it in the appropriate frequency list in your config file:

```python
SITE_FREQUENCY_DICT = {
    'long': ['NANT', 'BLCK', 'AMAG', 'MRCH', 'HEMP', 'HOOK', 'LOVE', 'BRIG', 'WILD', 'YOUR_SITE'],
    'medium': ['SEAB', 'BRAD', 'SPRK', 'HLGT', 'BRMR', 'RATH', 'WOOD', 'CMPT'],
    'short': ['SILD', 'OLDB', 'PORT', 'HLPN', 'LEWE', 'CAPE']
}
```

## Output

The tool creates a timestamped output directory containing:

```
{DRIFTER_ID}_{TIMESTAMP}_hfr_drifter_comparison/
├── Radials/
│   ├── SITE1/
│   │   ├── track_map.png
│   │   ├── velocity_timeseries.png
│   │   ├── difference_timeseries.png
│   │   ├── velocity_with_distance.png
│   │   └── correlation_scatter.png
│   └── SITE2/
│       └── ...
└── Totals/
    ├── track_map.png
    ├── current_roses.png
    ├── uv_timeseries.png
    ├── correlation_scatter.png
    ├── improved_uv_timeseries.png
    └── animation_frames/
        ├── frame_001.png
        ├── frame_002.png
        └── ...
```

## How It Works

### Radial Comparison Workflow

1. **Load Drifter Data**: Import drifter NetCDF file and calculate velocities from GPS positions
2. **Apply QC**: Filter drifter data using configured quality control thresholds
3. **Find Radial Files**: For each drifter position, locate corresponding HFR radial files based on time and proximity to radar sites
4. **Extract Radial Velocities**: Read radial files and find nearest velocity vector to drifter position
5. **Calculate Drifter Radials**: Project drifter velocity onto radar line-of-sight direction
6. **Compare & Visualize**: Generate comparison plots and statistics

### Totals Comparison Workflow

1. **Load Drifter Data**: Import and QC drifter data
2. **Load HFR Totals**: Access gridded total current fields via THREDDS
3. **Extract at Drifter Positions**: Get HFR velocities at drifter locations (nearest neighbor or bilinear interpolation)
4. **Compare & Visualize**: Generate U/V component comparisons, current roses, and animations

## Troubleshooting

### "ERROR: Default configuration file detected!"

You're using `config_default.py` without modification. Create your own config file and update the import in `run_drifter_comparison.py`.

### "No radial data found for comparison"

Check that:
- `RADIAL_SRC_DIR` points to the correct directory
- Radial files follow the naming convention: `RDLm_SITE_YYYY_MM_DD_HHMM.ruv`
- Drifter time range overlaps with available radial data
- `proximity` setting allows the drifter to be within range of at least one site

### "No totals data found for comparison"

Verify:
- THREDDS URL is accessible and correct
- Drifter time range overlaps with totals data coverage
- Drifter positions are within the totals grid domain

