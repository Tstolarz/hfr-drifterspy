#!/usr/bin/env python
"""
run_drifter_comparison.py
Drifter-HFR Comparison Tool
Main script for running drifter vs HF Radar comparisons

Usage:
    python run_drifter_comparison.py
    
This script allows users to:
1. Compare drifter data with HF Radar radial velocities
2. Compare drifter data with HF Radar total velocities
3. Apply various QC methods to drifter data
4. Generate customizable figures and analysis outputs
"""

import os
import sys
from drifter_hfr_functions import (
    process_radial_comparison,
    process_totals_comparison
)

# ============================================================================
# IMPORT YOUR CONFIGURATION FILE HERE
# ============================================================================

import config_default as config  # User configuration file defined here where 

# ============================================================================
# Script Runs From Here...
# ============================================================================

def check_default_config():
    """
    Check if user is running with the unchanged default config file.
    Returns True if default placeholders are detected, False otherwise.
    """
    default_indicators = [
        '/path/to/' in config.DRIFTER_FILE,
        '/path/to/' in config.OUTPUT_DIR,
        '/path/to/' in config.SITE_FILE,
        '/path/to/' in config.RADIAL_SRC_DIR,
        config.DRIFTER_ID == 'your_drifter_id',
        'SITE' in config.RADIAL_PARAMS.get('sites_to_process', []),
        'link-to-your-thredds-file-link' in str(config.TOTALS_URLS.get('long', ''))
    ]

    return any(default_indicators)

def main():
    """Main execution function"""

    print("=" * 70)
    print("DRIFTER-HFR COMPARISON TOOL")
    print("=" * 70)
    print("Started run at:", os.popen('date').read().strip())

    # Check if user is using unchanged default config
    if check_default_config():
        print("\n" + "!" * 70)
        print("ERROR: Default configuration file detected!")
        print("!" * 70)
        print("\nIt appears you are using the default config file with placeholder values.")
        print("Please edit your configuration file to specify your actual datasets and paths.")
        print("\nRequired updates:")
        print("  - DRIFTER_FILE: Path to your drifter data file")
        print("  - DRIFTER_ID: Your specific drifter ID")
        print("  - OUTPUT_DIR: Directory for output files")
        print("  - SITE_FILE: Path to HFR sites Excel file")
        print("  - RADIAL_SRC_DIR: Directory containing radial data")
        print("  - RADIAL_PARAMS['sites_to_process']: Specific sites to process")
        print("  - TOTALS_URLS: Update if using custom THREDDS URLs")
        print("\nPlease update your config file and try again :).")
        print("=" * 70)
        sys.exit(1)

    # Check if drifter file exists
    if not os.path.exists(config.DRIFTER_FILE):
        print(f"ERROR: Drifter file not found: {config.DRIFTER_FILE}")
        print("Please update DRIFTER_FILE path in the configuration section")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(config.OUTPUT_DIR, exist_ok=True)
    
    # Print configuration
    print("\nConfiguration:")
    print(f"  Drifter file: {config.DRIFTER_FILE}")
    print(f"  Drifter ID: {config.DRIFTER_ID}")
    print(f"  Output directory: {config.OUTPUT_DIR}")
    print(f"  Run radial comparison: {config.RUN_RADIAL_COMPARISON}")
    print(f"  Run totals comparison: {config.RUN_TOTALS_COMPARISON}")
    
    print("\nQC Parameters:")
    for key, value in config.QC_PARAMS.items():
        print(f"  {key}: {value}")
    
    # Run radial comparison
    if config.RUN_RADIAL_COMPARISON:
        print("\n" + "-" * 70)
        print("RUNNING RADIAL COMPARISON")
        print("-" * 70)
        
        # Check if required files exist
        if not os.path.exists(config.SITE_FILE):
            print(f"ERROR: Site file not found: {config.SITE_FILE}")
            print("Please update SITE_FILE path in the configuration section")
        elif not os.path.exists(config.RADIAL_SRC_DIR):
            print(f"ERROR: Radial directory not found: {config.RADIAL_SRC_DIR}")
            print("Please update RADIAL_SRC_DIR path in the configuration section")
        else:
            try:
                results = process_radial_comparison(
                    config.DRIFTER_FILE,
                    config.QC_PARAMS,
                    config.RADIAL_PARAMS,
                    config.RADIAL_FIGURES,
                    config.OUTPUT_DIR,
                    config.DRIFTER_ID
                )
                
                if results['status'] == 'success':
                    print(f"\nRadial comparison completed successfully!")
                    print(f"Sites processed: {', '.join(results['sites_processed'])}")
                    print(f"Output directory: {results['output_dir']}")
                    
                    # Print figures created per site
                    for site, figures in results['figures_created'].items():
                        print(f"\n  {site}:")
                        for fig in figures:
                            print(f"    - {fig}")
                else:
                    print("\nNo radial data found for comparison")
                    print(f"Example radial that wasn't found: {results.get('example_missing_radial', 'N/A')}")
                    
            except Exception as e:
                print(f"\nERROR in radial comparison: {str(e)}")
                import traceback
                traceback.print_exc()
    
    # Run totals comparison
    if config.RUN_TOTALS_COMPARISON:
        print("\n" + "-" * 70)
        print("RUNNING TOTALS COMPARISON")
        print("-" * 70)
        
        try:
            results = process_totals_comparison(
                config.DRIFTER_FILE,
                config.QC_PARAMS,
                config.TOTALS_PARAMS,
                config.TOTALS_FIGURES,
                config.OUTPUT_DIR,
                config.DRIFTER_ID
            )
            
            if results['status'] == 'success':
                print(f"\nTotals comparison completed successfully!")
                print(f"Output directory: {results['output_dir']}")
                print(f"Figures created:")
                for fig in results['figures_created']:
                    print(f"  - {fig}")
            else:
                print("\nNo totals data found for comparison")
                
        except Exception as e:
            print(f"\nERROR in totals comparison: {str(e)}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "=" * 70)
    print("PROCESSING COMPLETE")
    print("=" * 70)
    print("Finished run at:", os.popen('date').read().strip())

if __name__ == "__main__":
    main()