#!/usr/bin/env python

import argparse
import sys
import yt
import glob
import numpy as np
import matplotlib.pyplot as plt

# Command-line arguments
parser = argparse.ArgumentParser(description="Plot the log density profile along X for particles.")

parser.add_argument('-p', action='store', required=False, type=str, dest='prefix',
                    help='path prefix [%(default)s]', default='..')
parser.add_argument('-s', action='store', required=True, type=int, dest='idx_start',
                    help='first data index')
parser.add_argument('-e', action='store', required=True, type=int, dest='idx_end',
                    help='last data index')
parser.add_argument('-d', action='store', required=False, type=int, dest='didx',
                    help='delta data index [%(default)d]', default=1)
parser.add_argument('--nbins', type=int, default=30,
                    help='number of bins for the histogram [%(default)d]')

args = parser.parse_args()

# Echo arguments
print('\nCommand-line arguments:')
print('-------------------------------------------------------------------')
print(' '.join(map(str, sys.argv)))
print('-------------------------------------------------------------------\n')

idx_start = args.idx_start
idx_end = args.idx_end
didx = args.didx
prefix = args.prefix
nbins = args.nbins

yt.enable_parallelism()

# Load specific range of datasets
ts = yt.DatasetSeries([f"{prefix}/Data_{idx:06d}" for idx in range(idx_start, idx_end + 1, didx)])

for ds in ts.piter():
    print(f"Processing density profile: {ds}")
    ad = ds.all_data()
    x = ad["particle_position_x"].to("kpc").value

    # Bin edges and centers
    x_bins = np.linspace(x.min(), x.max(), nbins + 1)
    x_centers = 0.5 * (x_bins[:-1] + x_bins[1:])

    # Count particles in each X bin
    counts, _ = np.histogram(x, bins=x_bins)
    counts = np.where(counts == 0, 0.1, counts)  # Avoid log(0)

    # Plot
    plt.figure(figsize=(8, 5))
    plt.plot(x_centers, counts, drawstyle="steps-mid", color="darkred")
    plt.xlabel("X [kpc]")
    plt.ylabel("Particle Count (log scale)")
    plt.yscale("log")
    plt.ylim(10, 3e2)  # You can make this an argument too if needed
    plt.title(f"Log Particle Density Profile Along X\n{ds.basename}, Time = {ds.current_time.to('Myr').value:.2f} Myr")
    plt.tight_layout()
    plt.savefig(f"{ds.basename}_density_profile_x_log.png", dpi=300)
    plt.close()
