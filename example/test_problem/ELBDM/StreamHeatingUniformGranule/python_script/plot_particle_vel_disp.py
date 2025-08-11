import argparse
import sys
import yt
import numpy as np
import matplotlib.pyplot as plt

# ------------------ CLI ------------------
parser = argparse.ArgumentParser(description='Velocity dispersion components vs time in subplots')
parser.add_argument('-p', type=str, default='..', dest='prefix', help='path prefix [%(default)s]')
parser.add_argument('-s', type=int, required=True, dest='idx_start', help='first data index')
parser.add_argument('-e', type=int, required=True, dest='idx_end', help='last data index')
parser.add_argument('-d', type=int, default=1, dest='didx', help='delta data index [%(default)d]')
parser.add_argument('--out', type=str, default='sigma_components_subplots.png',
                    help='output filename for the dispersion plot')
args = parser.parse_args()

print('\nCommand-line arguments:')
print('-------------------------------------------------------------------')
print(' '.join(map(str, sys.argv)))
print('-------------------------------------------------------------------\n')

# Dataset series
yt.enable_parallelism()
ts = yt.DatasetSeries([f"{args.prefix}/Data_{i:06d}" for i in range(args.idx_start, args.idx_end+1, args.didx)])

# ---- helpers ----
def weighted_std(v, w):
    wsum = w.sum()
    mean = (w * v).sum() / wsum
    var  = (w * (v - mean)**2).sum() / wsum
    return var**0.5

def component_sigmas(ad):
    vx = ad['particle_velocity_x'].to('km/s')
    vy = ad['particle_velocity_y'].to('km/s')
    vz = ad['particle_velocity_z'].to('km/s')
    m  = ad['particle_mass']
    return weighted_std(vx, m), weighted_std(vy, m), weighted_std(vz, m)

# storage
storage = {}

# ---- main loop ----
for sto, ds in ts.piter(storage=storage):
    ad = ds.all_data()
    sx, sy, sz = component_sigmas(ad)
    t = ds.current_time.to('Myr')
    sto.result = dict(time=t, sigma_x=sx, sigma_y=sy, sigma_z=sz)

# ---- gather & plot subplots ----
if yt.is_root():
    recs = sorted(storage.values(), key=lambda r: r['time'].to_value())
    times = np.array([r['time'].to_value('Myr') for r in recs])
    sx    = np.array([r['sigma_x'].to_value('km/s') for r in recs])
    sy    = np.array([r['sigma_y'].to_value('km/s') for r in recs])
    sz    = np.array([r['sigma_z'].to_value('km/s') for r in recs])

    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

    axes[0].plot(times, sx, marker='o', color='r')
    axes[0].set_ylabel(r'$\sigma_x$ [km/s]')
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(times, sy, marker='s', color='g')
    axes[1].set_ylabel(r'$\sigma_y$ [km/s]')
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(times, sz, marker='^', color='b')
    axes[2].set_ylabel(r'$\sigma_z$ [km/s]')
    axes[2].set_xlabel('Time [Myr]')
    axes[2].grid(True, alpha=0.3)

    fig.suptitle('Number density-weighted Velocity Dispersion Components vs Time', fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(args.out, dpi=150)
    print(f"Saved subplot figure to {args.out}")

    # optional CSV
    try:
        import pandas as pd
        df = pd.DataFrame({
            'time_Myr': times,
            'sigma_x_kms': sx,
            'sigma_y_kms': sy,
            'sigma_z_kms': sz
        })
        csv_path = args.out.rsplit('.', 1)[0] + '.csv'
        df.to_csv(csv_path, index=False)
        print(f"Wrote data to {csv_path}")
    except Exception:
        pass
