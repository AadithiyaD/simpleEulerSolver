import pandas as pd
import xarray as xr
import glob
import matplotlib.pyplot as plt
from matplotlib import ticker
import sodshock
import re


# Get exact solution from sodshock module
# This can be cached
gamma = 1.4
dustFrac = 0.0
npts = 1000
t = 0.2
left_state = (1, 1, 0)
right_state = (0.1, 0.125, 0.)

positions, regions, exact_solution = sodshock.solve(left_state=left_state,
                                                    right_state=right_state, geometry=(0., 1., 0.5), t=t,
                                                    gamma=gamma, npts=npts, dustFrac=dustFrac)


def extract_pattern(fname):
    r''' Extract the timestep from the file name. Regex explanation
    _ => literal underscore
    (\d+) => () is a capture group \d means digits + means one or more
    \. means dot. We're escaping the . bc it normally means "everything"
    $ Means end of search pattern
        The entire thing is in single quotes '''

    match = re.search(r'_(\d+)\.csv$', fname)

    # 1 means just the capture group so we get just 000126 for example
    return int(match.group(1)) if match else -1


def load_solver_data(filename_pattern="solnData/solution_*_*.csv"):
    """Takes csv files and loads data into an xarray dataset,
    and returns a combined dataset with dimensions time and x

    default arg:
    filename_pattern= solnData/solution_*_*.csv
    looks at all files inside solnData that match the pattern solution_*_*.csv"""

    csv_files = sorted(glob.glob(filename_pattern))

    timesteps = []
    datasets = []

    for file_path in csv_files:
        # get timestep from filename
        timesteps.append(extract_pattern(file_path))

        # Read csv
        df = pd.read_csv(file_path)

        # Convert to xarray dataset
        ds = df.set_index('x').to_xarray()
        datasets.append(ds)

    # Combine along the time dimension so we have one big array where each layer \
    # Is a timestep
    combined = xr.concat(datasets, dim='time')
    combined['time'] = timesteps

    return combined

# Uncomment for first time, afterwards, use cached data
# solver_data = load_solver_data()

# Cache data so that you dont have to parse csv again and again
# solver_data.to_netcdf("solver_data_cache.nc")


# Load data from local cache
solver_data = xr.load_dataset('solver_data_cache.nc')

# Plot just the final timestep solutions


def visualize():
    """Function to plot our graphs and postprocess"""

    # Put everything into shorter variable for ease of use
    plot_timestep_data = solver_data.sel(time=solver_data.time.max().values)
    plot_x = plot_timestep_data.x
    plot_u = plot_timestep_data.u
    plot_rho = plot_timestep_data.rho
    plot_p = plot_timestep_data.p

    # Plot velocity
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(plot_x, plot_u, label="Osher Flux limiter",
            color="purple", linestyle='-')
    ax.plot(exact_solution['x'], exact_solution['u'],
            label="Exact solution", color='black', linestyle='-')

    ax.grid(True)
    ax.set(xlabel='x', ylabel=r'Velocity [m/s]', title='Velocity Plot')
    ax.legend()

    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))

    # Plot pressure
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(plot_x, plot_p, label="Osher Flux limiter",
            color="purple", linestyle='-')
    ax.plot(exact_solution['x'], exact_solution['p'],
            label="Exact solution", color='black', linestyle='-')

    ax.grid(True)
    ax.set(xlabel='x', ylabel=r'Pressure [Pa]', title='Pressure plot')
    ax.legend()

    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))

    # Plot desnity
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(plot_x, plot_rho, label="Osher Flux limiter",
            color="purple", linestyle='-')
    ax.plot(exact_solution['x'], exact_solution['rho'],
            label="Exact solution", color='black', linestyle='-')

    ax.grid(True)
    ax.set(xlabel='x', ylabel=r'Desnity [kg/m^-3]', title='Density plot')
    ax.legend()

    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))

    plt.show()

    return 0


visualize()
