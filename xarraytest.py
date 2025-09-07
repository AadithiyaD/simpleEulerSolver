import xarray as xr
from pythia_datasets import DATASETS

filepath = DATASETS.fetch('NARR_19930313_0000.nc')

ds = xr.open_dataset(filepath)

temperature_data = ds['Temperature_isobaric'].sel(x="-3292.0078")

print(temperature_data)