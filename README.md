# BaffinBaySIT
Code created to determine sea ice thickness in Baffin Bay (2003-2020) using different snow depth products and redistribution methods. Used for Glissenaar et al. (2021) available in The Cryosphere (https://doi.org/10.5194/tc-15-4909-2021).


## Repository structure:
```
BaffinBaySIT/
│
├─ src/
│  ├─ analyze/                                  | analyze sea ice thickness records Baffin Bay
│  │  ├─ march_asym_allmethods.py
│  │  ├─ march_mean_cs2_methods.py
|  |  ├─ march_mean_icesat_diff_methods.py
|  |  ├─ march_timeseries_allmethods.py
│  │  └─ mean_std_all_sats_all_methods.py
│  │     
│  ├─ data/                                     | preprocessing data
|  |  ├─ make_sit_netcdf.py
|  |  ├─ open_cs2_cci.py
│  │  └─ open_icesat2.py
|  |
|  ├─ functions/                                
|  |  ├─ func_easegrid.py
|  |  ├─ func_ieasegrid.py
|  |  ├─ func_reproj.py
|  |  ├─ func_scatteredInterpolant.py
|  |  ├─ func_smooth.py
|  |  ├─ func_snowdepth.py
│  │  └─ sit_functions.py
|  |
│  └─ sit/                                      | create SIT products
│     ├─ cryosat2_allmethods.py
│     ├─ cryosat2_cci_allmethods.py
│     ├─ envi_allmethods.py
│     ├─ freeboard_ICESat.py
│     ├─ icesat_allmethods.py
│     └─ icesat2_allmethods.py
│
├─ LICENSE
└─ README.md

```

