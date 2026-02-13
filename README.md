# SPA Benchmark Data

This repository contains the raw data and code for benchmarking the Separable Pair Ansatz (SPA) method on hydrogen chains ($H_n$).


## Requirements

Python 3.10. 
Install the necessary packages:

```bash
pip install tequila-basic sunrise spafastprototype matplotlib
```
Note: Some packages (like `sunrise` and `spafastprototype`) may require specialized installation. Refer to their documentation.

## Generating Data 

### 1. Dissociation Curve
Generate a CSV file for a given linear $H_n$ molecule.

```python
from main_Hn import run_dissociation

run_dissociation(n=4, max_iter=36, d_min=0.5, d_max=4.0, get_fci=True, get_var=True)
```

- `n` – Number of H atoms
- `max_iter` – Number of points along the dissociation curve
- `d_min`, `d_max` – Minimum and maximum interatomic distances (Å)

Output CSV file: `results_h{n}.csv`

Each row of the generated file contains:
- `distance` – Interatomic distance
- `spa` – Optimized SPA energy
- `fci` – Reference FCI ground-state energy (if `get_fci=True`)
- `fid` – Fidelity between SPA and FCI states (if `get_fci=True`)
- `var` – Energy variance of the SPA state (if `get_var=True`)

### 2. Single-Point Scaling
Generate data for increasing $H_n$ chain length at a fixed interatomic distance.

```python
from main_Hn import run_scaling

run_scaling(n_min=2, n_max=30, distance=1.0)
```
Output CSV files: 
- `results_vs_n.csv` – SPA/FCI/fidelity/variance vs. n
- `timing_vs_n.csv` – Total runtime and individual timing contributions vs. n

## Notes on Computational Cost

- **FCI calculations** become computationally demanding for larger systems.  
  For example, a single-point FCI calculation for $H_{14}$ may take approximately **30 minutes to 1 hour**, depending on hardware.  
  For larger chains, FCI becomes prohibitively expensive and is therefore not recommended.

- **Energy variance calculations** are also computationally intensive.  
  Even for $H_{10}$, the variance evaluation is very slow (**3/4 hours**) and is therefore also not recommended beyond this system size.

- In the provided scaling script, FCI and variance calculations are automatically disabled for larger systems to avoid excessive runtimes.


## Reproducing Plots

```python
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv(f"timing_vs_n.csv")
for col in data.columns:
    if col != "n" and col != "total_time":
        plt.plot(data["n"], data[col], marker="o", label=col)

plt.xlabel("Number of hydrogens (n)")
plt.ylabel("Time (s)")
plt.yscale("log")   
plt.legend()
plt.tight_layout()
plt.xticks(df["n"])
# plt.savefig("times_vs_n.pdf")
plt.show()
```
```python
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv(f"results_vs_n.csv")
data = data[data["n"] <= 14].copy()
plt.plot(data["n"], 1-data["fidelity"], marker="o", label="1 - Fidelity")
plt.plot(data["n"], data["variance"], marker="o", label="Variance (eH)")
plt.plot(data["n"], data["spa"]-data["fci"], marker="o", label="Error (eH)")
plt.xlabel("Number of hydrogens (n)")
plt.legend()
plt.tight_layout()
plt.xticks(df["n"])
# plt.savefig(f"results_vs_n.pdf")
plt.show()
```

```python
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv(f"results_h6.csv")
plt.plot(data["distance"], data["spa"], "o-", label=f"SPA")
plt.plot(data["distance"], data["fci"], "o-", label=f"FCI")
plt.xlabel("Interatomic distance (Å)")
plt.ylabel("Energy (eH)")
plt.legend()
plt.grid(True)
plt.tight_layout()
# plt.savefig("dissociation_curves_h6.pdf")
plt.show()
```

```python
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv(f"results_h6.csv")
plt.plot(data["distance"], data["spa"]-data["fci"], "o-")
plt.xlabel("Interatomic distance (Å)")
plt.ylabel("Error (eH)")
plt.legend()
plt.grid(True)
plt.tight_layout()
# plt.savefig("error.pdf")
plt.show()
```
Replace `data["spa"]-data["fci"]` by `data["fid"]` or `data["var"]`.
