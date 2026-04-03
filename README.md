# Demo: CH Link Prediction on Yeast DIP Network

Evaluates CH-based link prediction methods on the Yeast DIP network using 10% perturbation and AUCPR scoring. Both MATLAB and Python implementations are provided.

## Folder Structure

```
NSIforPPI/
├── legacy/  
├── matlab/                  # MATLAB implementation
│   ├── run_aucpr_from_precomputed.m
│   ├── CHA_linkpred_monopartite_final.m
│   ├── prediction_evaluation.m
│   ├── replace_inf_distances.m
│   └── CH_scores_mex_final.c
├── python/                  # Python implementation
│   ├── run_aucpr_from_precomputed.py
│   ├── ch_linkpred.py
│   ├── prediction_evaluation.py
│   ├── CH_scores_new_v2.cpp
│   ├── CH_scores_new_v2.h
│   ├── bindings.cpp
│   └── setup.py
├── matrix/                  # Input data
│   ├── Yeast_DIP_net.mat
│   ├── network_perturbed_10percent_GSP_GSN_Yeast_noconn.mat
│   └── list_pairs_10percent_GSP_GSN_Yeast_noconn.mat
└── README.md
README.md
```

Output is written to `table/table_auc_pr_CH_Yeast_DIP_net.xlsx`.

---

## MATLAB

### Requirements

- MATLAB R2019b or later
- A C compiler supported by MATLAB
- Statistics and Machine Learning Toolbox (for `tiedrank`)

### Setup: Compile MEX (once)

Navigate to `matlab/` and run in MATLAB:

```matlab
% macOS / Linux
mex CH_scores_mex_final.c CFLAGS='$CFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp'

% Windows (adjust MinGW path if needed)
mex C:\ProgramData\MATLAB\SupportPackages\R2020b\3P.instrset\mingw_w64.instrset\lib\gcc\x86_64-w64-mingw32\6.3.0\libgomp.a CH_scores_mex_final.c CFLAGS='$CFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp'
```

### Run

```matlab
cd /path/to/NSIforPPI/matlab
run_aucpr_from_precomputed()
```

---

## Python

### Requirements

- Python 3.8+
- A C++ compiler with OpenMP support
  - macOS: `brew install libomp llvm`
  - Linux: GCC (usually pre-installed)

### Setup: Install dependencies and compile C++ extension (once)

```bash
cd /path/to/NSIforPPI/python
pip install pybind11 scipy numpy pandas networkx openpyxl
python setup.py build_ext --inplace
```

### Run

```bash
cd /path/to/NSIforPPI/python
python run_aucpr_from_precomputed.py
```
