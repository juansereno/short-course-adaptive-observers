# Python code

| File | What it does |
|---|---|
| [`ex1_cT.py`](ex1_cT.py) | Example 1: estimates $\gamma = c_T S_T$ and reconstructs the T cells |
| [`ex2_cT_r.py`](ex2_cT_r.py) | Example 2: estimates $\gamma$ and $r$ |
| [`ex1_cT.ipynb`](ex1_cT.ipynb), [`ex2_cT_r.ipynb`](ex2_cT_r.ipynb) | the same code as executed notebooks, one cell per step, plots saved |
| [`make_notebooks.py`](make_notebooks.py) | rebuilds the notebooks from the scripts |
| [`requirements.txt`](requirements.txt) | NumPy, SciPy, Matplotlib (+ JupyterLab for the notebooks) |

```sh
python -m pip install -r requirements.txt
python ex1_cT.py
python ex2_cT_r.py
```

Each script is a single file organised in numbered steps (parameters, patient model,
measurements, state-affine matrices, tuning, observer ODE, initial conditions, solve,
post-processing, PE meter, plots). It uses explicit column vectors and `@` products,
so it reads line by line like the MATLAB version and gives the same numbers.

The `.py` scripts are the source. After editing them, rebuild the notebooks with
`python make_notebooks.py --execute`.
