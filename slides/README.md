# Slides

**[Parameter Monitoring with Adaptive Observers (PDF, 26 pages)](MiniCurso_AdaptiveObservers.pdf)**
&nbsp;·&nbsp; source [`MiniCurso_AdaptiveObservers.tex`](MiniCurso_AdaptiveObservers.tex)

| Section | What it covers |
|---|---|
| Motivation | why an adaptive observer, the 4-step recipe |
| Step 1 · Model | virus–T-cell model (1)–(2), the "mean" patient, what is measured |
| Step 2 · Can we? | observability ($\det = -c_T V$) and identifiability |
| Step 3 · Rewrite | change of coordinates (13), model (14), the state-affine recipe (15), Examples 1 and 2 |
| Step 4 · Observer | the adaptive observer (17), how to read it, assumptions A1–A3 |
| Persistent excitation | PE as a minimal-eigenvalue condition (Aranovskiy et al., 2017), the PE meter $\lambda_{\min}(\Gamma)$ |
| Practice | tuning knobs, from equations to code, results, noise and sampling |
| Wrap-up | exercises, take-home recipe, references, open projects (PAE / theses) |

## Build

```sh
python make_figures.py                      # figures/ from the Python examples
latexmk -pdf MiniCurso_AdaptiveObservers.tex
```

Needs Beamer with the metropolis theme, TikZ and `listings`.
