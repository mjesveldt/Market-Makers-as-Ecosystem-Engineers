# Market Makers as Ecosystem Engineers: an ODE analysis of survival and collapse

J.T. van Dijk, S.A.M. van Diermen, M.J. Esveldt
Group project for *Modelling with ODEs and PDEs*, Utrecht University, 2026.

📄 **Read the report: [report.pdf](report.pdf)**

## Summary

Market makers earn the bid–ask spread, but by competing with each other they make that spread smaller. A smaller spread attracts more trading, but leaves less profit per trade. Like ecosystem engineers in biology, market makers change their own environment in a way that changes the availability of resources.

We model trading interest $I$ and market-making presence $M$ with a system of two ODEs, and ask two questions:

1. When can a market-making industry sustain itself?
2. Can it collapse abruptly?

The main results:

- A dimensionless survival number $R$ decides whether a first market maker in an empty market grows ($R>1$) or shrinks ($R<1$).
- For $R>1$, a liquid market is the only stable state.
- If trading volume responds strongly enough to tighter spreads ($c>1+b$), an empty and a liquid market can both be stable for $R<1$. The industry then survives only if it is already large enough, and it can collapse abruptly, either when conditions slowly get worse or after a sudden shock. Restoring the conditions afterwards does not bring it back.

## The model

After nondimensionalization (Section 2 of the report) the system is

$$
\frac{d\tilde I}{d\tilde t} = 1 - b\,(1+c\tilde M)\,\tilde I\tilde M - \tilde I,
\qquad
\frac{d\tilde M}{d\tilde t} = \varepsilon\,\tilde M\left(R\,\frac{1+c\tilde M}{1+\tilde M}\,\tilde I - 1\right),
$$

| Symbol | Meaning |
|---|---|
| $R$ | Survival number: does a first market maker in an empty market grow or shrink? |
| $c$ | Liquidity response: how strongly trading rises when the spread tightens |
| $b$ | How fast market makers use up trading interest, relative to how fast it expires |
| $\varepsilon$ | How fast market makers leave, relative to how fast trading interest expires ($\varepsilon \ll 1$) |

All figures use $b = 0.1$ and $\varepsilon = 0.01$.

## Files

| File | What it is |
|---|---|
| `report.pdf` | The report |
| `sliders.py` | Phase portraits and bifurcation diagrams (Figure 1), interactive or saved to `figures/` |
| `2paramdiagram.py` | Two-parameter diagram in the $(c,R)$ plane (Figure 2) |
| `requirements.txt` | Python packages needed |

The `figures/` folder is not included; the two scripts create it when you run them.

## How to reproduce the figures

You need Python 3. Install the packages with

```
pip install -r requirements.txt
```

Then run

```
python sliders.py --save     # saves the 12 panels of Figure 1 to figures/
python 2paramdiagram.py      # saves Figure 2 to figures/2paramdiagram.pdf
```

Running `python sliders.py` without `--save` opens an interactive window with sliders for the original (dimensional) parameters. It shows the phase portrait, the bifurcation diagram and a table of the equilibria with their type.

File names of Figure 1 give $R$ and $c$ without the decimal point: `bif075.pdf` is the bifurcation diagram for $R = 0.75$, and `phaseR075c30.pdf` is the phase portrait for $R = 0.75$ and $c = 3$.

## Use of AI

Claude Opus 5.5 (Anthropic) was used during this project. It helped us to write and debug the Python code for the phase portraits and bifurcation diagrams of Figure 1, as well as the (c, R) plane of Figure 2. We also used Claude for formatting in LaTeX as well as grammar and sentence clarity checks throughout. No AI tool was used for the model derivation, the mathematical analysis, or the interpretation of results. The code made by AI was tested extensively.

## Licence

MIT, see [LICENSE](LICENSE).
