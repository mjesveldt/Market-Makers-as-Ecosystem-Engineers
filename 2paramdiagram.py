"""
Two-parameter diagram in the (c, R) plane for the market-maker model (Figure 2).
 
Usage
    python 2paramdiagram.py    saves figures/2paramdiagram.pdf
 
Scaled model (see the report):
    dI/dt = 1 - b (1 + cM) I M - I
    dM/dt = eps M ( R (1+cM)/(1+M) I - 1 )
 
Positive equilibria satisfy R = Phi(M), equation (12) in the report, with
    Phi(M) = (1+M) (1 + b M (1 + c M)) / (1 + c M).
 
Phi(0) = 1, and for c > 1 + b the curve Phi first dips below 1 before rising.
Its minimum R_f(c) is the saddle-node point where the liquid market and the
threshold appear or disappear together.
 
The diagram distinguishes three different regions:
    - R > 1            : liquid market only; empty market (1,0) unstable.
    - R_f(c) < R < 1   : empty and liquid market both stable, threshold in between.
    - R < R_f(c)       : empty market only.
 
eps changes neither the position nor (for eps <= 1) the stability of the
equilibria, so it does not appear in this diagram.
"""
 
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
 
# ============================================================
# 1. GENERAL PLOT SETTINGS
# ============================================================
 
# Use serif fonts and Computer Modern mathematics to make the
# figure visually consistent with the LaTeX report.
plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'cm'})
 
# Create the output directory relative to the location of this script.
# This ensures that the PDF is saved in the correct folder regardless
# of the directory from which the script is executed.
FIG_DIR = Path(__file__).parent / "figures"
FIG_DIR.mkdir(exist_ok=True)
 
# ============================================================
# 2. MODEL PARAMETERS
# ============================================================
 
# Fixed dimensionless parameter b.
# In terms of the original model parameters:
# b = beta0 / (mu * alpha)
b = 0.1
 
# Plotting ranges for the two parameters shown in the diagram.
# C_MAX determines the maximum value of c on the horizontal axis.
# R_LO and R_HI determine the minimum and maximum values of R.
C_MAX, R_LO, R_HI = 8.0, 0.5, 1.5
 
 
# ============================================================
# 3. EQUILIBRIUM CONDITION
# ============================================================
 
def Phi(M, c):
    """
    Calculate the value of R for which a given positive M
    represents an interior (positive) equilibrium.
 
    Positive equilibria satisfy:
        R = Phi(M)
 
    The function follows from setting both differential equations
    equal to zero and eliminating the variable I.
 
    M represents the dimensionless market-making presence,
    while c measures the strength of the liquidity response.
    """
    return (1 + M) * (1 + b * M * (1 + c * M)) / (1 + c * M)
 
 
# ============================================================
# 4. SADDLE-NODE BIFURCATION CURVE
# ============================================================
 
def R_fold(c):
    """
    Calculate the critical value R_f(c) at which a saddle-node
    bifurcation occurs.
 
    The saddle-node corresponds to the minimum of Phi(M, c)
    over positive values of M.
 
    Economic interpretation:
    At this critical value, a stable liquid-market equilibrium
    and an unstable threshold equilibrium appear or disappear
    together.
 
    For c <= 1 + b, Phi does not dip below 1, so there is
    no saddle-node below R = 1. In this case, R_f(c) is set to 1.
    """
    if c <= 1 + b:
        return 1.0
 
    # First, evaluate Phi over a logarithmically spaced grid.
    # A logarithmic grid allows us to investigate both small
    # and large values of M efficiently.
    M = np.logspace(-6, 3, 20000)
 
    # Find the grid point at which Phi reaches its smallest value.
    i = Phi(M, c).argmin()
 
    # Select the neighbouring grid points around this minimum.
    # These provide a smaller interval for a more accurate estimate.
    lo, hi = M[max(i - 1, 0)], M[min(i + 1, len(M) - 1)]
 
    # Refine the minimum by evaluating Phi over a denser linear grid.
    return Phi(np.linspace(lo, hi, 20000), c).min()
 
 
# ============================================================
# 5. CALCULATE THE BIFURCATION CURVE
# ============================================================
 
# Create a grid of c-values over the chosen plotting interval.
# For each value of c, calculate the corresponding saddle-node value.
c = np.linspace(0, C_MAX, 800)
Rf = np.array([R_fold(x) for x in c])
 
 
# ============================================================
# 6. COLOUR AND STYLE DEFINITIONS
# ============================================================
 
# Define separate colours for the three parameter regions.
# Green: liquid market only.
# Orange: both empty and liquid markets can be stable.
# Blue: empty market only.
GREEN, ORANGE, BLUE = "#0B7A4B", "#C2700F", "#3B6FD4"
 
# INK is used for the main bifurcation curves and annotations.
# MUTED is used for secondary elements such as axis labels and ticks.
INK, MUTED = "#1A1A1A", "#5B5B5B"
 
 
# ============================================================
# 7. CREATE THE FIGURE AND AXES
# ============================================================
 
fig, ax = plt.subplots(figsize=(6.4, 4.2))
 
# Set a white background for the entire figure.
fig.patch.set_facecolor("white")
 
 
# ============================================================
# 8. SHADE THE THREE PARAMETER REGIONS
# ============================================================
 
# Region 1: R > 1.
# The empty market is unstable, so only the liquid market is stable.
ax.fill_between(c, 1.0, R_HI, color=GREEN, alpha=0.13, lw=0)
 
# Region 2: R_f(c) < R < 1.
# Both the empty market and the liquid market are stable.
# An unstable threshold equilibrium separates their basins of attraction.
ax.fill_between(c, Rf, 1.0, color=ORANGE, alpha=0.18, lw=0)
 
# Region 3: R < R_f(c).
# No positive equilibrium exists, and the empty market is stable.
ax.fill_between(c, R_LO, Rf, color=BLUE, alpha=0.10, lw=0)
 
 
# ============================================================
# 9. PLOT THE BIFURCATION CURVES
# ============================================================
 
# First bifurcation boundary: R = 1.
# At this value, the empty-market equilibrium changes stability.
ax.axhline(
    1.0,
    color=INK,
    lw=1.6,
    zorder=3,
    label=r"empty market changes stability ($R=1$)"
)
 
# Second bifurcation boundary: R = R_f(c).
# This curve represents the saddle-node bifurcation.
# It is only present for c > 1 + b.
ax.plot(
    c[c > 1 + b],
    Rf[c > 1 + b],
    color=INK,
    lw=1.6,
    ls="--",
    zorder=3,
    label=r"liquid market appears or disappears ($R=R_{\mathrm{f}}(c)$)"
)
 
 
# ============================================================
# 10. MARK THE POINT WHERE THE TWO BOUNDARIES MEET
# ============================================================
 
# The two bifurcation curves meet at:
#     c = 1 + b, R = 1.
#
# This is the critical point at which the saddle-node curve
# emerges from the stability boundary of the empty market.
ax.plot(
    1 + b, 1.0,
    "o",
    ms=7,
    mfc="white",
    mec=INK,
    mew=1.6,
    zorder=4
)
 
# Add an annotation identifying the coordinates of this point.
ax.annotate(
    r"$(c,R)=(1+b,\,1)$",
    xy=(1 + b, 1.0),
    xytext=(1.9, 1.11),
    fontsize=9,
    color=INK,
    arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8)
)
 
 
# ============================================================
# 11. ADD ECONOMIC INTERPRETATIONS TO THE REGIONS
# ============================================================
 
# Green region:
# The empty market is unstable and a stable liquid market exists.
ax.text(
    5.6, 1.22, "liquid market only",
    ha="center", va="center",
    fontsize=9.5, color=GREEN, linespacing=1.4
)
 
# Orange region:
# Both stable equilibria exist. Which one the system approaches
# depends on the initial market-making presence.
ax.text(
    5.0, 0.83,
    "empty or liquid market\n(depends on starting size)",
    ha="center", va="center",
    fontsize=9.5, color=ORANGE, linespacing=1.4
)
 
# Blue region:
# No positive equilibrium exists, so the system converges
# towards the empty-market equilibrium.
ax.text(
    1.4, 0.68, "empty market only",
    ha="center", va="center",
    fontsize=9.5, color=BLUE
)
 
 
# ============================================================
# 12. FORMAT THE AXES
# ============================================================
 
# Add labels describing the two dimensionless parameters.
ax.set_xlabel(r"liquidity response $c$")
ax.set_ylabel(r"survival number $R$")
 
# Set the displayed parameter ranges.
ax.set_xlim(0, C_MAX)
ax.set_ylim(R_LO, R_HI)
 
# Adjust tick colours and font size for a consistent appearance.
ax.tick_params(colors=MUTED, labelsize=9)
 
# Remove the upper and right borders for a cleaner figure.
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
 
# Style the remaining axis borders.
for side in ("left", "bottom"):
    ax.spines[side].set_color(MUTED)
    ax.spines[side].set_linewidth(0.8)
 
# Add a legend explaining the two bifurcation boundaries.
ax.legend(
    loc="upper left",
    fontsize=8.5,
    frameon=True,
    framealpha=0.95,
    edgecolor="#D8D8D8"
)
 
 
# ============================================================
# 13. SAVE THE FIGURE
# ============================================================
 
# Automatically adjust the layout to avoid overlapping elements.
fig.tight_layout()
 
# Export the diagram as a PDF for inclusion in the report.
fig.savefig(
    FIG_DIR / "2paramdiagram.pdf",
    bbox_inches="tight",
    facecolor="white"
)
 
print("Saved to", FIG_DIR / "2paramdiagram.pdf")
 
 
# ============================================================
# 14. PRINT SELECTED NUMERICAL RESULTS
# ============================================================
 
# Print the fixed value of b and the point where the two
# bifurcation boundaries meet.
print("b =", b, " lines meet at c = 1+b =", 1 + b)
 
# Calculate and display selected saddle-node values.
# These values can be used as numerical references in the report.
for x in (1.2, 1.5, 2.0, 2.4, 3.0, 5.0):
    print(f"  c = {x:7.4f}   R_fold = {R_fold(x):.4f}")
 
 
# ============================================================
# 15. DISPLAY THE FIGURE
# ============================================================
 
# Open the figure in an interactive Matplotlib window.
plt.show()