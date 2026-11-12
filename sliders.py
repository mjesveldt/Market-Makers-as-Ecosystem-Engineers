"""
Market makers as ecosystem engineers: interactive phase portrait and
bifurcation diagram, plus the figures used in the report.
 
Usage
    python sliders.py          interactive window with sliders
    python sliders.py --save   save the 12 report figures (PDF) to ./figures
 
Report figure names (R and c without the decimal point):
    bif075.pdf                 bifurcation diagram for R = 0.75
    phaseR075c30.pdf           phase portrait for R = 0.75 and c = 3
 
Scaled model (see the report):
    dI/dt = 1 - b (1 + cM) I M - I
    dM/dt = eps M ( R (1+cM)/(1+M) I - 1 )
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
 
# LaTeX-like fonts, so the figures match the report
plt.rcParams.update({'font.size': 18, 'font.family': 'serif',
                     'mathtext.fontset': 'cm'})
 
# ---------- Parameters: (name, min, max, initial value) ----------
# delta = 0.01 with mu = 1 gives eps = 0.01, as in the report.
PARS = [("I0", 0, 2, 1), ("mu", 0, 2, 1), ("beta0", 0, 1, 0.2),
        ("c", 0, 5, 2), ("s0", 0, 0.5, 0.1), ("alpha", 0, 5, 2),
        ("gamma", 0, 5, 0.8), ("delta", 0, 0.05, 0.01), ("Mmax", 1, 10, 5)]  # Mmax = plot range only
 
# Settings for the report figures
FIG_DIR = Path(__file__).parent / "figures"
B_REPORT, EPS_REPORT = 0.1, 0.01
R_VALUES = {0.75: "075", 1.0: "100", 1.25: "125"}  # value of R -> tag in file name
C_VALUES = {0.5: "05", 2.0: "20", 3.0: "30"}       # value of c -> tag in file name
MMAX_REPORT, IMAX_REPORT, CMAX_REPORT = 3, 1.2, 30
 
# ---------- Model ----------
def nondim(I0, mu, beta0, c, s0, alpha, gamma, delta):
    """Dimensionless groups b, eps and R, as defined in the report."""
    b   = beta0 / (mu * alpha)
    eps = delta / mu
    R   = gamma * s0 * beta0 * I0 / (mu * delta)
    return b, eps, R
 
def rhs(I, M, b, c, eps, R):
    dI = 1 - I - b * (1 + c * M) * I * M
    dM = eps * M * (R * (1 + c * M) / (1 + M) * I - 1)
    return dI, dM
 
def fixed_points(b, c, R):
    """(1,0) + positive roots of the cubic equation for M*:
       (1+M)(1+bM(1+cM)) = R(1+cM)  -> bc M^3 + b(c+1) M^2 + (b+1-Rc) M + (1-R) = 0"""
    pts = [(1.0, 0.0)]
    for M in np.roots([b * c, b * (c + 1), b + 1 - R * c, 1 - R]):
        if abs(M.imag) < 1e-9 and M.real > 1e-9:
            M = M.real
            pts.append((1 / (1 + b * M * (1 + c * M)), M))   # I from the I-nullcline
    return pts
 
def classify(I, M, b, c, eps, R):
    # Special case: the empty market (1,0) at R = 1. One eigenvalue is zero
    # (det J = 0), so the linearization says nothing. Close to (1,0),
    # dM/dt ~ eps (c - (1+b)) M^2, so the sign of c - (1+b) decides.
    if M == 0 and abs(R - 1) < 1e-9:
        return "stable node" if c < 1 + b else "unstable node"
    h, dh = (1 + c * M) / (1 + M), (c - 1) / (1 + M)**2
    J = np.array([[-1 - b * (1 + c * M) * M, -b * I * (1 + 2 * c * M)],
                  [eps * M * R * h,           eps * (R * h * I - 1) + eps * M * R * dh * I]])
    tr, det = np.trace(J), np.linalg.det(J)
    if det < 0: return "saddle"
    return ("stable" if tr < 0 else "unstable") + (" spiral" if tr**2 < 4 * det else " node")
 
def trdet(I, M, b, c, eps, R):
    """Trace and determinant of the Jacobian (also works with arrays)."""
    h, dh = (1 + c*M) / (1 + M), (c - 1) / (1 + M)**2
    a  = -1 - b*(1 + c*M)*M
    bb = -b*I*(1 + 2*c*M)
    cc = eps*M*R*h
    d  = eps*(R*h*I - 1) + eps*M*R*dh*I
    return a + d, a*d - bb*cc
 
# ---------- Drawing functions (used by the sliders and by --save) ----------
def draw_bifurcation(ax, b, eps, R, Mmax, cmax=CMAX_REPORT):
    """Bifurcation diagram of M against c for fixed R (no title: R goes in the caption)."""
    ax.clear()
    M = np.linspace(1e-3, Mmax, 4000)
    with np.errstate(all="ignore"):
        c = ((1 + M)*(1 + b*M) - R) / (M*(R - b*M*(1 + M)))
    ok = (c >= 0) & (c <= cmax)
    M, c = M[ok], c[ok]
    I = 1 / (1 + b*M*(1 + c*M))
    tr, det = trdet(I, M, b, c, eps, R)
    st = (det > 0) & (tr < 0)
    ax.scatter(c[st],  M[st],  s=3, color="green", label="stable")
    ax.scatter(c[~st], M[~st], s=3, color="black", label="unstable")
 
    # --- line M=0 (empty market) ---
    tol = 1e-9
    if R < 1 - tol:                      # eigenvalue eps(R-1) < 0: stable for all c
        ax.plot([0, cmax], [0, 0], color="green", lw=2, label=r"$\tilde M=0$ stable")
    elif R > 1 + tol:                    # eigenvalue > 0: saddle for all c
        ax.plot([0, cmax], [0, 0], color="black", ls="--", lw=2, label=r"$\tilde M=0$ unstable")
    else:                                # R = 1: second-order term eps(c-1-b)M^2 decides
        cb = b + 1
        ax.plot([0, min(cb, cmax)], [0, 0], color="green", lw=2,
                label=r"$\tilde M=0$ attracting ($c<1+b$)")
        ax.plot([min(cb, cmax), cmax], [0, 0], color="black", ls="--", lw=2,
                label=r"$\tilde M=0$ repelling ($c>1+b$)")
        ax.axvline(cb, color="orange", ls=":", lw=2, label=r"$c=1+b$")
 
    ax.set(xlim=(0, cmax), ylim=(-0.1, Mmax), xlabel=r"$c$", ylabel=r"$\tilde M$")
    ax.legend(markerscale=4, loc="upper right", fontsize=12)
 
def draw_phase(ax, b, c, eps, R, Mmax, Imax=IMAX_REPORT):
    """Phase portrait in the (M, I) plane. Returns the table rows with the equilibria."""
    ax.clear()
 
    # arrows (axes reversed: M horizontal, I vertical)
    Ig, Mg = np.meshgrid(np.linspace(0, Imax, 22), np.linspace(0, Mmax, 22))
    dI, dM = rhs(Ig, Mg, b, c, eps, R)
    k = 0.035 / (np.hypot(dI / Imax, dM / Mmax) + 1e-12)
    ax.quiver(Mg, Ig, dM * k, dI * k, angles="xy", scale_units="xy", scale=1,
              pivot="mid", color="gray", width=0.006)
 
    # nullclines (axes reversed)
    M = np.linspace(0, Mmax, 400)
    ax.plot(M, 1 / (1 + b * M * (1 + c * M)), "b", lw=2, label=r"$\dot{\tilde I}=0$")
    ax.plot(M, (1 + M) / (R * (1 + c * M)), "r", lw=2, label=r"$\dot{\tilde M}=0$")
    # thicker and not clipped, because it lies exactly on the left axis
    ax.axvline(0, color="r", lw=4, ls="--", clip_on=False, zorder=4,
               label=r"$\dot{\tilde M}=0$ ($\tilde M=0$)")
 
    # equilibria
    rows = []
    for (Is, Ms) in fixed_points(b, c, R):
        typ = classify(Is, Ms, b, c, eps, R)
        stab = typ.startswith("stable")
        ax.plot(Ms, Is, "o", ms=11, mfc="green" if stab else "white",
                mec="green" if stab else "black", mew=2, zorder=5, clip_on=False)
        rows.append([f"{Is:.3f}", f"{Ms:.3f}", typ])
 
    ax.set(xlim=(0, Mmax), ylim=(0, Imax), xlabel=r"$\tilde M$", ylabel=r"$\tilde I$")
    ax.legend(loc="upper right", fontsize=14)
    return rows
 
# ---------- Report figures ----------
def save_report_figures():
    """Save the 12 panels of Figure 1 as PDF, with b = 0.1 and eps = 0.01."""
    FIG_DIR.mkdir(exist_ok=True)
    for R, r_tag in R_VALUES.items():
        fig_b, ax_b = plt.subplots(figsize=(7, 5))
        draw_bifurcation(ax_b, B_REPORT, EPS_REPORT, R, MMAX_REPORT)
        fig_b.savefig(FIG_DIR / f"bif{r_tag}.pdf", bbox_inches="tight")
        plt.close(fig_b)
        for c, c_tag in C_VALUES.items():
            fig_p, ax_ph = plt.subplots(figsize=(8, 7.8))
            draw_phase(ax_ph, B_REPORT, c, EPS_REPORT, R, MMAX_REPORT)
            fig_p.savefig(FIG_DIR / f"phaseR{r_tag}c{c_tag}.pdf", bbox_inches="tight")
            plt.close(fig_p)
    print("Saved 12 figures to", FIG_DIR)
 
if "--save" in sys.argv:
    save_report_figures()
    sys.exit()
 
# ---------- Interactive figure ----------
fig2, ax2 = plt.subplots(figsize=(7, 5))
 
fig = plt.figure(figsize=(15, 8))
ax    = fig.add_axes([0.06, 0.08, 0.48, 0.87])
ax_p  = fig.add_axes([0.58, 0.70, 0.40, 0.26]); ax_p.axis("off")
ax_fp = fig.add_axes([0.58, 0.45, 0.40, 0.22]); ax_fp.axis("off")
 
sliders = []
for i, (n, lo, hi, v0) in enumerate(PARS):
    s = Slider(fig.add_axes([0.65, 0.5 - i * 0.04, 0.25, 0.025]), n, lo, hi, valinit=v0)
    s.label.set_fontsize(10)        # parameter name
    s.valtext.set_fontsize(10)      # current value
    sliders.append(s)
 
def show_table(a, rows, cols):
    a.clear(); a.axis("off")
    t = a.table(cellText=rows, colLabels=cols, loc="upper center", cellLoc="center")
    t.auto_set_font_size(False); t.set_fontsize(10); t.scale(1, 1.5)
 
def update(_=None):
    v = [s.val for s in sliders]
    I0, mu, beta0, c, s0, alpha, gamma, delta, Mmax = v
    b, eps, R = nondim(I0, mu, beta0, c, s0, alpha, gamma, delta)
 
    rows = draw_phase(ax, b, c, eps, R, Mmax)
    draw_bifurcation(ax2, b, eps, R, Mmax)
    ax2.set_title(f"$R={R:.3g}$")            # title only in the interactive window
    fig2.canvas.draw_idle()
 
    show_table(ax_p, [[r"$b$", f"{b:.4g}", r"$\beta_0/(\mu\alpha)$"],
                      [r"$c$", f"{c:.4g}", "c"],
                      [r"$\varepsilon$", f"{eps:.4g}", r"$\delta/\mu$"],
                      [r"$R$", f"{R:.4g}", r"$\gamma s_0\beta_0 I_0/(\mu\delta)$"]],
               ["parameter", "value", "definition"])
    show_table(ax_fp, rows, ["I*", "M*", "type"])
    fig.canvas.draw_idle()
 
for s in sliders: s.on_changed(update)
update()
 
plt.show()