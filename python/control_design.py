"""Feedback design for the inverted cart-pole: LQR vs. pole placement.

Both controllers are full-state feedback laws u = -K z, where z is the
deviation from the upright equilibrium (theta = pi). Each design is checked
on the linearized model and then applied to the full nonlinear model to
confirm it actually holds the pendulum up from a perturbed start.

Outputs (written to ../figures/):
    closed_loop_linear.svg    angle deviation, LQR vs. pole placement
    closed_loop_nonlinear.svg angle under both controllers, nonlinear model

Running this file also prints the gains, closed-loop poles, and response
metrics (settling time, peak |u|) for both designs.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_continuous_are
from scipy.signal import place_poles
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from cartpole import dynamics, linearize

FIGDIR = Path(__file__).resolve().parent.parent / "figures"
FIGDIR.mkdir(exist_ok=True)

A, B = linearize(np.pi)          # inverted equilibrium
B = B.ravel()
Z0 = np.array([0.0, 0.0, -0.2, 0.0])   # start 0.2 rad off upright, as in the free-response study
T = 6.0


def lqr(A, B, Q, R):
    """State-feedback gain from the continuous algebraic Riccati equation."""
    P = solve_continuous_are(A, B.reshape(-1, 1), Q, R)
    return (np.linalg.solve(R, B.reshape(-1, 1).T @ P)).ravel()


def simulate_linear(K, z0=Z0, T=T):
    Acl = A - np.outer(B, K)
    sol = solve_ivp(lambda t, z: Acl @ z, (0, T), z0,
                    t_eval=np.linspace(0, T, 1201), rtol=1e-9, atol=1e-11)
    u = -sol.y.T @ K
    return sol.t, sol.y, u


def simulate_nonlinear(K, z0=Z0, T=T):
    eq = np.array([0.0, 0.0, np.pi, 0.0])
    X0 = eq + z0
    def rhs(t, X):
        return dynamics(X, u=float(-K @ (X - eq)))
    sol = solve_ivp(rhs, (0, T), X0, t_eval=np.linspace(0, T, 1201),
                    rtol=1e-9, atol=1e-11)
    z = sol.y - eq.reshape(-1, 1)
    u = -z.T @ K
    return sol.t, z, u


def settling_time(t, theta, band=0.02):
    """Last time |theta deviation| exceeds `band` rad (2% of a radian scale)."""
    outside = np.where(np.abs(theta) > band)[0]
    return t[outside[-1]] if len(outside) else 0.0


def report(name, K, t, z, u):
    poles = np.linalg.eigvals(A - np.outer(B, K))
    print(f"{name}")
    print("  K =", np.array2string(K, precision=3, suppress_small=True))
    print("  closed-loop poles:", ", ".join(f"{p:.3f}" for p in
          sorted(poles, key=lambda p: -p.real)))
    print(f"  settling time (|theta| < 0.02 rad): {settling_time(t, z[2]):.2f} s")
    print(f"  peak cart excursion: {np.max(np.abs(z[0])):.3f} m")
    print(f"  peak control force: {np.max(np.abs(u)):.1f} N\n")


if __name__ == "__main__":
    # --- LQR: penalize angle deviation hardest, then cart position ---
    Q = np.diag([10.0, 1.0, 100.0, 1.0])
    R = np.array([[1.0]])
    K_lqr = lqr(A, B, Q, R)

    # --- Pole placement: two slow real poles, one faster complex pair ---
    desired = [-1.5, -2.0, -3.0 + 1.0j, -3.0 - 1.0j]
    K_pp = place_poles(A, B.reshape(-1, 1), desired).gain_matrix.ravel()

    t, z_lqr, u_lqr = simulate_linear(K_lqr)
    _, z_pp, u_pp = simulate_linear(K_pp)
    report("LQR (Q = diag(10, 1, 100, 1), R = 1)", K_lqr, t, z_lqr, u_lqr)
    report("Pole placement at {-1.5, -2, -3 +/- j}", K_pp, t, z_pp, u_pp)

    plt.figure(figsize=(7.5, 4.2))
    plt.plot(t, np.degrees(z_lqr[2]), label="LQR")
    plt.plot(t, np.degrees(z_pp[2]), "--", label="Pole placement")
    plt.axhline(0, color="k", lw=0.6)
    plt.xlabel("Time (s)")
    plt.ylabel("Angle deviation from upright (deg)")
    plt.title("Closed-loop response, linearized model (0.2 rad initial offset)")
    plt.legend()
    plt.grid(alpha=0.35)
    plt.tight_layout()
    plt.savefig(FIGDIR / "closed_loop_linear.svg")
    plt.close()

    # --- Same gains on the nonlinear model ---
    tn, zn_lqr, _ = simulate_nonlinear(K_lqr)
    _, zn_pp, _ = simulate_nonlinear(K_pp)
    print("Nonlinear model, final angle deviation (deg): "
          f"LQR {np.degrees(zn_lqr[2, -1]):.4f}, "
          f"pole placement {np.degrees(zn_pp[2, -1]):.4f}")

    plt.figure(figsize=(7.5, 4.2))
    plt.plot(tn, np.degrees(zn_lqr[2]), label="LQR")
    plt.plot(tn, np.degrees(zn_pp[2]), "--", label="Pole placement")
    plt.axhline(0, color="k", lw=0.6)
    plt.xlabel("Time (s)")
    plt.ylabel("Angle deviation from upright (deg)")
    plt.title("Closed-loop response, nonlinear model (same gains)")
    plt.legend()
    plt.grid(alpha=0.35)
    plt.tight_layout()
    plt.savefig(FIGDIR / "closed_loop_nonlinear.svg")
    plt.close()
    print("Figures written to", FIGDIR)
