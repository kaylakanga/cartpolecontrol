"""Cart-pole (inverted pendulum on a cart): model, linearization, and analysis.

The state is X = [x, v, theta, omega]:
    x      cart position (m)
    v      cart velocity (m/s)
    theta  pendulum angle from the downward (hanging) position (rad)
    omega  pendulum angular velocity (rad/s)

theta = 0 is the hanging equilibrium and theta = pi is the inverted
(upright) equilibrium. The only input is a horizontal force u on the cart.

Running this file prints the eigenvalues of both linearized models and the
controllability rank of each.
"""

import numpy as np

# Physical parameters
M = 4.0      # cart mass (kg)
m = 0.6      # pendulum mass (kg)
L = 1.2      # pendulum arm length (m)
d = 0.4      # viscous damping coefficient of the cart
g = 9.81     # gravitational acceleration (m/s^2)


def dynamics(X, u=0.0):
    """Nonlinear equations of motion, dX/dt = dynamics(X, u)."""
    _, v, theta, omega = X
    s, c = np.sin(theta), np.cos(theta)
    D = m * L**2 * (M + m * s**2)
    dx = np.empty(4)
    dx[0] = v
    dx[1] = (m**2 * L**2 * g * c * s
             + m * L**2 * (m * L * omega**2 * s - d * v)
             + m * L**2 * u) / D
    dx[2] = omega
    dx[3] = (-(m + M) * m * g * L * s
             - m * L * c * (m * L * omega**2 * s - d * v)
             - m * L * c * u) / D
    return dx


def linearize(theta_eq):
    """Analytic linearization (A, B) about an equilibrium angle theta_eq.

    Equilibria have v = omega = 0 and sin(theta_eq) = 0, so only the
    cos(theta_eq) = +/-1 terms survive in the Jacobian.
    """
    c = np.cos(theta_eq)
    A = np.array([
        [0.0, 1.0, 0.0,             0.0],
        [0.0, -d / M, m * g / M,       0.0],
        [0.0, 0.0, 0.0,             1.0],
        [0.0, c * d / (M * L), -c * (m + M) * g / (M * L), 0.0],
    ])
    B = np.array([[0.0], [1.0 / M], [0.0], [-c / (M * L)]])
    return A, B


def controllability_rank(A, B):
    """Rank of the controllability matrix [B, AB, A^2 B, A^3 B]."""
    n = A.shape[0]
    ctrb = np.hstack([np.linalg.matrix_power(A, k) @ B for k in range(n)])
    return np.linalg.matrix_rank(ctrb), ctrb


def describe(name, A):
    eigs = np.linalg.eigvals(A)
    order = np.argsort([-e.real for e in eigs])
    print(f"{name} equilibrium")
    print("  A =\n", A)
    print("  eigenvalues:", ", ".join(f"{e:.4f}" for e in eigs[order]))
    return eigs


if __name__ == "__main__":
    for name, theta_eq in [("Hanging (theta = 0)", 0.0),
                           ("Inverted (theta = pi)", np.pi)]:
        A, B = linearize(theta_eq)
        describe(name, A)
        print("  B =", B.ravel())
        rank, _ = controllability_rank(A, B)
        print(f"  controllability rank: {rank} of {A.shape[0]}\n")
