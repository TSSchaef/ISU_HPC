import numpy as np
import matplotlib.pyplot as plt


def GaussElimination(A, b):
    """
    Solves Ax = b using Gaussian elimination for any size of A and b.
    Arguments:
        A : np.ndarray
            Coefficient matrix of shape (n, n)
        b : np.ndarray
            Right-hand side vector or matrix of shape (n,)
    Returns:
        x : np.ndarray
            Solution vector or matrix
    """
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    n = A.shape[0]

    # Construct the augmented matrix
    if b.ndim == 1:
        b = b.reshape(-1, 1)
    Ab = np.hstack([A, b])

    # Forward elimination
    for i in range(n):
        max_row = np.argmax(np.abs(Ab[i:, i])) + i
        if i != max_row:
            Ab[[i, max_row]] = Ab[[max_row, i]]

        # Eliminate below
        for j in range(i+1, n):
            factor = Ab[j, i] / Ab[i, i]
            Ab[j, i:] -= factor * Ab[i, i:]

    # Back substitution
    x = np.zeros((n, b.shape[1]))
    for k in range(b.shape[1]):
        for i in reversed(range(n)):
            x[i, k] = (Ab[i, n+k] - np.dot(Ab[i, i+1:n], x[i+1:, k])) / Ab[i, i]

    # Return as 1D for single rhs
    if x.shape[1] == 1:
        return x.flatten()
    return x


if __name__ == "__main__":
    # Matrix from last lab
    #matrix = np.array([
        #[-2, 0, 1, -4],
        #[-1, 7, 1, -50],
        #[5, -1, 1, -26]
    #])

    # Interpolation points
    xs = np.array([-0.1, -0.02, 0.02, 0.1])
    ys = np.cos(xs)

    # Construct Vandermonde matrix for cubic polynomial
    A = np.vstack([xs**3, xs**2, xs, np.ones_like(xs)]).T
    b = ys

    # Solve for coefficients [a, b, c, d]
    coeffs = GaussElimination(A, b)
    a, b_, c, d = coeffs

    # Polynomial p(x)
    def p(x):
        return a*x**3 + b_*x**2 + c*x + d

    # Plot
    x_plot = np.linspace(xs.min()-0.05, xs.max()+0.05, 200)
    plt.plot(x_plot, np.cos(x_plot), label='f(x)=cos(x)', lw=2)
    plt.plot(x_plot, p(x_plot), label='p(x) cubic interpolation', lw=2)
    plt.scatter(xs, ys, color='red', label='Interpolation points')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Cubic Polynomial Interpolation of cos(x)')
    plt.grid(True)
    plt.show()




