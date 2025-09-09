import numpy as np

def gaussian_elimination_3x4(matrix):
    """
    Performs Gaussian elimination on a 3x4 augmented matrix.
    The result is the matrix in row echelon form.
    """
    A = matrix.astype(float)

    # Make sure A[0,0] != 0
    if A[0, 0] == 0:
        for i in range(1, 3):
            if A[i, 0] != 0:
                A[[0, i]] = A[[i, 0]]
                break

    # Eliminate first column below first row
    for i in range(1, 3):
        factor = A[i, 0] / A[0, 0]
        A[i] = A[i] - factor * A[0]

    # Make sure A[1,1] != 0
    if A[1, 1] == 0:
        if A[2, 1] != 0:
            A[[1, 2]] = A[[2, 1]]

    # Eliminate second column below second row
    factor = A[2, 1] / A[1, 1] if A[1, 1] != 0 else 0
    A[2] = A[2] - factor * A[1]

    return A

def back_substitution(A):
    """
    Solve the upper triangular system for x, y, z.
    Assumes last column is the augmented values.
    """
    x = np.zeros(3)
    # z
    x[2] = A[2, 3] / A[2, 2]
    # y
    x[1] = (A[1, 3] - A[1, 2] * x[2]) / A[1, 1]
    # x
    x[0] = (A[0, 3] - A[0, 1] * x[1] - A[0, 2] * x[2]) / A[0, 0]
    return x

# Augmented matrix
matrix = np.array([
    [-2, 0, 1, -4],
    [-1, 7, 1, -50],
    [5, -1, 1, -26]
])

copy_matrix = matrix.copy()

print("Original Matrix:\n")
print(matrix)

print("\nMy solution:")
echelon = gaussian_elimination_3x4(matrix)
print("Row Echelon Form:")
print(echelon)
solution = back_substitution(echelon)
print("Solution [x, y, z]:")
print(solution)

print("\n\nCompare to NumPy's built-in solver:")
print(np.linalg.solve(copy_matrix[:, :3], copy_matrix[:, 3]))

