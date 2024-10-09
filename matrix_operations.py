import numpy as np

class MatrixOperations:

    def add_matrices(A, B):
        return np.add(A, B)

    def multiply_matrices(A,B):
        return np.dot(A,B)

    def determinant(A):
        return np.linalg.det(A)

    def transpose(A):
        return np.transpose(A)

    def eigenvector(A):
        values, vectors = np.linalg.eig(A)
        return values, vectors

   def power_iter(A, num_simulations:int):
       b_k = np.random.rand(A.shape[0])

        for _ in range(num_simulations):
            b_k1 = np.dot(A, b_k)
            b_k1norm = np.linalg.norm(b_k1)
            b_k = b_k1 / b_k1norm

        eigenvalue = np.dot(b_k.T, np.dot(A, b_k)) / np.dot(b_k.T, b_k)
        return eigenvalue, b_k
