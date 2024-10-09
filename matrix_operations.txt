#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    int rws;
    int cols;
    double **data;
} Matrix;

Matrix create_matrix(int rows, int cols) {
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.data = (double**) malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        mat.data[i] = (double*) malloc(cols * sizeof(double));
    }
    return mat;
}

void free_matrix(Matrix *mat) {
    for (int i = 0; i < mat->rows; i++) {
        free(mat->data[i]);
    }
    free(mat->data);
}

void print_matrix(Matrix mat) {
    for (int i = 0; i < mat.rows; i++) {
        for (int j = 0; j < mat.cols; j++) {
            printf("%lf ", mat.data[i][j]);
        }
        printf("\n");
    }
}

Matrix add_matrices(Matrix A, Matrix B) {
    if (A.rows != B.rows || A.cols != B.cols) {
        printf("ERROR: Matrices dimensions do not match\n");
        exit(EXIT_FAILURE);
    }

    Matrix result = create_matrix(A.rows, A.cols);
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < A.cols; j++) {
            result.data[i][j] = A.data[i][j] + B.data[i][j];
        }
    }
    return result;
}

Matrix multiply_matrices(Matrix A, Matrix B) {
    if (A.cols != B.rows) {
        printf("ERROR: Matrices dimensions are not compatible\n");
        exit(EXIT_FAILURE);
    }

    Matrix result = create_matrix(A.rows, B.cols);
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < B.cols; j++) {
            result.data[i][j] = 0;
            for (int k = 0; k < A.cols; k++) {
                result.data[i][j] += A.data[i][k] * B.data[k][j];
            }
        }
    }
    return result;
}

double determinant(Matrix mat) {
    if (mat.rows != mat.cols) {
        printf("ERROR: Not square matrices\n");
        exit(EXIT_FAILURE);
    }

    if (mat.rows == 1) {
        return mat.data[0][0];
    }

    if (mat.rows == 2) {
        return mat.data[0][0] * mat.data[1][1] - mat.data[0][1] * mat.data[1][0];
    }

    double det = 0.0;
    for (int k = 0; k < mat.cols; k++) {
        Matrix submat = create_matrix(mat.rows - 1, mat.cols - 1);

        for (int i = 1; i < mat.rows; i++) {
            int subi = 0;
            for (int j = 0; j < mat.cols; j++) {
                if (j == k) continue;
                submat.data[i-1][subi] = mat.data[i][j];
                subi++;
            }
        }

        double sign = (k % 2 == 0) ? 1 : -1;
        det += sign * mat.data[0][k] * determinant(submat);

        free_matrix(&submat);
    }
    return det;
}

Matrix inverse_matrix(Matrix A) {
    if (A.rows != A.cols) {
        printf("ERROR: Not square matrices\n");
        exit(EXIT_FAILURE);
    }

    int n = A.rows;
    Matrix augmented = create_matrix(n, 2 * n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmented.data[i][j] = A.data[i][j];
        }
        augmented.data[i][n + i] = 1.0;
    }


    for (int i = 0; i < n; i++) {
        if (augmented.data[i][i] == 0) {
            printf("ERROR: Singular matrix\n");
            free_matrix(&augmented);
            exit(EXIT_FAILURE);
        }

    
        double pivot = augmented.data[i][i];
        for (int j = 0; j < 2 * n; j++) {
            augmented.data[i][j] /= pivot;
        }

        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = augmented.data[k][i];
                for (int j = 0; j < 2 * n; j++) {
                    augmented.data[k][j] -= factor * augmented.data[i][j];
                }
            }
        }
    }

    Matrix inverse = create_matrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inverse.data[i][j] = augmented.data[i][n + j];
        }
    }

    free_matrix(&augmented);
    return inverse;
}

void normalize_vector(double *vector, int size) {
    double norm = 0.0;
    for (int i = 0; i < size; i++) {
        norm += vector[i] * vector[i];
    }
    norm = sqrt(norm);
    for (int i = 0; i < size; i++) {
        vector[i] /= norm;
    }
}

void power_iteration(Matrix A, double *eigenvector, int num_iterations) {
    double *b_k = (double *) malloc(A.rows * sizeof(double));

    // Initialize eigenvector to a random non-zero vector
    for (int i = 0; i < A.rows; i++) {
        b_k[i] = 1.0;  // Can randomize if needed
    }

    for (int i = 0; i < num_iterations; i++) {
        // Multiply A by the vector b_k
        double *b_k1 = (double *) calloc(A.rows, sizeof(double));
        for (int row = 0; row < A.rows; row++) {
            for (int col = 0; col < A.cols; col++) {
                b_k1[row] += A.data[row][col] * b_k[col];
            }
        }
        
         normalize_vector(b_k1, A.rows);

         for (int j = 0; j < A.rows; j++) {
            b_k[j] = b_k1[j];
        }

        free(b_k1);
    }

    for (int i = 0; i < A.rows; i++) {
        eigenvector[i] = b_k[i];
    }

    free(b_k);
}

