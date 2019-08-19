package Jama;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;

/**
 * A nice sparse matrix. We assume the dimensions are not changing.Use it
 * carefully!!!!
 *
 * Only the product and "set" are rewritten
 *
 * @author aeromero
 */
public class SparseMatrix extends Matrix {

    /**
     * For each row stores the indices of non-zero columns
     */
    ArrayList<HashSet<Integer>> nonZeroIds;

    public SparseMatrix(int m, int n) {
        super(m, n);
        if (m > 0) {
            nonZeroIds = new ArrayList<HashSet<Integer>>();
            for (int i = 0; i < m; ++i) {
                this.nonZeroIds.add(new HashSet<Integer>());
            }
        }
    }

    @Override
    public void set(int i, int j, float s) {
        super.set(i, j, s);
        if (s != 0.0f) {
            this.nonZeroIds.get(i).add(j);
        }
    }

    @Override
    public Matrix times(Matrix B) {
        if (B.m != n) {
            throw new IllegalArgumentException("Matrix inner dimensions must agree.");
        }
        final int n_b = B.n;
        Matrix X = new Matrix(m, n_b);

        if (nonZeroIds != null) {
            float[][] _C = X.getArray();
            final float[][] _A = this.getArray();
            final float[][] _B = B.getArray();

            for (int i = 0; i < m; i++) {
                for (int k : this.nonZeroIds.get(i)) {
                    for (int j = 0; j < n_b; ++j) {
                        _C[i][j] += _A[i][k] * _B[k][j];
                    }
                }
            }
        }
        return X;
    }

    public static Matrix getRandom(int n, int m) {
        SparseMatrix rand = new SparseMatrix(n, m);
        Random r = new Random();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; ++j) {
                rand.set(i, j, r.nextFloat());
            }
        }
        return rand;
    }

    public static void main(String args[]) {
        int size_m = 2000;
        int size_n = 3000;


        Matrix A = SparseMatrix.getRandom(size_m, size_n);
    }
}
