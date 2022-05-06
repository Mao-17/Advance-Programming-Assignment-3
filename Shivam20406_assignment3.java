import java.util.Random;
import java.util.Scanner;
import java.lang.Math;
import java.util.ArrayList;
import java.util.Stack;
import java.util.Arrays;

class Matrix {

    int nRows;
    int nCols;
    double[][] mat;
    String name;

    Matrix(int rows, int cols, String name) {
        this.nRows = rows;
        this.nCols = cols;
        this.mat = new double[rows][cols];
        this.name = name;
    }

    Matrix(double[][] m) {

        this.nRows = m.length;
        this.nCols = m[0].length;
        this.mat = m;
    }

    Matrix(double[] m) {
        this.nRows = 1;
        this.nCols = m.length;

        this.mat = new double[1][m.length];
        this.mat[0] = m;

    }

    // Copy Constructor
    Matrix(Matrix m) {
        this.mat = deepCopy(m.mat);
        this.nRows = m.nRows;
        this.nCols = m.nCols;
    }

    public String getName() {
        return this.name;
    }

    public double[][] deepCopy(double[][] matrix) {
        return java.util.Arrays.stream(matrix).map(el -> el.clone()).toArray($ -> matrix.clone());
    }

    public double[][] toArray() {
        return this.mat;
    }

    public int[] shape() {
        int[] s = { this.nRows, this.nCols };
        return s;
    }

    public static boolean equals(Matrix m1, Matrix m2) {
        return Arrays.deepEquals(m1.mat, m2.mat);
    }

    public void printMatrix() {
        for (int i = 0; i < nRows; i++) {
            System.out.print("[ ");
            for (int j = 0; j < this.nCols; j++) {
                System.out.print(this.mat[i][j] + " ");
            }
            System.out.println(" ]");
        }
        System.out.println();
    }

    public static Matrix ones(Matrix m) {

        for (int i = 0; i < m.nRows; i++)
            for (int j = 0; j < m.nCols; j++)
                m.mat[i][j] = (double) 1;

        return m;
    }

    public static Matrix ones(int rows, int cols) {

        Matrix m = new Matrix(rows, cols);

        for (int i = 0; i < m.nRows; i++)
            for (int j = 0; j < m.nCols; j++)
                m.mat[i][j] = (double) 1;

        return m;
    }

    public static Matrix eye(int n) {
        Matrix m = new Matrix(n, n);

        for (int i = 0; i < m.nRows; i++)
            m.mat[i][i] = (double) 1;

        return m;
    }

    // Returns Matrix with randomly initialized integers between min and max
    public static Matrix rand(int rows, int cols, int min, int max) {
        Matrix m = new Matrix(rows, cols);

        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                m.mat[i][j] = (int) (Math.random() * (max - min)) + min;

        return m;
    }

}

public class Shivam20406_assignment3 {

    public static ArrayList<Matrix> matrixArray = new ArrayList<Matrix>();
    public static Scanner sc = new Scanner(System.in);

    public static void menu() {
        System.out.println("1. Take matrices as input and label them with appropriate matrix-types.");
        System.out
                .println("2. Create matrices of requested matrix-types and label them with appropriate matrix-types.");
        System.out.println("3. Change the elements of a matrix as long as the fixed matrix-type labels remain valid.");
        System.out.println("4. Display all the matrix-type labels of a requested matrix.");
        System.out.println("5. Perform addition, subtraction, multiplication & division.");
        System.out.println("6. Perform element-wise operations.");
        System.out.println("7. Transpose matrices.");
        System.out.println("8. Inverse matrices.");
        System.out.println("9. Compute means: row-wise mean, column-wise mean, mean of all the elements.");
        System.out.println("10. Compute determinants.");
        System.out.println("11. Use singleton matrices as scalars, if requested.");
        System.out.println("12. Compute A+AT for a matrix A.");
        System.out.println("13. Compute Eigen vectors and values.");
        System.out.println("14. Solve sets of linear equations using matrices.");
        System.out.println(
                "15. Retrieve all the existing matrices (entered or created) having requested matrix-type labels.");
        System.out.println("16. Exit");
    }

    public static Matrix getMatrixByName() {
        System.out.println("Enter name of the matrix: ");
        String name2 = sc.next();
        for (int j = 0; j < matrixArray.size(); j++) {
            Matrix matrix = matrixArray.get(j);
            String name = matrix.getName();
            if (name.compareTo(name2) == 0) {
                {
                    return matrix;
                }
            }
        }
        return null;
    }

    public static void main(String args[]) {
        int flag = 1;
        int start;
        while (flag == 1) {
            menu();
            start = sc.nextInt();
            if (start == 1) {
                System.out.println("Please enter the name of the matrix : ");
                String name_mat = sc.next();
                System.out.println("Please enter number of matrix rows : ");
                int row = sc.nextInt();
                System.out.println("Please enter number of matrix columns : ");
                int col = sc.nextInt();
                // defining two dimensional array java
                Matrix m = new Matrix(row, col, name_mat);
                // filling java matrix
                fillingMatrix(sc, m, row, col);
                // printing 2d array in matrix form in java
                m.printMatrix();
                matrixArray.add(m);
            } else if (start == 2) {

            } else if (start == 3) {

            } else if (start == 4) {
                Matrix matrix1 = getMatrixByName();
                if (isSquare(matrix1))
                    System.out.println("Square Matrix");
                if (isRectangular(matrix1))
                    System.out.println("Rectangular Matrix");
                if (isRow(matrix1))
                    System.out.println("Row Matrix");
                if (isColumn(matrix1))
                    System.out.println("Column Matrix");
                if (isSymmetric(matrix1))
                    System.out.println("Symmetric Matrix");
                if (isSkewSymmetric(matrix1))
                    System.out.println("Skew-symmetric Matrix");
                if (isUpperTriangular(matrix1))
                    System.out.println("Upper Triangular Matrix");
                if (isLowerTriangular(matrix1))
                    System.out.println("Lower Triangular Matrix");
                if (isSingleton(matrix1))
                    System.out.println("Singleton Matrix");
                if (isSingular(matrix1))
                    System.out.println("Singular Matrix");
                if (isDiagonal(matrix1))
                    System.out.println("Diagonal Matrix");
                if (isScalar(matrix1))
                    System.out.println("Scalar Matrix");
                if (isIdentity(matrix1))
                    System.out.println("Identity Matrix");
                if (isOnes(matrix1))
                    System.out.println("Ones Matrix");
                if (isNull(matrix1))
                    System.out.println("Null Matrix");
            } else if (start == 5) {
                System.out.println("Which operation would you like to perform (add/sub/mul): ");
                String opt = sc.next();
                System.out.println("\n" + opt);
                if (opt.compareTo("add") == 0) {
                    Matrix matrix1 = getMatrixByName();
                    Matrix matrix2 = getMatrixByName();
                    Matrix matrix3 = null;
                    if (matrix1 != null && matrix2 != null)
                        matrix3 = add(matrix1, matrix2);
                    if (matrix3 != null)
                        matrix3.printMatrix();
                } else if (opt.compareTo("sub") == 0) {
                    Matrix matrix1 = getMatrixByName();
                    Matrix matrix2 = getMatrixByName();
                    Matrix matrix3 = null;
                    if (matrix1 != null && matrix2 != null)
                        matrix3 = sub(matrix1, matrix2);
                    if (matrix3 != null)
                        matrix3.printMatrix();
                }
                if (opt.compareTo("mul") == 0) {
                    Matrix matrix1 = getMatrixByName();
                    Matrix matrix2 = getMatrixByName();
                    Matrix matrix3 = null;
                    if (matrix1 != null && matrix2 != null)
                        matrix3 = mul(matrix1, matrix2);
                    if (matrix3 != null)
                        matrix3.printMatrix();
                }
            } else if (start == 6) {
                System.out.println("Which operation would you like to perform (add/sub/mul/div): ");
                String opt = sc.next();
                if (opt.compareTo("add") == 0) {
                    Matrix matrix1 = getMatrixByName();
                    System.out.println("Please enter the element: ");
                    double ele = sc.nextDouble();
                    Matrix matrix2 = null;
                    if (matrix1 != null)
                        matrix2 = add(ele, matrix1);
                    if (matrix2 != null)
                        matrix2.printMatrix();
                }
                if (opt.compareTo("sub") == 0) {
                    Matrix matrix1 = getMatrixByName();
                    System.out.println("Please enter the element: ");
                    double ele = sc.nextDouble();
                    Matrix matrix2 = null;
                    if (matrix1 != null)
                        matrix2 = sub(ele, matrix1);
                    if (matrix2 != null)
                        matrix2.printMatrix();
                }
                if (opt.compareTo("mul") == 0) {
                    Matrix matrix1 = getMatrixByName();
                    System.out.println("Please enter the element: ");
                    double ele = sc.nextDouble();
                    Matrix matrix2 = null;
                    if (matrix1 != null)
                        matrix2 = mul(ele, matrix1);
                    if (matrix2 != null)
                        matrix2.printMatrix();
                }
                if (opt.compareTo("div") == 0) {
                    Matrix matrix1 = getMatrixByName();
                    System.out.println("Please enter the element: ");
                    double ele = sc.nextDouble();
                    Matrix matrix2 = null;
                    if (matrix1 != null)
                        matrix2 = div(ele, matrix1);
                    if (matrix2 != null)
                        matrix2.printMatrix();
                }

            } else if (start == 7) {

                Matrix matrix1 = getMatrixByName();
                Matrix matrix2 = null;
                if (matrix1 != null)
                    matrix2 = transpose(matrix1);
                if (matrix2 != null)
                    matrix2.printMatrix();
            } else if (start == 8) {
                Matrix matrix1 = getMatrixByName();
                Matrix matrix2 = null;
                if (matrix1 != null)
                    matrix2 = inv(matrix1);
                if (matrix2 != null)
                    matrix2.printMatrix();
            } else if (start == 9) {
                Matrix matrix1 = getMatrixByName();
                System.out.println("Which type of mean would you like to calculate(row/column/all): ");
                String opt = sc.next();
                if (opt.compareTo("row") == 0) {
                    if (matrix1 != null)
                        findrowMean(matrix1);
                }
                if (opt.compareTo("column") == 0) {
                    if (matrix1 != null)
                        findcolMean(matrix1);
                }
                if (opt.compareTo("all") == 0) {
                    if (matrix1 != null)
                        System.out.println(findMean(matrix1));
                }
            } else if (start == 10) {
                Matrix matrix1 = getMatrixByName();
                double deter = 0;
                if (matrix1 != null)
                    deter = matrixDeterminant(matrix1);
                System.out.println(deter);
            } else if (start == 11) {
                Matrix matrix1;
                double deter = 0;
                for (int i = 0; i < matrixArray.size(); i++) {
                    matrix1 = matrixArray.get(i);
                    if (isSingleton(matrix1)) {
                        deter = matrixDeterminant(matrix1);
                        System.out.println(deter);
                    }
                }
            } else if (start == 12) {

                Matrix matrix1 = getMatrixByName();
                Matrix matrix2 = null;
                Matrix matrix3 = null;
                if (matrix1 != null)
                    matrix2 = transpose(matrix1);
                if (matrix2 != null)
                    matrix3 = add(matrix1, matrix2);
                if (matrix3 != null)
                    matrix3.printMatrix();

            } else if (start == 13) {

            } else if (start == 14) {
                Matrix matrix1 = getMatrixByName();
                Matrix matrix2 = getMatrixByName();
                Matrix matrix3 = null;
                if (matrix1 != null && matrix2 != null)
                    matrix3 = mul(inv(matrix1), matrix2);
                if (matrix3 != null)
                    matrix3.printMatrix();
            } else if (start == 15) {
                System.out.println(
                        "Which type of matrices would you like to print(rec/square/column/row/symm/skews/upptri/lowtri/sing/diag/scalar/identity/singleton/ones/null): ");
                String opt = sc.next();
                Matrix matrix1;
                if (opt.compareTo("rec") == 0) {
                    for (int i = 0; i < matrixArray.size(); i++) {
                        matrix1 = matrixArray.get(i);
                        if (isRectangular(matrix1)) {
                            System.out.println(matrix1.getName());
                            matrix1.printMatrix();
                        }
                    }
                } else if (opt.compareTo("square") == 0) {
                    for (int i = 0; i < matrixArray.size(); i++) {
                        matrix1 = matrixArray.get(i);
                        if (isSquare(matrix1)) {
                            System.out.println(matrix1.getName());
                            matrix1.printMatrix();
                        }
                    }
                } else if (opt.compareTo("column") == 0) {
                    for (int i = 0; i < matrixArray.size(); i++) {
                        matrix1 = matrixArray.get(i);
                        if (isColumn(matrix1)) {
                            System.out.println(matrix1.getName());
                            matrix1.printMatrix();
                        }
                    }
                } else if (opt.compareTo("row") == 0) {
                    for (int i = 0; i < matrixArray.size(); i++) {
                        matrix1 = matrixArray.get(i);
                        if (isRow(matrix1)) {
                            System.out.println(matrix1.getName());
                            matrix1.printMatrix();
                        }
                    }
                } else if (opt.compareTo("symm") == 0) {
                    for (int i = 0; i < matrixArray.size(); i++) {
                        matrix1 = matrixArray.get(i);
                        if (isSymmetric(matrix1)) {
                            System.out.println(matrix1.getName());
                            matrix1.printMatrix();
                        }
                    }
                } else if (opt.compareTo("skews") == 0) {
                    for (int i = 0; i < matrixArray.size(); i++) {
                        matrix1 = matrixArray.get(i);
                        if (isSkewSymmetric(matrix1)) {
                            System.out.println(matrix1.getName());
                            matrix1.printMatrix();
                        }
                    }
                } else if (opt.compareTo("upptri") == 0) {
                    for (int i = 0; i < matrixArray.size(); i++) {
                        matrix1 = matrixArray.get(i);
                        if (isUpperTriangular(matrix1)) {
                            System.out.println(matrix1.getName());
                            matrix1.printMatrix();
                        }
                    }
                } else if (opt.compareTo("lowtri") == 0) {
                    for (int i = 0; i < matrixArray.size(); i++) {
                        matrix1 = matrixArray.get(i);
                        if (isLowerTriangular(matrix1)) {
                            System.out.println(matrix1.getName());
                            matrix1.printMatrix();
                        }
                    }
                } else if (opt.compareTo("sing") == 0) {
                    for (int i = 0; i < matrixArray.size(); i++) {
                        matrix1 = matrixArray.get(i);
                        if (isSingular(matrix1)) {
                            System.out.println(matrix1.getName());
                            matrix1.printMatrix();
                        }
                    }
                } else if (opt.compareTo("diag") == 0) {
                    for (int i = 0; i < matrixArray.size(); i++) {
                        matrix1 = matrixArray.get(i);
                        if (isDiagonal(matrix1)) {
                            System.out.println(matrix1.getName());
                            matrix1.printMatrix();
                        }
                    }
                } else if (opt.compareTo("scalar") == 0) {
                    for (int i = 0; i < matrixArray.size(); i++) {
                        matrix1 = matrixArray.get(i);
                        if (isScalar(matrix1)) {
                            System.out.println(matrix1.getName());
                            matrix1.printMatrix();
                        }
                    }
                } else if (opt.compareTo("identity") == 0) {
                    for (int i = 0; i < matrixArray.size(); i++) {
                        matrix1 = matrixArray.get(i);
                        if (isIdentity(matrix1)) {
                            System.out.println(matrix1.getName());
                            matrix1.printMatrix();
                        }
                    }
                } else if (opt.compareTo("singleton") == 0) {
                    for (int i = 0; i < matrixArray.size(); i++) {
                        matrix1 = matrixArray.get(i);
                        if (isSingleton(matrix1)) {
                            System.out.println(matrix1.getName());
                            matrix1.printMatrix();
                        }
                    }
                } else if (opt.compareTo("ones") == 0) {
                    for (int i = 0; i < matrixArray.size(); i++) {
                        matrix1 = matrixArray.get(i);
                        if (isOnes(matrix1)) {
                            System.out.println(matrix1.getName());
                            matrix1.printMatrix();
                        }
                    }
                } else if (opt.compareTo("null") == 0) {
                    for (int i = 0; i < matrixArray.size(); i++) {
                        matrix1 = matrixArray.get(i);
                        if (isNull(matrix1)) {
                            System.out.println(matrix1.getName());
                            matrix1.printMatrix();
                        }
                    }
                }
            } else if (start == 16) {
                flag = 0;
            }
        }
    }

    public static Matrix diag(Matrix m) {

        if (m.nRows != m.nCols)
            System.out.println("Incompatibe");

        Matrix diag_m = new Matrix(1, m.nRows);

        for (int i = 0; i < m.nRows; i++)
            diag_m.mat[0][i] = m.mat[i][i];

        return diag_m;
    }

    public static Matrix getColumn(Matrix m, int index) {
        return new Matrix(Matrix.getColumn(m.mat, index));
    }

    public static double[] getColumn(double[][] array, int index) {
        double[] column = new double[array[0].length];

        for (int i = 0; i < column.length; i++)
            column[i] = array[i][index];

        return column;
    }

    public static Matrix getColVec(Matrix m, int index) {
        Matrix ret = new Matrix(m.mat[0].length, 1);

        for (int i = 0; i < ret.nRows; i++)
            ret.mat[i][0] = m.mat[i][index];

        return ret;
    }

    public static double matrixDeterminant(Matrix m) {
        Matrix temporary;
        double result = 0;

        if (m.nCols == m.nRows) {
            if (m.nCols == 1) {
                result = m.mat[0][0];
                return result;
            }
            if (m.nCols == 2) {
                result = ((m.mat[0][0] * m.mat[1][1]) - (m.mat[0][1] * m.mat[1][0]));
                return result;
            }
            for (int i = 0; i < m.nRows; i++) {
                temporary = new Matrix(m.nRows - 1, m.nCols - 1, m.name);
                for (int j = 1; j < m.nCols; j++) {
                    for (int k = 0; k < m.nCols; k++) {
                        if (k < i) {
                            temporary.mat[j - 1][k] = m.mat[j][k];
                        } else if (k > i) {
                            temporary.mat[j - 1][k - 1] = m.mat[j][k];
                        }
                    }
                }
                result += m.mat[0][i] * Math.pow(-1, i) * matrixDeterminant(temporary);
            }
        } else {
            System.out.println("The matrix must be a square!");
        }
        return result;
    }

    public static Matrix transpose(Matrix m) {

        Matrix tran = new Matrix(m.nCols, m.nRows, m.name);

        for (int i = 0; i < m.nRows; i++) {

            for (int j = 0; j < m.nCols; j++) {

                tran.mat[j][i] = m.mat[i][j];

            }
        }

        return tran;
    }

    static boolean isSquare(Matrix m) {
        int N = m.nRows;
        int M = m.nCols;
        if (M == N)
            return true;
        else
            return false;
    }

    static boolean isRectangular(Matrix m) {
        int N = m.nRows;
        int M = m.nCols;
        if (M != N)
            return true;
        else
            return false;
    }

    static boolean isRow(Matrix m) {
        int N = m.nRows;
        if (N == 1)
            return true;
        else
            return false;
    }

    static boolean isColumn(Matrix m) {
        int M = m.nCols;
        if (M == 1)
            return true;
        else
            return false;
    }

    static boolean isSymmetric(Matrix m) {
        int N = m.nRows;
        int M = m.nCols;
        Matrix tr = transpose(m);
        if (N == M) {
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    if (m.mat[i][j] != tr.mat[i][j])
                        return false;
            return true;
        } else {
            return false;
        }
    }

    static boolean isSkewSymmetric(Matrix m) {
        int N = m.nRows;
        int M = m.nCols;
        Matrix tr = transpose(m);
        if (N == M) {
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    if (m.mat[i][j] != -tr.mat[i][j])
                        return false;
            return true;
        } else {
            return false;
        }
    }

    static boolean isLowerTriangular(Matrix m) {
        int N = m.nRows;
        int M = m.nCols;
        if (N == M) {
            for (int i = 0; i < N; i++)
                for (int j = i + 1; j < N; j++)
                    if (m.mat[i][j] != 0)
                        return false;
            return true;
        } else {
            return false;
        }
    }

    static boolean isUpperTriangular(Matrix m) {
        int N = m.nRows;
        int M = m.nCols;
        if (N == M) {
            for (int i = 1; i < N; i++)
                for (int j = 0; j < i; j++)
                    if (m.mat[i][j] != 0)
                        return false;
            return true;
        } else {
            return false;
        }
    }

    static boolean isSingular(Matrix m) {
        int N = m.nRows;
        int M = m.nCols;
        if (N == M) {
            if (matrixDeterminant(m) == 0) {
                return true;
            } else
                return false;
        } else {
            return false;
        }
    }

    static boolean isDiagonal(Matrix m) {
        int N = m.nRows;
        int M = m.nCols;
        if (N == M) {
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    if ((i != j) && (m.mat[i][j] != 0))
                        return false;
            return true;
        } else {
            return false;
        }
    }

    static boolean isScalar(Matrix m) {
        int N = m.nRows;
        int M = m.nCols;
        if (N == M) {
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    if ((i != j) && (m.mat[i][j] != 0))
                        return false;
            for (int i = 0; i < N - 1; i++)
                if (m.mat[i][i] != m.mat[i + 1][i + 1])
                    return false;
            return true;
        } else {
            return false;
        }
    }

    static boolean isIdentity(Matrix m) {
        if (isScalar(m)) {
            if (m.mat[0][0] == 1)
                return true;
            else
                return false;
        } else {
            return false;
        }
    }

    static boolean isSingleton(Matrix m) {
        int N = m.nRows;
        int M = m.nCols;
        if (N == M) {
            if (N == 1)
                return true;
            else
                return false;
        } else {
            return false;
        }
    }

    static boolean isOnes(Matrix m) {
        int N = m.nRows;
        int M = m.nCols;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                if (m.mat[i][j] != 1)
                    return false;
            }
        }
        return true;
    }

    static boolean isNull(Matrix m) {
        int N = m.nRows;
        int M = m.nCols;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                if (m.mat[i][j] != 0)
                    return false;
            }
        }
        return true;
    }

    public static Matrix mul(double a, Matrix m) {

        Matrix ret = new Matrix(m);

        for (int i = 0; i < m.nRows; i++)
            for (int j = 0; j < m.nCols; j++)
                ret.mat[i][j] *= a;

        return ret;
    }

    // Element-wise divison to m
    public static Matrix div(double a, Matrix m) {

        Matrix ret = new Matrix(m);

        for (int i = 0; i < m.nRows; i++)
            for (int j = 0; j < m.nCols; j++)
                ret.mat[i][j] /= a;

        return ret;
    }

    // Element-wise addition to m
    public static Matrix add(double a, Matrix m) {

        Matrix ret = new Matrix(m);

        for (int i = 0; i < m.nRows; i++)
            for (int j = 0; j < m.nCols; j++)
                ret.mat[i][j] += a;

        return ret;
    }

    // Element-wise subtraction to m
    public static Matrix sub(double a, Matrix m) {

        Matrix ret = new Matrix(m);

        for (int i = 0; i < m.nRows; i++)
            for (int j = 0; j < m.nCols; j++)
                ret.mat[i][j] -= a;

        return ret;
    }

    public static Matrix add(Matrix... args) {

        Matrix res = new Matrix(args[0]);

        for (int m = 1; m < args.length; m++) {

            if (res.nRows != args[m].nRows || res.nCols != args[m].nCols) {
                // throw new IncompatibleDimensionsException("Incompatible Dimensions: (" +
                // res.nRows + ", " + res.nCols
                // + ") and (" + args[m].nRows + ", " + args[m].nCols + ")");
                System.out.println("Incompatibe");
                return null;
            }
            for (int i = 0; i < res.nRows; i++)
                for (int j = 0; j < res.nCols; j++)
                    res.mat[i][j] += args[m].mat[i][j];
        }

        return res;

    }

    // Subtract two or more matrices from first matrix argument.
    public static Matrix sub(Matrix... args) {

        Matrix res = new Matrix(args[0]);

        for (int m = 1; m < args.length; m++) {

            if (res.nRows != args[m].nRows || res.nCols != args[m].nCols)
                // throw new IncompatibleDimensionsException("Incompatible Dimensions: (" +
                // res.nRows + ", " + res.nCols
                // + ") and (" + args[m].nRows + ", " + args[m].nCols + ")");
                System.out.println("Incompatibe");

            for (int i = 0; i < res.nRows; i++)
                for (int j = 0; j < res.nCols; j++)
                    res.mat[i][j] -= args[m].mat[i][j];
        }

        return res;

    }

    public static Matrix mul(Matrix... args) {
        Matrix ret = args[0];

        for (int m = 1; m < args.length; m++) {

            if (ret.nCols != args[m].nRows)
                // throw new IncompatibleDimensionsException("Incompatible Dimensions: (" +
                // ret.nRows + ", " + ret.nCols
                // + ") and (" + args[m].nRows + ", " + args[m].nCols + ")");
                System.out.println("Incompatibe");

            Matrix res = new Matrix(ret.nRows, args[m].nCols, ret.name);

            for (int i = 0; i < ret.nRows; i++) {

                for (int j = 0; j < args[m].nCols; j++) {

                    for (int k = 0; k < args[m].nRows; k++) {

                        res.mat[i][j] += ret.mat[i][k] * args[m].mat[k][j];

                    }

                }

            }

            ret = res;

        }

        return ret;
    }

    public static double findMean(Matrix m) {
        int sum = 0;
        int N = m.nRows;
        int M = m.nCols;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < M; j++)
                sum += m.mat[i][j];

        return (double) sum / (N * M);
    }

    public static void findcolMean(Matrix m) {
        int sum = 0;
        int N = m.nRows;
        int M = m.nCols;
        for (int i = 0; i < M; i++) {
            sum = 0;
            for (int j = 0; j < N; j++) {
                sum += m.mat[j][i];

            }
            System.out.println((double) sum / (N));
        }
    }

    public static void findrowMean(Matrix m) {
        int sum = 0;
        int N = m.nRows;
        int M = m.nCols;
        for (int i = 0; i < N; i++) {
            sum = 0;
            for (int j = 0; j < M; j++) {
                sum += m.mat[i][j];

            }
            System.out.println((double) sum / (M));
        }

    }

    private static Matrix minor(Matrix m, int row, int col) {

        Matrix min = new Matrix(m.nRows - 1, m.nCols - 1, m.name);
        int k = 0, l = 0;

        for (int i = 0; i < m.nRows; i++) {

            if (i == row)
                continue;

            for (int j = 0; j < m.nCols; j++) {

                if (j == col)
                    continue;

                min.mat[k][l] = m.mat[i][j];

                l = ++l % min.nCols;

                if (l == 0)
                    k++; // l becomes 0 when we've filled all columns in a row.

            }

        }

        return min;

    }

    public static Matrix tranCof(Matrix m) {
        Matrix c = new Matrix(m.nRows, m.nCols, m.name);

        for (int i = 0; i < m.nRows; i++) {

            for (int j = 0; j < m.nCols; j++) {

                c.mat[j][i] = Math.pow(-1, i + j) * matrixDeterminant(minor(m, i, j));

            }
        }

        return c;
    }

    public static Matrix inv(Matrix m) {
        return mul(1 / matrixDeterminant(m), tranCof(m));
    }

    public static void fillingMatrix(Scanner scan, Matrix m, int rows, int columns) {
        System.out.println("Please enter elements in matrix : ");
        for (int a = 0; a < rows; a++) {
            for (int b = 0; b < columns; b++) {
                m.mat[a][b] = scan.nextInt();
            }
        }
    }
}