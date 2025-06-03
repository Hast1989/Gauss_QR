#include <iostream>
#include <fstream>
#include<cmath>
#include <iomanip>
int indikG = 0;
int indikQR = 0;
template<typename T>
T NormVec1(T* x, int n)
{
    T res = 0;
    for (int i = 0; i < n; i++)
    {
        res = res + abs(x[i]);
    }
    return res;
}
template<typename T>
T NormVecinf(T* x, int n)
{
    T max = 0;
    for (int i = 0; i < n; i++)
    {
        if (max < abs(x[i]))
        {
            max = abs(x[i]);
        }
    }
    return max;
}
template<typename T>
void MultM(T** Matrixl, T** Matrixr, T** Mresult, int n)
{
    T res;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            res = 0;
            for (int k = 0; k < n; k++)
            {
                res = res + Matrixl[i][k] * Matrixr[k][j];
            }
            Mresult[i][j] = res;
        }
    }
}
template<typename T>
T NormVect(T* x, int n)
{
    T res = 0;
    for (int i = 0; i < n; i++)
    {
        res = res + x[i] * x[i];
    }
    res = sqrt(res);
    return res;
}
template<typename T>
void MultWV(T** Matrix, T* x, T* result, int n)
{
    for (int i = 0; i < n; i++)
    {
        result[i] = 0;
        for (int j = 0; j < n; j++)
        {
            result[i] = result[i] + Matrix[i][j] * x[j];
        }
    }
}
template<typename T>
int FindLeader(T** Matrix, int n, int j)
{
    int maxint = j;
    T max = 0;
    for (int i = j; i < n; i++)
    {
        if (abs(Matrix[i][j]) > max)
        {
            max = abs(Matrix[i][j]);
            maxint = i;
        }
    }
    if (max == 0)
    {
        indikG = 1;
    }
    return maxint;
}
template<typename T>
void SwichLines(T** Matrix, T* rightb, int i, int j)
{
    T resb = rightb[i];
    T* resline = Matrix[i];
    Matrix[i] = Matrix[j];
    Matrix[j] = resline;
    rightb[i] = rightb[j];
    rightb[j] = resb;

}
template<typename T>
void Gauss(T** Matrix, T* rightb, int n, T* x)
{
    T d, s;
    for (int k = 0; k < n - 1; k++)
    {
        SwichLines(Matrix, rightb, FindLeader(Matrix, n, k), k);
        if (indikG == 1)
        {
            return;
        }
        for (int j = k + 1; j < n; j++)
        {
            if (abs(Matrix[j][k]) != 0)
            {
                d = Matrix[j][k] / Matrix[k][k];
                Matrix[j][k] = 0;
                for (int i = k + 1; i < n; i++)
                {
                    Matrix[j][i] = Matrix[j][i] - d * Matrix[k][i];
                }
                rightb[j] = rightb[j] - d * rightb[k];
            }
        }
    }
    if (fabs(Matrix[n - 1][n - 1]) < 0.00000000000001)
    {
        indikG = 1;
        return;
    }
    for (int k = n; k > 0; k--)
    {
        d = 0;
        for (int j = k; j < n; j++)
        {
            s = Matrix[k - 1][j] * x[j];
            d = d + s;
        }
        x[k - 1] = (rightb[k - 1] - d) / Matrix[k - 1][k - 1];
    }
}
template<typename T>
void Invers(T** Matrix, T** Invers, int n)
{
    T* righte;
    T* xe;
    T** ResM;
    ResM = new T * [n];
    xe = new T[n];
    righte = new T[n];
    for (int i = 0; i < n; i++)
    {
        ResM[i] = new T[n];
    }
    for (int i = 0; i < n; i++)
    {
        for (int o = 0; o < n; o++)
        {
            righte[o] = 0;
            for (int p = 0; p < n; p++)
            {
                ResM[o][p] = Matrix[o][p];
            }
        }
        righte[i] = 1;
        Gauss(ResM, righte, n, xe);
        for (int k = 0; k < n; k++)
        {
            Invers[k][i] = xe[k];
        }
    }
    for (int i = 0; i < n; i++) {
        delete[] ResM[i];
    }
    delete[] ResM;
    delete[] righte;
    delete[] xe;

}
template<typename T>
void Transpose(T** Matrix, int n)
{
    T res;
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            res = Matrix[i][j];
            Matrix[i][j] = Matrix[j][i];
            Matrix[j][i] = res;
        }
    }
}
template<typename T>
T NormM1(T** Matrix, int n)
{
    T res;
    T max = 0;
    for (int i = 0; i < n; i++)
    {
        res = 0;
        for (int j = 0; j < n; j++)
        {
            res = res + abs(Matrix[j][i]);
        }
        if (max < res)
        {
            max = res;
        }
    }
    return max;
}
template<typename T>
T NormMInf(T** Matrix, int n)
{
    T res;
    T max = 0;
    for (int i = 0; i < n; i++)
    {
        res = 0;
        for (int j = 0; j < n; j++)
        {
            res = res + abs(Matrix[i][j]);
        }
        if (max < res)
        {
            max = res;
        }
    }
    return max;
}
template<typename T>
void OrtM(T** Matrix, int k, int n, int m, T c, T s)
{
    T* res;

    res = new T[k];
    for (int i = 0; i < k; i++)
    {
        res[i] = c * Matrix[n][i] + s * Matrix[m][i];
        Matrix[m][i] = -s * Matrix[n][i] + c * Matrix[m][i];
    }
    for (int j = 0; j < k; j++)
    {
        Matrix[n][j] = res[j];
    }
    delete[] res;
}
template<typename T>
void QR(T** ET, T** Matrix, T* rightb, int n, T* x)
{
    T c, s, d, rh;
    T r;
    for (int k = 0; k < n - 1; k++)
    {
        for (int j = k + 1; j < n; j++)
        {
            if (abs(Matrix[j][k]) != 0)
            {
                rh = (sqrt(Matrix[k][k] * Matrix[k][k] + Matrix[j][k] * Matrix[j][k]));
                c = Matrix[k][k] / rh;
                s = Matrix[j][k] / rh;
                OrtM(Matrix, n, k, j, c, s);
                OrtM(ET, n, k, j, c, s);
                Matrix[j][k] = 0;
                r = c * rightb[k] + s * rightb[j];
                rightb[j] = -s * rightb[k] + c * rightb[j];
                rightb[k] = r;
            }
        }
    }
    if (abs(Matrix[n - 1][n - 1]) < 0.00000000000001)
    {
        indikQR = 1;
        return;
    }
    for (int k = n; k > 0; k--)
    {
        d = 0;
        for (int j = k; j < n; j++)
        {
            s = Matrix[k - 1][j] * x[j];
            d = d + s;
        }
        x[k - 1] = (rightb[k - 1] - d) / Matrix[k - 1][k - 1];
    }
}
template <typename T>
void test(int set)
{
    T* result;
    T* result1;
    T** ResM;
    T* rightres;
    T** Inv;
    T** ET;
    T** Matrix;
    T* x;
    T* x1;
    T* rightb;
    T** Matrix1;
    T* rightb1;
    T res, max;
    int n, i, j;
    std::ifstream file;
    file.open("test.txt");
    file >> n;
    rightres = new T[n];
    result = new T[n];
    result1 = new T[n];
    Matrix = new T * [n];
    Matrix1 = new T * [n];
    ResM = new T * [n];
    ET = new T * [n];
    Inv = new T * [n];
    x = new T[n];
    x1 = new T[n];
    rightb = new T[n];
    rightb1 = new T[n];
    for (i = 0; i < n; i++) {
        Matrix[i] = new T[n];
        Matrix1[i] = new T[n];
        ET[i] = new T[n];
        Inv[i] = new T[n];
        ResM[i] = new T[n];
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
        {
            file >> Matrix[i][j];
            Matrix1[i][j] = Matrix[i][j];
            ResM[i][j] = Matrix[i][j];
            if (i == j)
            {
                ET[i][j] = 1.0;
            }
            else
            {
                ET[i][j] = 0.0;
            }
        }
        file >> rightb[i];
        rightb1[i] = rightb[i];
        rightres[i] = rightb[i];
    }
    file.close();
    std::ofstream ans;
    ans.open("answerG.txt");
    ans << "Matrix" << std::endl;
    ans << std::setprecision(set);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            ans << Matrix[i][j] << ' ';
        }
        ans << rightb[i] << std::endl;
    }
    Gauss(Matrix, rightb, n, x);
    if (indikG == 1)
    {
        ans << "generate Matrix Gauss" << std::endl;
    }
    else
    {
        ans << std::endl;
        ans << "XGaus" << std::endl;
        for (j = 0; j < n; j++)
        {
            ans << x[j] << ' ';
        }
        ans << std::endl;
        ans << "NormV XGaus" << std::endl;
        ans << NormVect(x, n);
        ans << std::endl;
        ans << "Treg Matrix after Gaus" << std::endl;
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                ans << Matrix[i][j] << ' ';
            }
            ans << rightb[i] << std::endl;
        }
        MultWV(ResM, x, result, n);
        res = 0;
        for (int p = 0; p < n; p++)
        {
            res = res + (rightres[p] - result[p]) * (rightres[p] - result[p]);
        }
        ans << std::endl;
        ans << "discrepancyG=" << sqrt(res) << std::endl;
    }
    QR(ET, Matrix1, rightb1, n, x1);
    if (indikQR == 1)
    {
        ans << "generate Matrix QR" << std::endl;
    }
    else
    {
        ans << std::endl;
        ans << "XQR" << std::endl;
        for (j = 0; j < n; j++)
        {
            ans << x1[j] << ' ';
        }
        ans << std::endl;
        ans << "NormV XQR" << std::endl;
        ans << NormVect(x1, n);
        ans << std::endl;
        ans << "Treg Matrix after QR" << std::endl;
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                ans << Matrix1[i][j] << ' ';
            }
            ans << rightb1[i] << std::endl;
        }
        MultWV(ResM, x1, result1, n);
        res = 0;
        for (int p = 0; p < n; p++)
        {
            res = res + (rightres[p] - result1[p]) * (rightres[p] - result1[p]);
        }
        ans << std::endl;
        ans << "discrepancyQR=" << sqrt(res) << std::endl;
    }
    if (indikQR == 0)
    {
        Transpose(ET, n);
        ans << std::endl;
        ans << "Matrix Q" << std::endl;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                ans << ET[i][j] << ' ';
            }
            ans << std::endl;
        }
        ans << std::endl;
        ans << "Matrix R" << std::endl;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                ans << Matrix1[i][j] << ' ';
            }
            ans << std::endl;
        }
    }
    ans << std::endl;
    if ((indikG == 0) && (indikQR == 0))
    {
        Invers(ResM, Inv, n);
        ans << "Cond1=" << std::fixed << NormM1(ResM, n) * NormM1(Inv, n) << "   Condinf=" << std::fixed << NormMInf(ResM, n) * NormMInf(Inv, n) << std::endl;
        MultM(ResM, Inv, ET, n);
        ans << std::endl;
        ans << "Matrix A^(-1)" << std::endl;
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++) {
                ans << Inv[i][j] << ' ';
            }
            ans << std::endl;
        }
        ans << std::endl;
        ans << "Matrix A*A^(-1)" << std::endl;
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++) {
                ans << ET[i][j] << ' ';
            }
            ans << std::endl;
        }
        res = 0;
        max = 0;
        for (int t = 0; t < n; t++)
        {
            rightres[t] = rightres[t] + 0.01;
            MultWV(Inv, rightres, result, n);
            for (int u = 0; u < n; u++)
            {
                result[u] = result[u] - x[u];
            }
            if (0.01 / NormVect(rightres, n) != 0)
            {
                res = (NormVect(result, n) / NormVect(x, n)) / (0.01 / NormVect(rightres, n));
            }
            if (max < res)
            {
                i = t + 1;
                max = res;
            }
            rightres[t] = rightres[t] - 0.01;
        }
        ans << std::endl;
        ans << "Max deltax/deltaright=" << ' ' << max << ' ' << "Number" << ' ' << i << std::endl;
    }
    else
    {
        ans << "Cond=inf" << std::endl;
    }

    ans.close();
    for (i = 0; i < n; i++) {
        delete[] Inv[i];
        delete[] ET[i];
        delete[] Matrix[i];
        delete[] Matrix1[i];
        delete[] ResM[i];
    }
    delete[] rightres;
    delete[] result;
    delete[] result1;
    delete[] ResM;
    delete[] Inv;
    delete[] ET;
    delete[] x;
    delete[] Matrix;
    delete[] rightb;
    delete[] Matrix1;
    delete[] rightb1;
}
template <typename T>
void test2(int set)
{
    T** result;
    T** ET;
    T** Matrix;
    T** ResM;
    T* x;
    T* rightb;
    T res, max;
    int n, i, j;
    std::ifstream file;
    file.open("test.txt");
    file >> n;
    result = new T * [n];
    Matrix = new T * [n];
    ET = new T * [n];
    ResM = new T * [n];
    x = new T[n];
    rightb = new T[n];
    for (i = 0; i < n; i++)
    {
        Matrix[i] = new T[n];
        ET[i] = new T[n];
        result[i] = new T[n];
        ResM[i] = new T[n];
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
        {
            file >> Matrix[i][j];
            ResM[i][j] = Matrix[i][j];
            if (i == j)
            {
                ET[i][j] = 1.0;
            }
            else
            {
                ET[i][j] = 0.0;
            }
        }
        file >> rightb[i];
    }
    file.close();
    std::ofstream ans;
    ans.open("answer.txt");
    ans << "Matrix" << std::endl;
    ans << std::setprecision(set);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            ans << Matrix[i][j] << ' ';
        }
        ans << std::endl;
    }
    QR(ET, Matrix, rightb, n, x);
    if (indikQR == 0)
    {
        Transpose(ET, n);
        MultM(ET, Matrix, result, n);
        ans << std::endl;
        ans << "Matrix Q" << std::endl;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                ans << ET[i][j] << ' ';
            }
            ans << std::endl;
        }
        ans << std::endl;
        ans << "Matrix R" << std::endl;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                ans << Matrix[i][j] << ' ';
            }
            ans << std::endl;
        }
        ans << std::endl;
        ans << "Matrix Q*R" << std::endl;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                ans << result[i][j] << ' ';
                result[i][j] = result[i][j] - ResM[i][j];
            }
            ans << std::endl;
        }
        ans << std::endl;
        ans << "Matrix A-Q*R" << std::endl;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                ans << result[i][j] << ' ';
            }
            ans << std::endl;
        }
    }
    ans << std::endl;
    ans << "Norm1 A-Q*R" << std::endl;
    ans << NormM1(result, n) << std::endl;
    ans << std::endl;
    ans << "Norminf A-Q*R" << std::endl;
    ans << NormMInf(result, n) << std::endl;
    ans.close();
    for (i = 0; i < n; i++) {
        delete[] ET[i];
        delete[] Matrix[i];
        delete[] result[i];
        delete[] ResM[i];
    }
    delete[] ResM;
    delete[] result;
    delete[] ET;
    delete[] x;
    delete[] Matrix;
    delete[] rightb;
}
int main()
{
    test<double>(15);
    test2<double>(15);
}
