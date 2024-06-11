#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <chrono>

const int n = 12000;//размерность матрицы
double e = 10e-10;
double A[n][n];
double ipsilon[n];
double b[n];
double r0[n],x0[n],z0[n],alpha,beta,Az[n];
double temp[n];

double scalMul(double u[n], double v[n]) {
    double result = 0;
    for (int i = 0; i < n; i++) {
        result += u[i] * v[i];
    }
    return result;
}

void MatrixMulVector(double result[], double matrix[n][n], double vector[]) {
    double tmp;
   for (int str = 0; str < n; str++) {
        tmp = 0;
        for (int col = 0; col < n; col++)
            tmp += matrix[str][col]*vector[col];
        result[str] = tmp;
    }
}


void printMatrix(double matrix[n][n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%lf\t", matrix[i][j]);
        }
    }
    printf("\n");
}
void printVector(double vector[]) {
    for (int j = 0; j < n; j++) {
        printf("%lf  ", vector[j]);

    }
    printf("\n");
}

void clsVec(double vec[n]) {
    for (int i = 0; i < n; i++)vec[i] = 0;
}

void cpyVec(double source[n], double dest[n]) {
    for (int i = 0; i < n; i++) dest[i] = source[i];
}

int main(int argc, char** argv)
{
    double sum1 , sum2=0, rScal ;
    int iter=0;
    //заполняем матрицу
    auto beginTime = std::chrono::steady_clock::now();
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            if (i == j)
                A[i][j] = 2.;
            else A[i][j] = 1;
        }
    //заполняем вектор u
    for (int i = 0; i < n; i++)
        ipsilon[i] = sin(2 * M_PI * i / n);
    //заполняем свободные коэффициенты
    MatrixMulVector(b,A,ipsilon);
    //r0 = b-Ax0
    clsVec(temp);
    MatrixMulVector(temp, A, x0);
    for (int i = 0; i < n; i++) 
        r0[i] = b[i] - temp[i];
    cpyVec(r0, z0);
    //считаем знаменатель останова до цикла
    for (int i = 0; i < n; i++) 
        sum2 += b[i]*b[i];
    sum2 = sqrt(sum2);
    
    //основной цикл
    do{
        clsVec(Az);
        MatrixMulVector(Az, A, z0);
        rScal = scalMul(r0, r0);
        alpha = rScal / scalMul(Az, z0);
        //xn+1 = xn+a*zn
        for (int i = 0; i < n; i++) {
            x0[i] = x0[i] + alpha * z0[i];
        }         
        //rn+1=rn-alpha*Az
        for (int i = 0; i < n; i++)
            r0[i] = r0[i] - alpha * Az[i];
        //beta=(rn+1,rn+1)/(rn,rn)
        beta = scalMul(r0, r0) /rScal;
        //zn+1=rn+1+beta*z
        for (int i = 0; i < n; i++)
            z0[i] = r0[i] + beta * z0[i];
        sum1 = 0;
        //критерий завершения
        for (int i = 0; i < n; i++) {
            sum1 += r0[i]*r0[i];
        }
        sum1 = sqrt(sum1);
        ++iter;
    }while(sum1 / sum2 > e);
    auto endTime = std::chrono::steady_clock::now();
    auto elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - beginTime);
    printf("Iteration: %d\t from %lf sec\n", iter - 1,elapsed_ns.count() / 1000000000.);
    //printVector(ipsilon);
    //printVector(x0);
    
}


