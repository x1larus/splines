#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct coords
{
    double x;
    double y;
};

void print_coords(struct coords *dots, int n)
{
    for (int i = 0; i < n; i++)
    {
        printf("%lf %lf\n", dots[i].x, dots[i].y);
    }
}

void print_matrix(double **matrix, double *x_col, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf("%5.2lf ", matrix[i][j]);
        }
        printf("| %5.2lf\n", x_col[i]);
    }
}

int main()
{
    int n;
    scanf("%d", &n);
    struct coords *base_dots = (struct coords*)malloc(n*sizeof(struct coords));
    for (int i = 0; i < n; i++)
    {
        scanf("%lf %lf", &base_dots[i].x, &base_dots[i].y);
    }

    // solution
    const int matrix_size = (n - 1) * 4;
    double **matrix = (double **)malloc(matrix_size * sizeof(double *));
    for (int i = 0; i < ((n - 1) * 4); i++)
    {
        matrix[i] = (double *)malloc(matrix_size * sizeof(double));
        memset(matrix[i], 0, matrix_size * sizeof(double));
    }
    double *x_column = (double *)malloc(matrix_size * sizeof(double));
    memset(x_column, 0, matrix_size * sizeof(double));
    
    int curr_row = 0;
    
    // 1. Cплайны проходят через узловые точки
    // S[i] = a[i] + b[i](x-x[i]) + c[i](x-x[i])^2 + d[i](x-x[i])^3
    for (int i = 0; i < n-1; i++)
    {
        // start dot in spline
        x_column[curr_row] = base_dots[i].y;
        matrix[curr_row][i*4] = 1; // a[i]
        curr_row++;
        
        // end dot in spline
        x_column[curr_row] = base_dots[i+1].y;
        matrix[curr_row][i*4] = 1; // a[i]
        matrix[curr_row][i*4+1] = base_dots[i+1].x - base_dots[i].x; // b[i]
        matrix[curr_row][i*4+2] = pow(base_dots[i+1].x - base_dots[i].x, 2); // c[i]
        matrix[curr_row][i*4+3] = pow(base_dots[i+1].x - base_dots[i].x, 3); // d[i]
        curr_row++;
    }

    // 2. Плавность в стыке сплайнов
    // S'[i] = b[i] + 2c[i](x-x[i]) + 3d[i](x-x[i])^2
    // S''[i] = 2c[i] + 6d[i](x-x[i])
    for (int i = 0; i < n-2; i++)
    {
        // S'[i] = S'[i+1] <=> b[i] + 2c[i](x[i+1]-x[i]) + 3d[i](x[i+1]-x[i])^2 -
        // - b[i+1] - 2c[i+1](x[i+1]-x[i+1]) - 3d[i+1](x[i+1]-x[i+1])^2
        matrix[curr_row][i*4+1] = 1; // b[i]
        matrix[curr_row][i*4+2] = 2*(base_dots[i+1].x - base_dots[i].x); // c[i]
        matrix[curr_row][i*4+3] = 3*pow(base_dots[i+1].x - base_dots[i].x, 2); // d[i]
        matrix[curr_row][(i+1)*4+1] = -1; // b[i+1]
        curr_row++;

        // S''[i] = S''[i+1] <=> 2c[i] + 6d[i](x[i+1]-x[i]) - 2c[i+1] - 6d[i+1](x[i+1]-x[i+1])
        matrix[curr_row][i*4+2] = 2; // c[i]
        matrix[curr_row][i*4+3] = 6*(base_dots[i+1].x - base_dots[i].x); // d[i]
        matrix[curr_row][(i+1)*4+2] = -2; // c[i+1]
        curr_row++;
    }

    // 3. Нулевая кривизна в крайних точках
    // S''[0] = 0 в точке x[0] <=> 2c[0] + 6d[i](x[0]-x[0])
    matrix[curr_row][2] = 2; // c[0]
    curr_row++;

    // S''[n-2] = 0 в точке x[n-1] <=> 2c[n-2] + 6d[n-2](x[n-1]-x[n-2])
    matrix[curr_row][matrix_size - 4 + 2] = 2; // c[n-2]
    matrix[curr_row][matrix_size - 4 + 3] = 6*(base_dots[n-1].x - base_dots[n-2].x); // d[n-2]
    curr_row++;

    print_matrix(matrix, x_column, matrix_size);
    return 0;
}