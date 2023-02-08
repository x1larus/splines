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

    print_matrix(matrix, x_column, matrix_size);
    return 0;
}