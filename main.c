#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define DEBUG 1

struct coords
{
    double x;
    double y;
};

struct spline
{
    int dots_count;
    struct coords *base_dots;
    // int pieces_count;
    double *equation_coefficients;
};

// calculate functions
void fill_matrix(double **matrix, double *b_column, struct coords *base_dots, const int n, const int matrix_size);
double *gauss_method(double **matrix, double *b_column, const int matrix_size);
void calcuate_spline_coefs(struct spline *spline1);
double calculate_point(double *coefs, double x, double x0);

// print functions
void print_coords(struct coords *dots, int n);
void print_matrix(double **matrix, double *b_column, int size, char comment[]);
void print_spline_coefs(struct spline a, int number);
void print_spline(struct spline *spline1, int number, double step);

// service functions
void swap(double *a, double *b);

int main()
{
    // Input
    int splines_count;
    printf("Enter the number of splines: ");
    scanf("%d", &splines_count);

    struct spline *splines = (struct spline *)malloc(splines_count * sizeof(struct spline));

    for (int i = 0; i < splines_count; i++)
    {
        printf("\tSpline %d:\nEnter the number of base dots: ", i);
        scanf("%d", &splines[i].dots_count);

        splines[i].base_dots = (struct coords *)malloc(splines[i].dots_count * sizeof(struct coords));

        printf("Enter base dots: x y\n");
        for (int j = 0; j < splines[i].dots_count; j++)
        {
            scanf("%lf %lf", &splines[i].base_dots[j].x, &splines[i].base_dots[j].y);
        }
    }
    // input end

    // solution
    for (int i = 0; i < splines_count; i++)
    {
        calcuate_spline_coefs(&splines[i]);
    }

    // output    
    for (int i = 0; i < splines_count; i++)
    {
        print_spline_coefs(splines[i], i);
        print_spline(&splines[i], i, 0.1);
    }
    return 0;
}

// ------------CALCULATE FUNCTIONS------------
void fill_matrix(double **matrix, double *b_column, struct coords *base_dots, const int n, const int matrix_size)
{
    int curr_row = 0;
    
    // 1. Cплайны проходят через узловые точки
    // S[i] = a[i] + b[i](x-x[i]) + c[i](x-x[i])^2 + d[i](x-x[i])^3
    for (int i = 0; i < n-1; i++)
    {
        // start dot in spline
        b_column[curr_row] = base_dots[i].y;
        matrix[curr_row][i*4] = 1; // a[i]
        curr_row++;
        
        // end dot in spline
        b_column[curr_row] = base_dots[i+1].y;
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
}

double *gauss_method(double **matrix, double *b_column, const int matrix_size)
{
    // https://prog-cpp.ru/gauss/
    double *x_column = (double *)malloc(matrix_size * sizeof(double));
    const double eps = 0.00001; // precision
    double max;
    int index_max;
    int curr = 0;

    while (curr < matrix_size)
    {
        // Поиск строки с максимальным matrix[i][curr]
        max = fabs(matrix[curr][curr]);
        index_max = curr;
        for (int i = curr+1; i < matrix_size; i++)
        {
            if (fabs(matrix[i][curr]) > max)
            {
                max = fabs(matrix[i][curr]);
                index_max = i;
            }
        }

        // Перестановка строк
        for (int j = 0; j < matrix_size; j++)
            swap(&matrix[curr][j], &matrix[index_max][j]);
        swap(&b_column[curr], &b_column[index_max]);

        // Нормализация уравнений
        for (int i = curr; i < matrix_size; i++)
        {
            double divider = matrix[i][curr];
            if (fabs(divider) < eps)
                continue; // На ноль делить нельзя :)
            
            // Деление строк
            for (int j = 0; j < matrix_size; j++)
                matrix[i][j] /= divider;
            b_column[i] /= divider;

            // Вычитание и приведение к ступенчатому виду
            if (i == curr)
                continue; // Не вычитать уравнение само из себя
            for (int j = 0; j < matrix_size; j++)
                matrix[i][j] -= matrix[curr][j];
            b_column[i] -= b_column[curr];
        }

        // DEBUG
        if (DEBUG)
        {
            char comment[50];
            sprintf(comment, "Matrix after %d step", curr);
            print_matrix(matrix, b_column, matrix_size, comment);
            getchar();
        }

        curr++;
    }
    // Обратная подстановка
    for (curr = matrix_size-1; curr >= 0; curr--) 
    {
        x_column[curr] = b_column[curr];
        for (int i = 0; i < curr; i++)
            b_column[i] -= matrix[i][curr] * x_column[curr];
    }

    return x_column;
}

void calcuate_spline_coefs(struct spline *spline1)
{
    const int matrix_size = (spline1->dots_count - 1) * 4;
    double **matrix = (double **)malloc(matrix_size * sizeof(double *));
    for (int i = 0; i < matrix_size; i++)
    {
        matrix[i] = (double *)malloc((matrix_size+1) * sizeof(double));
        memset(matrix[i], 0, (matrix_size+1) * sizeof(double)); // sets array values to 0
    }
    double *b_column = (double *)malloc(matrix_size * sizeof(double));
    memset(b_column, 0, matrix_size * sizeof(double)); // sets array values to 0

    // solution
    fill_matrix(matrix, b_column, spline1->base_dots, spline1->dots_count, matrix_size);
    if (DEBUG)
        print_matrix(matrix, b_column, matrix_size, "Calculated matrix");

    spline1->equation_coefficients = gauss_method(matrix, b_column, matrix_size);
}

double calculate_point(double *coefs, double x, double x0)
{
    return coefs[0] + coefs[1]*(x - x0) + coefs[2]*pow(x - x0, 2) + coefs[3]*pow(x - x0, 3);
}
// ----------CALCULATE FUNCTIONS END----------

// ------------PRINT FUNCTIONS------------
void print_coords(struct coords *dots, int n)
{
    for (int i = 0; i < n; i++)
    {
        printf("[DEBUG]: %lf %lf\n", dots[i].x, dots[i].y);
    }
}

void print_matrix(double **matrix, double *b_column, int size, char comment[])
{
    
    if (strlen(comment) == 0)
        comment = "Matrix";

    printf("[DEBUG]: ------------%s------------\n", comment);
    for (int i = 0; i < size; i++)
    {
        printf("[DEBUG]: ");
        for (int j = 0; j < size; j++)
        {
            printf("%8.5lf ", matrix[i][j]);
        }
        printf("| %8.5lf\n", b_column[i]);
    }
    printf("[DEBUG]: ----------%s end----------\n\n", comment);
}

void print_spline_coefs(struct spline a, int number)
{
    printf("------------Spline %d coefs------------\n", number);
    for (int i = 0; i < a.dots_count - 1; i++)
    {
        printf("Piece %d: ", i);
        for (int j = 0; j < 4; j++)
            printf("%8.5lf ", a.equation_coefficients[i*4+j]);
        printf("\n");
    }
    printf("----------Spline %d coefs end----------\n\n", number);
}

void print_spline(struct spline *spline1, int number, double step)
{
    printf("------------Spline %d graph------------\n", number);
    int current_piece = 0;
    for (double x = spline1->base_dots[0].x; x <= spline1->base_dots[spline1->dots_count-1].x; x += step)
    {
        if (x >= spline1->base_dots[current_piece+1].x)
            current_piece++;
        
        printf("%5.2lf => %8.5lf\n", x, calculate_point(&spline1->equation_coefficients[current_piece*4], x, spline1->base_dots[current_piece].x));
    }
    printf("----------Spline %d graph end----------\n", number);
}
// ----------PRINT FUNCTIONS END------------

// ------------SERVICE FUNCTIONS------------
void swap(double *a, double *b)
{
    double temp = *a;
    *a = *b;
    *b = temp;
}
// ----------SERVICE FUNCTIONS END----------