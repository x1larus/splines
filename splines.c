#include "splines.h"

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define M_PI (3.141592653589793)
#define M_2PI (2.*M_PI)

#define DEBUG 0


//===========================================FUNCTIONS DECLARATIONS===========================================

// ------------CALCULATE FUNCTIONS------------

// Filling the matrix of the system of linear algebraic equations
void fill_matrix(double **matrix, double *b_column, Coords *base_dots, const int n, const int matrix_size);

// Solve the system of linear algebraic equations (Gauss)
double *solve_matrix(double **matrix, double *b_column, const int matrix_size);

// Calculate f(x)
double calculate_point(double *coefs, double x, double x0);

// Solve x^3 + ax^2 + bx + c = 0 (Cardano)
// Returns the number of roots
int solve_cubic(double *x, double a, double b, double c);

// Solve ax^2 + bx + c = 0
// Returns the number of roots
int solve_quadratic(double *x, double a, double b, double c);

// ----------CALCULATE FUNCTIONS END----------


// ------------SERVICE FUNCTIONS------------

// Print the matrix of the system of linear algebraic equations
void print_matrix(double **matrix, double *b_column, int size, char comment[]);

void swap(double *a, double *b);

double max(double a, double b);

double min(double a, double b);

// ----------SERVICE FUNCTIONS END----------

//=========================================FUNCTIONS DECLARATIONS END=========================================


//===========================================FUNCTIONS DEFENITIONS===========================================

// ------------GLOBAL FUNCTIONS------------

void calcuate_spline(Spline *spline1)
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

    spline1->coefs = solve_matrix(matrix, b_column, matrix_size);

    // Delete
    for (int i = 0; i < matrix_size; i++)
        free(matrix[i]);
    free(matrix);
    free(b_column);
}

void print_spline(Spline a, int number)
{
    printf("------------Spline %d coefs------------\n", number);
    for (int i = 0; i < a.dots_count - 1; i++)
    {
        printf("Piece %d: ", i);
        for (int j = 0; j < 4; j++)
            printf("%8.5lf ", a.coefs[i*4+j]);
        printf("\n");
    }
    printf("----------Spline %d coefs end----------\n\n", number);
}

void print_graph(Spline *spline1, int number, double step)
{
    printf("------------Spline %d graph------------\n", number);
    int current_piece = 0;
    for (double x = spline1->base_dots[0].x; x <= spline1->base_dots[spline1->dots_count-1].x; x += step)
    {
        if (x >= spline1->base_dots[current_piece+1].x)
            current_piece++;
        
        printf("x = %5.2lf; y = %8.5lf\n", x, calculate_point(&spline1->coefs[current_piece*4], x, spline1->base_dots[current_piece].x));
    }
    printf("----------Spline %d graph end----------\n\n", number);
}

int get_intersection_points(Coords *x, Spline s1, Spline s2)
{
    int curr = 0;
    double *roots = (double *)malloc(3*sizeof(double));
    int roots_count;

    for (int i = 0; i < s1.dots_count-1; i++)
    {
        for (int j = 0; j < s2.dots_count-1; j++)
        {
            double left_border = max(s1.base_dots[i].x, s2.base_dots[j].x);
            double right_border = min(s1.base_dots[i+1].x, s2.base_dots[j+1].x);
            if (left_border <= right_border)
            {
                // ax^3 + bx^2 + cx + d
                double a = s1.coefs[i*4+3];
                a -= s2.coefs[j*4+3];

                double b = (s1.coefs[i*4+2] - 3*s1.coefs[i*4+3]*s1.base_dots[i].x); 
                b -= (s2.coefs[j*4+2] - 3*s2.coefs[j*4+3]*s2.base_dots[j].x);
                
                double c = s1.coefs[i*4+1] - 2*s1.coefs[i*4+2]*s1.base_dots[i].x + 3*s1.coefs[i*4+3]*pow(s1.base_dots[i].x, 2);
                c -= (s2.coefs[j*4+1] - 2*s2.coefs[j*4+2]*s2.base_dots[j].x + 3*s2.coefs[j*4+3]*pow(s2.base_dots[j].x, 2));

                double d = s1.coefs[i*4] - s1.coefs[i*4+1]*s1.base_dots[i].x + s1.coefs[i*4+2]*pow(s1.base_dots[i].x, 2) - s1.coefs[i*4+3]*pow(s1.base_dots[i].x, 3);
                d -= (s2.coefs[j*4] - s2.coefs[j*4+1]*s2.base_dots[j].x + s2.coefs[j*4+2]*pow(s2.base_dots[j].x, 2) - s2.coefs[j*4+3]*pow(s2.base_dots[j].x, 3));

                if (a != 0)
                {
                    // кубический случай
                    // приводим уравнение для Кардано
                    b /= a;
                    c /= a;
                    d /= a;
                    roots_count = solve_cubic(roots, b, c, d);
                }
                else if (b != 0)
                {
                    // квадратичный случай
                    roots_count = solve_quadratic(roots, b, c, d);
                }
                else if (c != 0)
                {
                    // линейный случай
                    roots_count = 1;
                    roots[0] = -d / c;
                }
                else
                {
                    roots_count = 0;
                }

                for (int k = 0; k < roots_count; k++)
                {
                    if (roots[k] >= left_border && roots[k] <= right_border)
                    {
                        // проверка совпадающих решений на концах областей определения функций
                        if (curr > 0 && x[curr-1].x == roots[k])
                            continue;
                        
                        x[curr].x = roots[k];
                        x[curr++].y = calculate_point(&s1.coefs[i*4], roots[k], s1.base_dots[i].x);
                    }
                }
            }
        }
    }

    free(roots);

    return curr > 0 ? curr : 0;
}

void print_real_graph(Spline s1)
{
    // Грязь))
    int cols;
    system("tput cols >> temp.txt"); // require ncurses library
    FILE *f = fopen("temp.txt", "r");
    fscanf(f, "%d", &cols);
    fclose(f);
    system("rm temp.txt");

    for (int i = 0; i < cols; i++)
        printf("-");
    printf("\n");

    double x_step = (s1.base_dots[s1.dots_count-1].x - s1.base_dots[0].x) / cols;
    Coords *graph = (Coords *)malloc(cols*sizeof(Coords));
    int current_piece = 0;
    int curr = 0;
    double max_y, min_y;
    short flag_first = 1;
    for (double x = s1.base_dots[0].x; x <= s1.base_dots[s1.dots_count-1].x; x += x_step)
    {
        if (x >= s1.base_dots[current_piece+1].x)
            current_piece++;
        graph[curr].x = x;
        graph[curr].y = calculate_point(&s1.coefs[current_piece*4], x, s1.base_dots[current_piece].x);
        
        if (flag_first)
        {
            flag_first = 0;
            max_y = graph[curr].y;
            min_y = graph[curr].y;
        }
        else
        {
            max_y = max(max_y, graph[curr].y);
            min_y = min(min_y, graph[curr].y);
        }

        curr++;
    }

    int rows = cols / 4;
    double y_step = (max_y - min_y) / rows;
    char screen_buffer[rows+1][cols+1];
    for (int i = 0; i < rows+1; i++)
    {
        memset(screen_buffer[i], ' ', cols);
        screen_buffer[i][cols] = '\0';
    }

    for (int i = 0; i < cols; i++)
    {
        screen_buffer[rows - (int)round((graph[i].y - min_y) / y_step)][i] = 'o';
    }

    for (int i = 0; i < rows+1; i++)
    {
        printf("%s\n", screen_buffer[i]);
    }

    for (int i = 0; i < cols; i++)
        printf("-");
    printf("\n");

    free(graph);
}

// ----------GLOBAL FUNCTIONS END----------


// ------------CALCULATE FUNCTIONS------------

void fill_matrix(double **matrix, double *b_column, Coords *base_dots, const int n, const int matrix_size)
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

double *solve_matrix(double **matrix, double *b_column, const int matrix_size)
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
            sprintf(comment, "Matrix after %d gauss step", curr);
            print_matrix(matrix, b_column, matrix_size, comment);
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

double calculate_point(double *coefs, double x, double x0)
{
    return coefs[0] + coefs[1]*(x - x0) + coefs[2]*pow(x - x0, 2) + coefs[3]*pow(x - x0, 3);
}

int solve_cubic(double *x, double a, double b, double c)
{
    // http://algolist.ru/maths/findroot/cubic.php
    double q, r, r2, q3;
    q = (a * a - 3. * b) / 9.;
    r = (a * (2. * a * a - 9. * b) + 27. * c) / 54.;
    r2 = r * r;
    q3 = q * q * q;
    if (r2 < q3)
    {
        double t = acos(r / sqrt(q3));
        a /= 3.;
        q = -2. * sqrt(q);
        x[0] = q * cos(t / 3.) - a;
        x[1] = q * cos((t + M_2PI) / 3.) - a;
        x[2] = q * cos((t - M_2PI) / 3.) - a;
        return (3);
    }
    else
    {
        double aa, bb;
        if (r <= 0.)
            r = -r;
        aa = -pow(r + sqrt(r2 - q3), 1. / 3.);
        if (aa != 0.)
            bb = q / aa;
        else
            bb = 0.;
        a /= 3.;
        q = aa + bb;
        r = aa - bb;
        x[0] = q - a;
        x[1] = (-0.5) * q - a;
        x[2] = (sqrt(3.) * 0.5) * fabs(r);
        if (x[2] == 0.)
            return (2);
        return (1);
    }
}

int solve_quadratic(double *x, double a, double b, double c)
{
    double d;
    /* the main case */
    d = b * b - 4. * a * c; /* the discriminant */
    /* one distinct root */
    if (d == 0.)
    {
        x[0] = x[1] = -b / (2. * a);
        return (1);
    }
    if (d < 0.)
    {
        return (0);
    }
    /* 2 real roots: avoid subtraction of 2 close numbers */
    if (b >= 0.)
        d = (-0.5) * (b + sqrt(d));
    else
        d = (-0.5) * (b - sqrt(d));
    x[0] = d / a;
    x[1] = c / d;
    return (2);
}

// ----------CALCULATE FUNCTIONS END----------


// ------------SERVICE FUNCTIONS------------

void swap(double *a, double *b)
{
    double temp = *a;
    *a = *b;
    *b = temp;
}

double max(double a, double b)
{
    return a > b ? a : b;
}

double min(double a, double b)
{
    return a > b ? b : a;
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
    printf("[DEBUG]: Press enter to continue...\n");
    getchar();
}

// ----------SERVICE FUNCTIONS END----------

//=========================================FUNCTIONS DEFENITIONS END=========================================