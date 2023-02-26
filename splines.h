#pragma once

// An object that stores point coordinates
typedef struct coords
{
    double x;
    double y;
} Coords;

// Spline object
typedef struct spline
{
    int dots_count;
    Coords *base_dots;
    double *coefs;
} Spline;

// Calculate equations
void calcuate_spline(Spline *spline1);

// Prints equations coefficients
// Param NUMBER uses only in comment
void print_spline(Spline *a, int number);

// Prints f(x) with a STEP difference
// Param NUMBER uses only in comment
void print_graph(Spline *spline1, int number, double step);

// Fill array of two splines intersection points 
// Returns the number of intersections
int get_intersection_points(Coords *x, Spline s1, Spline s2);

// Print spline graph in console
void print_real_graph(Spline s1);