#include "splines.h"
#include <stdio.h>
#include <stdlib.h>


int main()
{
    // Input
    int splines_count;
    printf("Enter the number of splines: ");
    scanf("%d", &splines_count);

    Spline *splines = (Spline*)malloc(splines_count * sizeof(Spline));

    for (int i = 0; i < splines_count; i++)
    {
        printf("\tSpline %d:\nEnter the number of base dots: ", i);
        scanf("%d", &splines[i].dots_count);

        splines[i].base_dots = (struct coords *)malloc(splines[i].dots_count * sizeof(struct coords));
        
        printf("Please enter dots in ascending order!\n");
        for (int j = 0; j < splines[i].dots_count; j++)
        {
            printf("Enter dot %d (x y): ", j);
            scanf("%lf %lf", &splines[i].base_dots[j].x, &splines[i].base_dots[j].y);
        }
        printf("\n");
    }

    // calculating
    for (int i = 0; i < splines_count; i++)
        calcuate_spline(&splines[i]);

    // printing   
    for (int i = 0; i < splines_count; i++)
    {
        print_spline(&splines[i], i);
        // print_graph(&splines[i], i, 0.1);
    }

    // get crossing points or min distance
    for (int i = 0; i < splines_count-1; i++)
    {
        for (int j = i+1; j < splines_count; j++)
        {
            printf("------------Cross between spline %d and %d------------\n", i, j);

            // максимум корней: кол-во уравнений 1 * кол-во уравнений 2 * 3
            Coords *x = (Coords *)malloc((splines[i].dots_count-1)*(splines[j].dots_count-1)*3*sizeof(Coords)); 
            int cross_count = get_intersection_points(x, &splines[i], &splines[j]);
            if (cross_count == 0)
            {
                printf("Splines do not intersect or equal\n");
                // Get min distance
                double dist = get_min_distance(&splines[i], &splines[j]);
                if (dist != -1)
                    printf("Min distance between splines is %lf\n", dist);
                else
                    printf("Min distance can't be found :(\n");
            }
            else
            {
                for (int k = 0; k < cross_count; k++)
                {
                    printf("x = %8.5lf; y = %8.5lf\n", x[k].x, x[k].y);
                }
            }

            free(x);
            printf("----------Cross between spline %d and %d end----------\n\n", i, j);
        }
    }

    // Print graphs in console
    // for (int i = 0; i < splines_count; i++)
    // {
    //     printf("\tSpline %d graph\n", i);
    //     print_real_graph(splines[i]);
    // }

    // freeing memory
    for (int i = 0; i < splines_count; i++)
    {
        free(splines[i].base_dots);
        free(splines[i].coefs);
    }
    free(splines);
    return 0;
}