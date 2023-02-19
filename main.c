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

    // solution
    for (int i = 0; i < splines_count; i++)
        calcuate_spline(&splines[i]);

    // output    
    for (int i = 0; i < splines_count; i++)
    {
        print_spline(splines[i], i);
        print_graph(&splines[i], i, 0.1);
    }

    for (int i = 0; i < splines_count; i++)
    {
        free(splines[i].base_dots);
        free(splines[i].coefs);
    }
    free(splines);
    return 0;
}