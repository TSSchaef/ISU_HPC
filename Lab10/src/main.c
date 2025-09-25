#include <stdio.h>
#include "quad.h"


int main(int argc, char *argv[]) {
    quad q;

    printf("Enter the coordinates of four points (x y) denoting a parallelogram:\n");
    for (int i = 0; i < 4; ++i) {
        printf("Point %d: ", i+1);
        if(scanf("%lf %lf", &q.corners[i].x, &q.corners[i].y) != 2){
            fprintf(stderr, "Error on scanf input, exiting\n");
            return 1;
        }
    }

    printf("The area of the quadrilateral is: %lf\n", quad_area(&q));
    printf("The perimeter of the quadrilateral is: %lf\n", quad_perimeter(&q));
    
    return 0;
}
