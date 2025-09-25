#ifndef QUAD_H
#define QUAD_H

#include <math.h>

typedef struct {
    double x;
    double y;
} point;

typedef struct {
    point corners[4];
    double perimeter;
    double area;
} quad;

double quad_perimeter(quad *q);
double quad_area(quad *q);

inline double distance(point *a, point *b) {
    double dx = b->x - a->x;
    double dy = b->y - a->y;
    return sqrt(dx * dx + dy * dy);
}

#endif
