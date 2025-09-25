#include "quad.h"

double quad_area(quad *q){
    double area = 0.0;
    for (int i = 0; i < 4; ++i) {
        int j = (i + 1) % 4;
        area += (q->corners[i].x * q->corners[j].y) - (q->corners[j].x * q->corners[i].y);
    }
    return fabs(area) / 2.0;
}
