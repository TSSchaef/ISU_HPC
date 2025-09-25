#include "quad.h"

double quad_perimeter(quad *q){
    double peri = 0.0;
    for (int i = 0; i < 4; ++i) {
        int j = (i + 1) % 4;
        peri += distance(&(q->corners[i]), &(q->corners[j]));
    }
    return peri;
}
