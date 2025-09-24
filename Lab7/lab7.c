#include <stdio.h>
#include <math.h>

int factorial(int n) {
    if( n < 0 ) {
        fprintf(stderr, "Error: Factorial of a negative number is undefined.\n");
        return -1;

    } else if(n == 0) {
        return 1;

    } else {
        int i;
        for(i = n - 1; i >= 1; i--) {
            n *= i;
        }
        return n;

    }
}


int main(int argc, char *argv[]) {
    printf("Factorial of 5 is %d\n", factorial(5));
    printf("Factorial of 12 is %d\n", factorial(12));

    printf("Exponential of 3.14 is %f\n", exp(3.14));
    printf("Logarithm of 3.14 is %f\n", log(3.14));

    return 0;
}
