#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define e 2.718281828459
// Keeping taylor depth low because the recursive factorial is inneficient
#define TAYLOR_DEPTH 20 

double *arr;
int arr_size;

int factorial(int n){
    if (n < 0) {
        fprintf(stderr, "Error: Factorial of a negative number is undefined.\n");
        return -1; 
    }
    if (n == 0) return 1;

    return n * factorial(n - 1);
}

double exponentiate(double x){
    int x0 = (int) x;
    double z = x - x0;

    double accumulator = 1 + z;

    int i;
    for( i = 2; i < TAYLOR_DEPTH; i++){
        accumulator += pow(z, i) / factorial(i);
    }

    return accumulator * pow(e, x0);
}

void get_input(){
    printf("Enter the size of the array: ");
    scanf("%d", &arr_size);
    arr = (double *) malloc(arr_size * sizeof(double));

    printf("Enter %d real numbers:\n", arr_size);
    for(int i = 0; i < arr_size; i++){
        scanf("%lf", &arr[i]);
    }

}

int main(int argc, char *argv[]){
    get_input();

    FILE *fp = fopen("data.txt", "w");

    if( !fp ){
        fprintf(stderr, "Error opening file\n");
        return 1;
    }

    printf("Array: \n[");
    fprintf(fp, "Array: \n[");

    int i;
    for (i = 0; i < arr_size; i++){
        printf("%lf", arr[i]);
        fprintf(fp, "%lf", arr[i]);

        if( i != arr_size - 1){

            printf(", ");
            fprintf(fp, ", ");
        } else {

            printf("]\n\n");
            fprintf(fp, "]\n\n");
        }
    }


    for(i = 0; i < arr_size; i++){
        printf("e^%lf = %.6f\n", arr[i], exponentiate((double) arr[i]));
        fprintf(fp, "e^%lf = %.6f\n", arr[i], exponentiate((double) arr[i]));
    }


    free(arr);
    fclose(fp);

    return 0;
}
