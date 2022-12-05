#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char *argv[]) {

    if (argc < 3) {
        printf("Expected two integer arguments: one for the amount of random numbers and one for the maximum value that can be generated.\n");
        return 0;
    }

    uint64_t count = atoi(argv[1]);
    uint64_t max = atoi(argv[2]);

    time_t t;
    srand(time(&t));

    FILE *file = NULL;
    file = fopen("data.in", "w");
    if (file == NULL) {
        printf("Failed to open output file.\n");
        return 0;
    }

    fprintf(file, "%d ", count);
    uint64_t r;
    while (count > 0) {
        r = rand() % max;
        fprintf(file, "%d ", r);
        count--;
    }

    fclose(file);
    file = NULL;
}