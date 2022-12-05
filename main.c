/**
 * Radix sort LSD parallelized over MPI implementation.
 * @author Dalton Caron
 * @file main.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <mpi.h>
#include <argp.h>

#define MAX(a, b) (a > b ? a : b)

const char *argp_program_version = "parallel radix sort LSD 1.0.0";
const char *argp_program_bug_address = "dpcaron@csu.fullerton.edu";
static char doc[] = "An MPI program for performing a parallelized radix sort least significant digit.";
static char args_doc[] = "[-b base] [-i] [-s serial_threshold]";
static struct argp_option options[] = {
    { "base", 'b', "BASE", 0, "Select the base for the radix sort." },
    { "print_output", 'i', 0, 0, "Prints the sorted list to stdout." },
    { "serial_threshold", 's', "SERIAL_THRES", 0, "Execute the radix sort serially if the input size is less than or equal to this value." },
    { 0 }
};

struct arguments {
    int base;
    char print_output;
    int serial_threshold;
};

struct arguments arguments;

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    struct arguments *arguments = state->input;
    switch(key) {
        case 'b': arguments->base = atoi(arg); break;
        case 'i': arguments->print_output = 1; break;
        case 's': arguments->serial_threshold = atoi(arg); break;
        case ARGP_KEY_ARG: return 0;
        default: return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0 };

/**
 * Swaps two elements.
 * @param i Index to swap with j.
 * @param j Index to swap with i.
 * @param a Data array.
 */
inline void swap(size_t i, size_t j, uint64_t *a) {
    uint64_t temp = a[i];
    a[i] = a[j];
    a[j] = temp;
}

/**
 * Computes the maximum of an array.
 * @param a Data array.
 * @param n Array size.
 */
uint64_t compute_max(const uint64_t *a, const int n) {
    uint64_t max = a[0];
    for (size_t i = 1; i < n; i++) {
        max = MAX(max, a[i]);
    }
    return max;
}

/**
 * Serial radix sort LSD implementation.
 * @param a Data array.
 * @param n Array size.
 */
void radix_sort(uint64_t *a, const int n) {
    const int max = compute_max(a, n);
    signed int digit_place = 1, i;
    int *r = (int *) malloc(sizeof(uint64_t) * n);
    while (max / digit_place > 0) {
        int count[arguments.base];
        memset(count, 0, sizeof(int) * arguments.base);

        for (i = 0; i < n; i++)
            count[(a[i] / digit_place) % arguments.base]++;

        for (i = 1; i < arguments.base; i++)
            count[i] += count[i - 1];

        for (i = n - 1; i >= 0; i--) {
            r[count[(a[i] / digit_place) % arguments.base] - 1] = a[i];
            count[(a[i] / digit_place) % arguments.base]--;
        }

        for (i = 0; i < n; i++)
            a[i] = r[i];

        digit_place *= arguments.base;
    }
    free(r);
}

/**
 * Prints a simple one-dimensional array of integers.
 * @param a Array to print.
 * @param n Size of the array.
 */ 
void print_array(uint64_t *a, const int n) {
    if (arguments.print_output) {
        for (size_t i = 0; i < n; i++) {
            printf("%d ", a[i]);
        }
        printf("\n");
    }
}

/**
 * Reads an array of numbers from the input stream and tracks the maximum.
 * @param arr_out The array to store the output in.
 * @param count_out The size of the resulting array.
 * @param max The maximum element in the resulting array.
 */
void read_input_array(uint64_t **arr_out, int *count_out, uint64_t *max) {
    scanf("%d", count_out);
    *arr_out = (uint64_t *) malloc(*count_out * sizeof(uint64_t));
    int c = *count_out;

    *max = 0;
    while (--c >= 0) {
        scanf("%d", (*arr_out) + c);
        *max = MAX(*max, (*arr_out)[c]);
    }
}

/**
 * Reads the input stream of numbers into an array in the leader (0) process.
 * @param rank The rank of this process.
 * @param a The array to store the data.
 * @param count The size of the resulting array.
 * @param max The largest element in the array. Needed for radix sort.
 */
void leader_read_input(const int rank, uint64_t **a, int *count, uint64_t *max) {
    if (rank == 0) {
        read_input_array(a, count, max);

        // Perform the computation serially if there are not much input.
        if (*count <= arguments.serial_threshold) {
            printf("Executing serial version\n");
            radix_sort(*a, *count);
            print_array(*a, *count);
        } else {
            printf("Executing in parallel\n");
        }
    }
}

/**
 * Not all arrays can be evenly divided among processes. This function calculates how much data each process should 
 * hold and the starting index of the portion held by each process in the displacements array.
 * @param size The size of the communicator.
 * @param send_counts Will contain how much data each process can hold.
 * @param displacements Contains the offset of the processes data relative to the main array.
 * @param count Size of the main array.
 */
void calculate_send_count_and_displacements(const int size, int *send_counts, int *displacements, const int count) {
    const int remainder = count % size;
    const int count_per_process = count / size;

    int increment = 0;
    for (int pid = 0; pid < size; pid++) {
        displacements[pid] = increment;
        // Assign extra iterations from left to right.
        send_counts[pid] = (pid + 1 <= remainder) ? count_per_process + 1 : count_per_process;
        increment += send_counts[pid];
    }
}

/**
 * Evenly divides the data among all processes and sends the data to each process.
 * send_counts and displacements are calculated by the function above this one.
 * @param data The data to scatter.
 * @param send_counts An array the maps each process to how much data it can hold.
 * @param displacements Holds the offset of a given process from the start of the data array.
 * @param local_numbers Where the resulting parition is to be stored for this process.
 * @param local_count The length of the local numbers array.
 */
void scatter_data(const uint64_t *const data, const int *const send_counts, 
    const int *const displacements, uint64_t *local_numbers, const int local_count) {
        MPI_Scatterv(data, send_counts, displacements, MPI_INT64_T, 
            local_numbers, local_count, MPI_INT64_T, 0, MPI_COMM_WORLD);
}

/**
 * Contains the procedure for the parallelized radix sort.
 * Very hard to modularize due to global MPI calls and 
 * I ran out of time to do so.
 */
int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {

        arguments.base = 2;
        arguments.serial_threshold = 300;
        arguments.print_output = 0;

        argp_parse(&argp, argc, argv, 0, 0, &arguments);

    }

    MPI_Bcast(&arguments, sizeof(struct arguments), MPI_BYTE, 0, MPI_COMM_WORLD);

    uint64_t *a;
    uint64_t max;
    int a_count;
    double start, end;

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    leader_read_input(rank, &a, &a_count, &max);

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(&a_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max, 1, MPI_INT64_T, 0, MPI_COMM_WORLD);

    // If input size is too small, algorithm is executed in serial.
    if (a_count <= arguments.serial_threshold) {
        MPI_Barrier(MPI_COMM_WORLD);
        end = MPI_Wtime();
        MPI_Finalize();

        if (rank == 0) {
            printf("Execution time: %lf\n", end - start);
        }

        return 0;
    }

    if (rank != 0) {
        a = (uint64_t *) malloc(sizeof(uint64_t) * a_count);
    }

    // Chunking the input into even sizes is not easy.
    int send_counts[size];
    int displacements[size];
    calculate_send_count_and_displacements(size, send_counts, displacements, a_count);
    const int local_count = send_counts[rank];
    uint64_t local_numbers[local_count];
    uint64_t *r = (uint64_t *) malloc(sizeof(uint64_t) * local_count);
    int64_t digit_place = 1;
    int i;

    scatter_data(a, send_counts, displacements, local_numbers, local_count);

    // This portion of the code starts getting disorganized.
    // I have no excuses and can do better.
    while (max / digit_place > 0) {
        // Sorting step.
        
        int count[arguments.base];
        memset(count, 0, sizeof(int) * arguments.base);

        for (i = 0; i < local_count; i++)
            count[(local_numbers[i] / digit_place) % arguments.base]++;

        for (i = 1; i < arguments.base; i++)
            count[i] += count[i - 1];

        for (i = local_count - 1; i >= 0; i--) {
            r[count[(local_numbers[i] / digit_place) % arguments.base] - 1] = local_numbers[i];
            count[(local_numbers[i] / digit_place) % arguments.base]--;
        }

        for (i = 0; i < local_count; i++) {
            local_numbers[i] = r[i];
            // Count how many digits of the base we have for the redistribution step.
        }

        digit_place *= arguments.base;

        // Data redistribution step.
        // This is the hard part.
        MPI_Allgatherv(local_numbers, local_count, MPI_INT64_T, a, send_counts, displacements, MPI_INT64_T, MPI_COMM_WORLD);

        int current_process = 0;
        int capacities[size];
        int local_count = 0;
        memcpy(capacities, send_counts, sizeof(int) * size);

        for (int base = 0; base < arguments.base; base++) {
            // Scan the list for this base.
            for (int pos = 0; pos < a_count; pos++) {

                // We found a number that needs to be placed. 
                if (((a[pos] / digit_place) % arguments.base) == base) {
                    // Place it in the current process.
                    capacities[current_process]--;

                    // If THIS thread is the current process, copy from the global array into local memory.
                    if (rank == current_process) {
                        local_numbers[local_count] = a[pos];
                        local_count++;
                    }

                    // Start filling the next process if we are out of room.
                    if (capacities[current_process] == 0) {
                        current_process++;
                    }
                }

            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();

    // The array is now sorted. Output.
    if (rank == 0) {
        printf("Execution time: %lf\n", end - start);
        print_array(a, a_count);
    }

    free(r);    
    free(a);

    MPI_Finalize();

    return 0;
}