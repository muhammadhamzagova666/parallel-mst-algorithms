// Language: C
// ----------------------------------------------------------------------------
// File: kruskal_omp.c
// Description: Implements a parallel version of Kruskal's algorithm using OpenMP
//              for finding the Minimum Spanning Tree (MST) of a weighted graph.
//              The algorithm sorts edges with a parallel merge sort and then uses a
//              union-find structure to select MST edges while avoiding cycles.
// Target Users: Developers and researchers interested in parallel graph algorithms.
// ----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <assert.h>

#define MAX_THREADS 8  ///< Maximum number of OpenMP threads to use.

/**
 * @brief Structure representing a weighted edge.
 *
 * Each edge connects two vertices (u and v) with a certain weight (w).
 */
typedef struct edge {
    int u, v, w;
} edge;

int n;                ///< Number of nodes in the graph.
edge* edges;          ///< Array to hold all generated edges.
edge* chosen_edges;   ///< Array to store the edges chosen for the MST.
int* par;             ///< Parent array for union-find operations.
int num_edge = 0;     ///< Total number of generated edges.

/**
 * @brief Finds the set representative for a vertex using path compression.
 *
 * Recursively finds the root of the set and compresses the path for
 * fast future lookups.
 *
 * @param u The vertex for which to find the set representative.
 * @return int The representative of the set containing vertex u.
 */
int find_set(int u) {
    if (par[u] == u)
        return u;
    return par[u] = find_set(par[u]);  // Apply path compression.
}

/**
 * @brief Merges the sets containing two vertices.
 *
 * Uses union-find to connect the sets if they are disjoint and prevents
 * forming cycles. Returns a non-zero value on successful merge.
 *
 * @param u First vertex.
 * @param v Second vertex.
 * @return int 1 if the merge was successful; otherwise, 0.
 */
int merge_set(int u, int v) {
    int pu = find_set(u), pv = find_set(v);
    if (pu == pv)
        return 0;
    par[pv] = pu;  // Merge by setting representative.
    return 1;
}

/**
 * @brief Comparator function based on edge weight with tie-breakers.
 *
 * If weights are equal, it compares source vertices; if they are also equal,
 * it compares destination vertices.
 *
 * @param x Pointer to the first edge.
 * @param y Pointer to the second edge.
 * @return int Non-zero if x should come before y, 0 otherwise.
 */
int comparison_weight(edge* x, edge* y) {
    if (x->w == y->w) {
        if (x->u == y->u)
            return x->v < y->v;
        return x->u < y->u;
    }
    return x->w < y->w;
}

/**
 * @brief Comparator function based solely on node order.
 *
 * Orders edges first by source vertex and then by destination vertex.
 *
 * @param x Pointer to the first edge.
 * @param y Pointer to the second edge.
 * @return int Non-zero if x should come before y, 0 otherwise.
 */
int comparison_node(edge* x, edge* y) {
    if (x->u == y->u)
        return x->v < y->v;
    return x->u < y->u;
}

/**
 * @brief Merges two sorted subarrays into one sorted array.
 *
 * This helper function merges two already sorted arrays (larr and rarr) into a
 * single sorted array using the provided comparison function.
 *
 * @param edges Destination array where merged result will be stored.
 * @param larr Left sorted subarray.
 * @param nl Number of elements in the left subarray.
 * @param rarr Right sorted subarray.
 * @param nr Number of elements in the right subarray.
 * @param comparison Function pointer to determine order.
 */
void merge(edge edges[], edge larr[], int nl, edge rarr[], int nr, int (*comparison)(edge*, edge*)) {
    int il = 0, ir = 0, j = 0;
    // Merge two subarrays until one is exhausted.
    while (il < nl && ir < nr) {
        if ((*comparison)(&larr[il], &rarr[ir])) {
            edges[j] = larr[il];
            il++;
        } else {
            edges[j] = rarr[ir];
            ir++;
        }
        j++;
    }
    // Copy any remaining elements from left subarray.
    while (il < nl) {
        edges[j++] = larr[il++];
    }
    // Copy any remaining elements from right subarray.
    while (ir < nr) {
        edges[j++] = rarr[ir++];
    }
}

/**
 * @brief Parallel merge sort for an array of edges.
 *
 * Recursively sorts the array of edges using merge sort in parallel with OpenMP.
 *
 * @param edges Array of edges to be sorted.
 * @param n Number of elements in the array.
 * @param comparison Function pointer to compare two edges.
 */
void merge_sort(edge edges[], int n, int (*comparison)(edge*, edge*)) {
    if (n > 1) {
        int m = n / 2;
        edge* larr = malloc(m * sizeof(edge));
        edge* rarr = malloc((n - m) * sizeof(edge));
        memcpy(larr, edges, m * sizeof(edge));
        memcpy(rarr, edges + m, (n - m) * sizeof(edge));

        // Parallelize recursive calls with OpenMP tasks.
        #pragma omp parallel
        {
            #pragma omp single
            {
                #pragma omp task shared(larr, m, comparison)
                merge_sort(larr, m, comparison);

                #pragma omp task shared(rarr, n, m, comparison)
                merge_sort(rarr, n - m, comparison);
                
                #pragma omp taskwait  // Ensure both tasks complete.
            }
        }

        // Merge the two sorted halves.
        merge(edges, larr, m, rarr, n - m, comparison);
        free(larr);
        free(rarr);
    }
}

/**
 * @brief Main function to compute the MST using parallel Kruskal's algorithm.
 *
 * Generates a random weighted graph, sorts the edges in parallel, and then
 * selects the MST edges using the union-find method. Outputs the total cost,
 * MST edges, and execution time.
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 * @return int EXIT_SUCCESS on successful execution.
 */
int main(int argc, char** argv) {
    int num_thread = MAX_THREADS;
    clock_t t = clock();

    // Prompt for input: number of nodes in the graph.
    printf("Enter the number of nodes: ");
    scanf("%d", &n);

    // Allocate memory for the maximum possible number of edges in a complete graph.
    edges = (edge*) malloc(n * (n + 1) / 2 * sizeof(edge));

    // Generate random weights and populate the edge list.
    for (int i = 0; i < n; i++) {
        printf("{");
        for (int j = 0; j < n; j++) {
            int x;
            // Generate a random weight for the edge.
            x = rand() % 99;
            printf("%d, ", x);

            // Avoid self-loops or duplicate edges in undirected graph.
            if (x == -1)
                continue;
            if (i >= j)
                continue;

            edges[num_edge].u = i;
            edges[num_edge].v = j;
            edges[num_edge].w = x;
            num_edge++;
        }
        printf("}\n");
    }
    // Ensure graph has enough edges to form an MST.
    assert(num_edge >= n - 1);

    // Sort the edges in parallel by weight.
    merge_sort(edges, num_edge, comparison_weight);

    // Initialize union-find parent array so that each node is its own parent.
    par = (int*) malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        par[i] = i;
    }

    long long total_cost = 0;  ///< Total weight of the MST.
    int num_chosen = 0;        ///< Counter for MST edges selected.
    chosen_edges = (edge*) malloc(num_edge * sizeof(edge));

    // Select edges for the MST using Kruskal's algorithm.
    for (int i = 0; i < num_edge; i++) {
        int u = edges[i].u;
        int v = edges[i].v;
        int w = edges[i].w;
        // Merge sets if the edge doesn't form a cycle.
        if (merge_set(u, v)) {
            total_cost += w;
            chosen_edges[num_chosen++] = edges[i];
            // Stop once MST is complete.
            if (num_chosen == n - 1)
                break;
        }
    }
    printf("MST Total Cost: %lld\n", total_cost);

    // Sort chosen edges by node labels for organized output.
    merge_sort(chosen_edges, num_chosen, comparison_node);
    for (int i = 0; i < num_chosen; i++) {
        printf("Edge: %d-%d\n", chosen_edges[i].u, chosen_edges[i].v);
    }

    // Calculate and display execution time (scaled for readability).
    double time_taken = ((double) (clock() - t)) / CLOCKS_PER_SEC;
    printf("EXECUTION TIME: %f s\n", (time_taken / 1000));

    // Free allocated memory.
    free(edges);
    free(chosen_edges);
    free(par);
    return 0;
}