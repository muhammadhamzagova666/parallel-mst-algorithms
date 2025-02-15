// Language: C
// ----------------------------------------------------------------------------
// File: kruskal_s.c
// Description: Implements a sequential version of Kruskal's algorithm to compute the
//              Minimum Spanning Tree (MST) of a weighted graph. The program generates
//              a random graph, sorts the edges by weight, and selects MST edges using
//              a union-find structure to prevent cycles.
// Target Users: Developers and researchers interested in graph algorithms.
// ----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>

/**
 * @brief Structure representing an edge in an undirected graph.
 *
 * Each edge connects two vertices (u and v) with an associated weight (w).
 */
typedef struct edge {
    int u, v, w;  // u: source vertex, v: destination vertex, w: edge weight
} edge;

int n;                      ///< Number of nodes in the graph.
edge* edges = NULL;         ///< Array of all generated edges.
edge* chosen_edges = NULL;  ///< Array to store edges selected for the MST.
int* par = NULL;            ///< Parent array for union-find operations.
int num_edge = 0;           ///< Total number of generated edges.

/**
 * @brief Recursively finds the representative (root) of the set to which vertex u belongs.
 *
 * This function applies path compression to flatten the structure of the union-find tree,
 * resulting in more efficient subsequent queries.
 *
 * @param u The vertex for which the set representative is sought.
 * @return int The set representative of vertex u.
 */
int find_set(int u) {
    if (par[u] == u)
        return u;
    return par[u] = find_set(par[u]);  // Path compression for efficiency.
}

/**
 * @brief Unions the sets containing vertices u and v.
 *
 * The function connects the two disjoint sets by linking the representative of one set 
 * to the representative of the other, effectively merging them. It prevents cycles by
 * ensuring that vertices already connected are not merged again.
 *
 * @param u First vertex.
 * @param v Second vertex.
 * @return int Returns 1 if the merge was successful (i.e., u and v were in separate sets),
 *             or 0 if they were already in the same set.
 */
int merge_set(int u, int v) {
    int pu = find_set(u), pv = find_set(v);
    if (pu == pv)
        return 0;  // u and v are already connected; merging would create a cycle.
    par[pv] = pu;  // Merge the sets by making pu the parent.
    return 1;
}

/**
 * @brief Comparator function for sorting edges by weight with tie-breakers.
 *
 * If two edges have equal weights, their source vertices are compared; if those
 * are equal, their destination vertices are compared.
 *
 * @param a Pointer to the first edge.
 * @param b Pointer to the second edge.
 * @return int Negative if *a should come before *b; positive if after; zero if equal.
 */
int comparison_weight(const void* a, const void* b) {
    edge* x = (edge*)a;
    edge* y = (edge*)b;
    if (x->w == y->w) {
        if (x->u == y->u)
            return x->v - y->v;
        return x->u - y->u;
    }
    return x->w - y->w;
}

/**
 * @brief Comparator function for sorting edges by node order.
 *
 * Orders edges first by the source vertex and then by the destination vertex.
 *
 * @param a Pointer to the first edge.
 * @param b Pointer to the second edge.
 * @return int Negative if *a should come before *b; positive if after; zero if equal.
 */
int comparison_node(const void* a, const void* b) {
    edge* x = (edge*)a;
    edge* y = (edge*)b;
    if (x->u == y->u)
        return x->v - y->v;
    return x->u - y->u;
}

/**
 * @brief Merges two sorted subarrays of edges into a single sorted array.
 *
 * This helper function is used by the custom merge sort algorithm.
 *
 * @param edges Destination array where the merged edges will be stored.
 * @param larr Sorted left subarray.
 * @param nl Number of elements in the left subarray.
 * @param rarr Sorted right subarray.
 * @param nr Number of elements in the right subarray.
 * @param comparison Pointer to a comparator function.
 */
void merge(edge edges[], edge larr[], int nl, edge rarr[], int nr, 
           int (*comparison)(const void*, const void*)) {
    int il = 0, ir = 0, j = 0;
    // Merge elements from both subarrays in order.
    while (il < nl && ir < nr) {
        if (comparison(&larr[il], &rarr[ir]) < 0) {
            edges[j++] = larr[il++];
        } else {
            edges[j++] = rarr[ir++];
        }
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
 * @brief Recursively sorts an array of edges using merge sort.
 *
 * This custom merge sort implementation splits the edge array into halves,
 * sorts each half recursively, and then merges them.
 *
 * @param edges Array of edges to sort.
 * @param n Number of elements in the array.
 * @param comparison Comparator function to determine order of edges.
 */
void merge_sort(edge edges[], int n, int (*comparison)(const void*, const void*)) {
    if (n > 1) {
        int m = n / 2;
        edge* larr = malloc(m * sizeof(edge));
        edge* rarr = malloc((n - m) * sizeof(edge));
        memcpy(larr, edges, m * sizeof(edge));
        memcpy(rarr, edges + m, (n - m) * sizeof(edge));

        // Recursively sort the left and right halves.
        merge_sort(larr, m, comparison);
        merge_sort(rarr, n - m, comparison);

        // Merge the two sorted halves into the original array.
        merge(edges, larr, m, rarr, n - m, comparison);
        free(larr);
        free(rarr);
    }
}

/**
 * @brief Entry point of the program.
 *
 * Generates a random weighted graph, computes its MST using Kruskal's algorithm,
 * and outputs the total cost of the MST along with the selected edges.
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 * @return int Returns 0 on successful execution.
 */
int main(int argc, char** argv) {
    clock_t t = clock();
    printf("Enter the number of nodes: ");
    scanf("%d", &n);

    // Allocate memory for the maximum number of edges in a complete graph.
    edges = (edge*)malloc(n * (n + 1) / 2 * sizeof(edge));

    // Generate a random graph: assign random weights between 0 and 98.
    for (int i = 0; i < n; i++) {
        printf("{");
        for (int j = 0; j < n; j++) {
            int x = rand() % 99;
            printf("%d, ", x);
            // Skip self-loops and duplicate edges (only consider edges where i < j).
            if (x == -1) continue;
            if (i >= j) continue;
            edges[num_edge].u = i;
            edges[num_edge].v = j;
            edges[num_edge].w = x;
            num_edge++;
        }
        printf("}\n");
    }
    // Ensure that the graph has enough edges to form a spanning tree.
    assert(num_edge >= n - 1);

    // Sort edges by ascending weight using the standard library quicksort.
    qsort(edges, num_edge, sizeof(edge), comparison_weight);

    // Initialize disjoint-set (union-find): each node is initially in its own set.
    par = (int*)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        par[i] = i;
    }

    long long total_cost = 0;  // Accumulate total weight of the MST.
    int num_chosen = 0;        // Count of edges included in the MST.
    chosen_edges = (edge*)malloc(num_edge * sizeof(edge));

    // Iterate through sorted edges and select the ones that do not form a cycle.
    for (int i = 0; i < num_edge; i++) {
        int u = edges[i].u;
        int v = edges[i].v;
        int w = edges[i].w;
        if (merge_set(u, v)) {  // If adding this edge does not form a cycle.
            total_cost += w;
            chosen_edges[num_chosen++] = edges[i];
            // Stop once the MST contains (n - 1) edges.
            if (num_chosen == n - 1)
                break;
        }
    }

    // Output the total cost of the MST.
    printf("MST Total Cost: %lld\n", total_cost);

    // Sort the chosen MST edges by node order for clear presentation.
    qsort(chosen_edges, num_chosen, sizeof(edge), comparison_node);
    for (int i = 0; i < num_chosen; i++) {
        printf("Edge: %d-%d\n", chosen_edges[i].u, chosen_edges[i].v);
    }

    double time_taken = ((double)(clock() - t)) / CLOCKS_PER_SEC;
    printf("EXECUTION TIME: %f ms\n", time_taken);

    // Free allocated memory.
    free(edges);
    free(par);
    free(chosen_edges);

    return 0;
}