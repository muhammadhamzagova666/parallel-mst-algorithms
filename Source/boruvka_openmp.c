// Language: C
// ----------------------------------------------------------------------------
// File: boruvka_openmp.c
// Description: Implements a parallel version of Boruvka's algorithm for finding
//              the Minimum Spanning Tree (MST) of a weighted graph using OpenMP.
//              This implementation reads in a graph file, computes the MST using
//              efficient parallel techniques, and outputs the final MST along with
//              performance metrics such as execution and communication time.
// Target Users: Developers and researchers interested in parallel algorithms and
//               graph processing.
// ----------------------------------------------------------------------------

#include <omp.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_THREADS 4

// Constant used to initialize unset parent in union-find structure.
const int UNSET_ELEMENT = -1;
// Global variable to accumulate time spent in inter-thread/process communication.
double commtime = 0;

/**
 * @brief Represents a disjoint-set (union-find) data structure.
 *
 * This structure aids in efficiently managing sets of vertices for cycle
 * detection during MST construction, employing path compression and union by rank.
 */
typedef struct Set {
    int elements;              ///< Number of elements in the set.
    int* canonicalElements;    ///< Array representing the parent of each element.
    int* rank;                 ///< Array representing the rank used in union-by-rank.
} Set;

/**
 * @brief Represents a weighted graph.
 *
 * This structure stores graph information including the number of vertices,
 * the number of edges, and the list of edges. Each edge is represented by a triplet:
 * source, destination, and weight.
 */
typedef struct WeightedGraph {
    int edges;     ///< Total number of edges in the graph.
    int vertices;  ///< Total number of vertices in the graph.
    int* edgeList; ///< Array storing edges in groups of three integers: [from, to, weight].
} WeightedGraph;

/**
 * @brief Initializes a WeightedGraph by allocating memory for its edge list.
 *
 * Memory allocation is performed for a total of (edges * 3) integers, which will
 * store all edge triplets.
 *
 * @param graph    Pointer to the graph to be initialized.
 * @param vertices Number of vertices in the graph.
 * @param edges    Number of edges in the graph.
 */
void newWeightedGraph(WeightedGraph* graph, const int vertices, const int edges) {
    graph->edges = edges;
    graph->vertices = vertices;
    graph->edgeList = (int*) calloc(edges * 3, sizeof(int));
}

/**
 * @brief Reads a graph from a given file and populates the WeightedGraph structure.
 *
 * Expects the graph file to have the first line indicating the number of vertices
 * and edges, followed by each edge represented as three integers: source, destination,
 * and weight.
 *
 * @param graph         Pointer to the graph structure to populate.
 * @param inputFileName Name of the file containing the graph data.
 */
void readGraphFile(WeightedGraph* graph, const char inputFileName[]) {
    // Open the file in read-only mode.
    FILE* inputFile;
    const char* inputMode = "r";
    inputFile = fopen(inputFileName, inputMode);
    if (inputFile == NULL) {
        fprintf(stderr, "Error: Couldn't open input file '%s', exiting!\n", inputFileName);
        exit(EXIT_FAILURE);
    }

    int fscanfResult;
    int vertices = 0, edges = 0;
    // Read graph dimensions from the first line.
    fscanfResult = fscanf(inputFile, "%d %d", &vertices, &edges);
    newWeightedGraph(graph, vertices, edges);

    int from, to, weight;
    // Read each edge and store the triplet in the edge list.
    for (int i = 0; i < edges; i++) {
        fscanfResult = fscanf(inputFile, "%d %d %d", &from, &to, &weight);
        graph->edgeList[i * 3]     = from;
        graph->edgeList[i * 3 + 1] = to;
        graph->edgeList[i * 3 + 2] = weight;

        // If an unexpected EOF is encountered, terminate with error.
        if (fscanfResult == EOF) {
            fprintf(stderr,"Error: Incomplete graph data in file, exiting!\n");
            fclose(inputFile);
            exit(EXIT_FAILURE);
        }
    }
    fclose(inputFile);
}

/**
 * @brief Deallocates memory for the disjoint-set.
 *
 * Frees the memory allocated for both the canonicalElements and rank arrays.
 *
 * @param set Pointer to the Set structure.
 */
void deleteSet(Set* set) {
    free(set->canonicalElements);
    free(set->rank);
}

/**
 * @brief Deallocates memory for the weighted graph.
 *
 * Frees the dynamically allocated edge list.
 *
 * @param graph Pointer to the WeightedGraph structure.
 */
void deleteWeightedGraph(WeightedGraph* graph) {
    free(graph->edgeList);
}

/**
 * @brief Prints the entire weighted graph.
 *
 * Displays all edges in a tabulated "source destination weight" format,
 * aiding in debugging and verification of the input graph.
 *
 * @param graph Pointer to the WeightedGraph structure.
 */
void printWeightedGraph(const WeightedGraph* graph) {
    printf("------------------------------------------------\n");
    for (int i = 0; i < graph->edges; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%d\t", graph->edgeList[i * 3 + j]);
        }
        printf("\n");
    }
    printf("------------------------------------------------\n");
}

/**
 * @brief Initializes a disjoint-set for union-find operations.
 *
 * Allocates memory for the canonicalElements and rank arrays, and initializes
 * the canonicalElements to UNSET_ELEMENT to indicate that they have no parent yet.
 *
 * @param set      Pointer to the Set structure.
 * @param elements Number of elements (vertices) to be managed.
 */
void newSet(Set* set, const int elements) {
    set->elements = elements;
    set->canonicalElements = (int*) malloc(elements * sizeof(int)); // Parent pointers.
    memset(set->canonicalElements, UNSET_ELEMENT, elements * sizeof(int));
    set->rank = (int*) calloc(elements, sizeof(int)); // Rank array for union-by-rank.
}

/**
 * @brief Finds and returns the representative (canonical) element of a given vertex.
 *
 * Includes path compression to flatten the tree structure for subsequent calls.
 *
 * @param set    Const pointer to the Set structure.
 * @param vertex The vertex whose set representative is required.
 * @return The representative element of the vertex.
 */
int findSet(const Set* set, const int vertex) {
    if (set->canonicalElements[vertex] == UNSET_ELEMENT) {
        return vertex;
    } else {
        // Recursively find the root and apply path compression.
        set->canonicalElements[vertex] = findSet(set, set->canonicalElements[vertex]);
        return set->canonicalElements[vertex];
    }
}

/**
 * @brief Merges the sets containing two given vertices using union by rank.
 *
 * This function prevents the formation of cycles when adding edges to the MST.
 *
 * @param set     Pointer to the Set structure.
 * @param parent1 A vertex in the first set.
 * @param parent2 A vertex in the second set.
 */
void unionSet(Set* set, const int parent1, const int parent2) {
    int root1 = findSet(set, parent1);
    int root2 = findSet(set, parent2);

    if (root1 == root2) {
        return; // Already in the same set; no action needed.
    }
    // Attaching the tree with a smaller rank under the tree with a larger rank.
    else if (set->rank[root1] < set->rank[root2]) {
        set->canonicalElements[root1] = root2;
    } else if (set->rank[root1] > set->rank[root2]) {
        set->canonicalElements[root2] = root1;
    } else {
        // If ranks are equal, choose one as new root and increment its rank.
        set->canonicalElements[root1] = root2;
        set->rank[root2] = set->rank[root1] + 1;
    }
}

/**
 * @brief Copies an edge (3 integers) from one array to another.
 *
 * Uses memcpy for efficiency, ensuring exactly 3 integers (from, to, weight)
 * are copied.
 *
 * @param to   Destination array where the edge will be copied.
 * @param from Source edge array.
 */
void copyEdge(int* to, int* from) {
    memcpy(to, from, 3 * sizeof(int));
}

/**
 * @brief Computes the Minimum Spanning Tree (MST) of a graph using Boruvka's algorithm.
 *
 * The algorithm runs in iterations, finding the minimum outgoing edge for each
 * component in parallel using OpenMP. The chosen edges are then added to the MST,
 * and the corresponding components are merged in the union-find structure.
 *
 * @param graph Pointer to the original graph.
 * @param mst   Pointer to the graph structure that will store the MST.
 */
void mstBoruvka(const WeightedGraph* graph, WeightedGraph* mst) {
    // Initialize union-find data structure for components.
    Set* set = &(Set){ .elements = 0, .canonicalElements = NULL, .rank = NULL };
    newSet(set, graph->vertices);

    int edgesMST = 0;
    // Array to store the best (minimum weight) edge for each component.
    int* closestEdge = (int*) malloc(graph->vertices * 3 * sizeof(int));
    
    double start;
    
    // Doubling iterations: each iteration potentially merges components.
    for (int i = 1; i < graph->vertices && edgesMST < graph->vertices - 1; i *= 2) {
        // Reset the closest edge for every component.
        double start = omp_get_wtime();
        #pragma omp parallel for
        for (int j = 0; j < graph->vertices; j++) {
            closestEdge[j * 3 + 2] = INT_MAX;
        }

        // Parallel search: For each edge, update the closest edge for each component.
        #pragma omp parallel for
        for (int j = 0; j < graph->edges; j++) {
            int* currentEdge = &graph->edgeList[j * 3];
            int canonicalElements[2] = { findSet(set, currentEdge[0]), findSet(set, currentEdge[1]) };

            // If edge connects two different components, test if it is the lightest so far.
            if (canonicalElements[0] != canonicalElements[1]) {
                for (int k = 0; k < 2; k++) {
                    bool closestEdgeNotSet = closestEdge[canonicalElements[k] * 3 + 2] == INT_MAX;
                    bool weightSmaller = currentEdge[2] < closestEdge[canonicalElements[k] * 3 + 2];
                    if (closestEdgeNotSet || weightSmaller) {
                        // Critical section to avoid data races when updating shared data.
                        #pragma omp critical
                        {
                            copyEdge(&closestEdge[canonicalElements[k] * 3], currentEdge);
                        }
                    }
                }
            }
        }
        // Accumulate communication time from edge scanning.
        commtime += omp_get_wtime() - start;

        // Add the selected minimum edges to the MST and merge components.
        double start2 = omp_get_wtime();
        #pragma omp parallel for
        for (int j = 0; j < graph->vertices; j++) {
            if (closestEdge[j * 3 + 2] != INT_MAX) {
                int from = closestEdge[j * 3];
                int to = closestEdge[j * 3 + 1];

                // Verify and prevent adding duplicate edges.
                if (findSet(set, from) != findSet(set, to)) {
                    #pragma omp critical
                    {
                        copyEdge(&mst->edgeList[edgesMST * 3], &closestEdge[j * 3]);
                        edgesMST++;
                    }
                    unionSet(set, from, to);
                }
            }
        }
        // Accumulate communication time for component merging.
        start += omp_get_wtime() - start2;
    }
    
    // Clean up allocated memory.
    deleteSet(set);
    free(closestEdge);
}

/**
 * @brief Main entry point of the Boruvka's MST program.
 *
 * Sets up the parallel environment, reads the input graph file, computes the MST,
 * and outputs the MST weight, execution time, and communication overhead.
 *
 * @param argc Number of command line arguments.
 * @param argv Array containing the input file name as its first argument.
 * @return EXIT_SUCCESS if the program completes successfully; otherwise EXIT_FAILURE.
 */
int main(int argc, char* argv[]) {
    omp_set_num_threads(MAX_THREADS);
    
    // Initialize graph and MST structures using designated initializers.
    WeightedGraph* graph = &(WeightedGraph){ .edges = 0, .vertices = 0, .edgeList = NULL };
    WeightedGraph* mst = &(WeightedGraph){ .edges = 0, .vertices = 0, .edgeList = NULL };

    // Read the input file to populate the graph.
    readGraphFile(graph, argv[1]);

    // Optionally print graph details for debugging.
    printf("Original Graph loaded successfully.\n");
    // printWeightedGraph(graph);

    // Prepare the MST graph container; MST will have (vertices - 1) edges.
    newWeightedGraph(mst, graph->vertices, graph->vertices - 1);

    double start = omp_get_wtime();
    // Compute the MST using the parallel Boruvka algorithm.
    mstBoruvka(graph, mst);

    // Optionally print the resulting MST for verification.
    printf("Minimum Spanning Tree (Boruvka):\n");
    // printWeightedGraph(mst);

    unsigned long weightMST = 0;
    // Compute total weight for performance and correctness verification.
    for (int i = 0; i < mst->edges; i++) {
        weightMST += mst->edgeList[i * 3 + 2];
    }

    // Display summary metrics.
    printf("MST weight: %lu\n", weightMST);
    printf("Time elapsed: %f s\n", omp_get_wtime() - start);
    printf("Communication Time: %f s\n", commtime);
    
    // Clean up allocated memory before exiting.
    deleteWeightedGraph(graph);
    deleteWeightedGraph(mst);

    return EXIT_SUCCESS;
}