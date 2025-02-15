// Language: C
// ----------------------------------------------------------------------------
// File: boruvka_s.c
// Description: Implements the sequential version of Boruvka's algorithm for
//              computing the Minimum Spanning Tree (MST) of a weighted graph.
//              Reads a graph file, processes the MST using optimal union-find 
//              operations, and outputs the MST along with execution metrics.
// Target Users: Developers and researchers in graph processing and algorithm
//               optimization.
// ----------------------------------------------------------------------------

#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Sentinel value for uninitialized parent in the union-find structure.
const int UNSET_ELEMENT = -1;

/**
 * @brief Structure representing a disjoint-set (union-find) for cycle detection.
 *
 * This data structure uses path compression and union-by-rank to efficiently
 * determine connectivity between graph vertices during MST formation.
 */
typedef struct Set {
    int elements;             ///< Number of elements (nodes) in the set.
    int* canonicalElements;   ///< Parent pointers for each node.
    int* rank;                ///< Rank array for union-by-rank optimization.
} Set;

/**
 * @brief Structure representing a weighted graph.
 *
 * The graph is defined by its vertices and edges. The edge list contains triplets 
 * (source, destination, weight) stored sequentially.
 */
typedef struct WeightedGraph {
    int edges;     ///< Total number of edges.
    int vertices;  ///< Total number of vertices.
    int* edgeList; ///< Array of edges, each represented by 3 consecutive integers.
} WeightedGraph;

/**
 * @brief Allocates memory for a weighted graph's edge list.
 *
 * @param graph Pointer to the WeightedGraph structure.
 * @param vertices Number of vertices in the graph.
 * @param edges Total number of edges in the graph.
 */
void newWeightedGraph(WeightedGraph* graph, const int vertices, const int edges) {
    graph->edges = edges;
    graph->vertices = vertices;
    // Allocate contiguous memory for storing edge triplets.
    graph->edgeList = (int*)calloc(edges * 3, sizeof(int));
}

/**
 * @brief Reads graph information from a file and populates the weighted graph.
 *
 * The input file is expected to have the first line with vertices and edges count,
 * followed by each line listing an edge as "from to weight".
 *
 * @param graph Pointer to the WeightedGraph structure.
 * @param inputFileName Name of the input file containing graph data.
 */
void readGraphFile(WeightedGraph* graph, const char inputFileName[]) {
    FILE* inputFile;
    const char* inputMode = "r";
    inputFile = fopen(inputFileName, inputMode);
    if (inputFile == NULL) {
        fprintf(stderr, "Error: Couldn't open input file '%s', exiting!\n", inputFileName);
        exit(EXIT_FAILURE);
    }

    int fscanfResult;
    int vertices = 0, edges = 0;
    // Read the graph dimensions.
    fscanfResult = fscanf(inputFile, "%d %d", &vertices, &edges);
    newWeightedGraph(graph, vertices, edges);

    int from, to, weight;
    // Read each edge and store in the edge list.
    for (int i = 0; i < edges; i++) {
        fscanfResult = fscanf(inputFile, "%d %d %d", &from, &to, &weight);
        graph->edgeList[i * 3]     = from;
        graph->edgeList[i * 3 + 1] = to;
        graph->edgeList[i * 3 + 2] = weight;

        // Check for unexpected EOF or read error.
        if (fscanfResult == EOF) {
            fprintf(stderr, "Error: Incomplete graph data in input file, exiting!\n");
            fclose(inputFile);
            exit(EXIT_FAILURE);
        }
    }
    fclose(inputFile);
}

/**
 * @brief Releases memory allocated for a disjoint-set.
 *
 * @param set Pointer to the Set structure.
 */
void deleteSet(Set* set) {
    free(set->canonicalElements);
    free(set->rank);
}

/**
 * @brief Releases memory allocated for the weighted graph's edge list.
 *
 * @param graph Pointer to the WeightedGraph structure.
 */
void deleteWeightedGraph(WeightedGraph* graph) {
    free(graph->edgeList);
}

/**
 * @brief Displays the weighted graph's edge list.
 *
 * Useful for debugging and verifying the input graph structure.
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
 * @brief Initializes the union-find disjoint-set structure.
 *
 * Allocates necessary arrays and marks all nodes as uninitialized.
 *
 * @param set Pointer to the Set structure.
 * @param elements Number of nodes in the graph.
 */
void newSet(Set* set, const int elements) {
    set->elements = elements;
    set->canonicalElements = (int*)malloc(elements * sizeof(int));
    // Initialize all parent pointers with UNSET_ELEMENT.
    memset(set->canonicalElements, UNSET_ELEMENT, elements * sizeof(int));
    set->rank = (int*)calloc(elements, sizeof(int));
}

/**
 * @brief Finds the representative set (root) for a given vertex.
 *
 * Implements path compression to flatten the union-find tree structure.
 *
 * @param set Pointer to the Set structure.
 * @param vertex Vertex to determine the representative for.
 * @return Representative element of the vertex's set.
 */
int findSet(const Set* set, const int vertex) {
    if (set->canonicalElements[vertex] == UNSET_ELEMENT) {
        return vertex;
    } else {
        // Recursively update the canonical element for path compression.
        set->canonicalElements[vertex] = findSet(set, set->canonicalElements[vertex]);
        return set->canonicalElements[vertex];
    }
}

/**
 * @brief Merges two subsets in the union-find structure via union by rank.
 *
 * Avoids cycles by merging tree with lower rank into one with higher rank.
 *
 * @param set Pointer to the Set structure.
 * @param parent1 A vertex in the first component.
 * @param parent2 A vertex in the second component.
 */
void unionSet(Set* set, const int parent1, const int parent2) {
    int root1 = findSet(set, parent1);
    int root2 = findSet(set, parent2);

    if (root1 == root2) {
        return; // Vertices already in the same set.
    } else if (set->rank[root1] < set->rank[root2]) {
        set->canonicalElements[root1] = root2;
    } else if (set->rank[root1] > set->rank[root2]) {
        set->canonicalElements[root2] = root1;
    } else {
        // Equal ranks: choose one root arbitrarily and increment its rank.
        set->canonicalElements[root1] = root2;
        set->rank[root2] = set->rank[root1] + 1;
    }
}

/**
 * @brief Copies an edge (triplet: source, destination, weight) from one array to another.
 *
 * @param to Destination array where the edge is copied.
 * @param from Source edge array.
 */
void copyEdge(int* to, int* from) {
    memcpy(to, from, 3 * sizeof(int));
}

/**
 * @brief Computes the Minimum Spanning Tree (MST) using Boruvka's algorithm.
 *
 * Iteratively finds the minimum valid edge for each component and merges them,
 * ensuring no cycles form during processing.
 *
 * @param graph Pointer to the input weighted graph.
 * @param mst Pointer to the weighted graph structure to store the MST.
 */
void mstBoruvka(const WeightedGraph* graph, WeightedGraph* mst) {
    Set set;
    newSet(&set, graph->vertices);  // Initialize union-find for all vertices.

    int edgesMST = 0;
    // Allocate an array to temporarily hold the minimum edge for each component.
    int* closestEdge = (int*)malloc(graph->vertices * 3 * sizeof(int));

    // Iterate over potential merging steps until MST is formed.
    for (int i = 1; i < graph->vertices && edgesMST < graph->vertices - 1; i *= 2) {
        // Reset closest edges by marking weight as maximum.
        for (int j = 0; j < graph->vertices; j++) {
            closestEdge[j * 3 + 2] = INT_MAX;
        }

        // Scan each edge to update the candidate minimum edge for relevant components.
        for (int j = 0; j < graph->edges; j++) {
            int* currentEdge = &graph->edgeList[j * 3];
            int canonicalElements[2] = {findSet(&set, currentEdge[0]), 
                                          findSet(&set, currentEdge[1])};

            // Process edges linking distinct components.
            if (canonicalElements[0] != canonicalElements[1]) {
                for (int k = 0; k < 2; k++) {
                    bool closestEdgeNotSet = (closestEdge[canonicalElements[k] * 3 + 2] == INT_MAX);
                    bool weightSmaller = (currentEdge[2] < closestEdge[canonicalElements[k] * 3 + 2]);
                    if (closestEdgeNotSet || weightSmaller) {
                        // Update candidate edge for the component.
                        copyEdge(&closestEdge[canonicalElements[k] * 3], currentEdge);
                    }
                }
            }
        }

        // Add valid candidate edges to the MST while merging components.
        for (int j = 0; j < graph->vertices; j++) {
            if (closestEdge[j * 3 + 2] != INT_MAX) {
                int from = closestEdge[j * 3];
                int to = closestEdge[j * 3 + 1];

                // Ensure the selected edge connects two different components.
                if (findSet(&set, from) != findSet(&set, to)) {
                    copyEdge(&mst->edgeList[edgesMST * 3], &closestEdge[j * 3]);
                    edgesMST++;  // Increase MST edge count.
                    unionSet(&set, from, to);  // Merge the components.
                }
            }
        }
    }

    // Free temporary structures.
    deleteSet(&set);
    free(closestEdge);
}

/* ---------------------------------------------------------------------------
 * Main Execution Entry Point
 * ---------------------------------------------------------------------------
 */
/**
 * @brief Main function to initialize, compute, and display the MST.
 *
 * Reads a graph from a file, runs Boruvka's algorithm to calculate the MST,
 * and prints metrics such as total weight and execution time.
 *
 * @param argc Number of input arguments.
 * @param argv Array containing the input file name as its first argument.
 * @return int EXIT_SUCCESS on completion.
 */
int main(int argc, char* argv[]) {
    // Define and initialize the graph structures
    WeightedGraph graph;
    newWeightedGraph(&graph, 0, 0);
    WeightedGraph mst;
    newWeightedGraph(&mst, 0, 0);

    // Read the graph from the specified file
    readGraphFile(&graph, argv[1]);
    printf("Original Graph loaded successfully.\n");

    // Allocate MST container: MST will have (vertices - 1) edges.
    newWeightedGraph(&mst, graph.vertices, graph.vertices - 1);

    // Measure execution time for MST computation
    clock_t start = clock();
    mstBoruvka(&graph, &mst);
    clock_t end = clock();

    printf("Minimum Spanning Tree (Boruvka):\n");
    // Uncomment next line to print full MST details:
    //printWeightedGraph(&mst);

    unsigned long weightMST = 0;
    // Sum the weights of all edges in the MST
    for (int i = 0; i < mst.edges; i++) {
        weightMST += mst.edgeList[i * 3 + 2];
    }

    // Display MST metrics
    printf("MST weight: %lu\n", weightMST);
    printf("Time elapsed: %f s\n", ((double)(end - start)) / CLOCKS_PER_SEC);

    // Free allocated memory for graphs
    deleteWeightedGraph(&graph);
    deleteWeightedGraph(&mst);

    return EXIT_SUCCESS;
}