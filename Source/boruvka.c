// Language: C
// ----------------------------------------------------------------------------
// File: boruvka.c
// Description: Implements Boruvka's algorithm using MPI to compute the Minimum
//              Spanning Tree (MST) of a weighted graph. The code distributes the
//              edge list among processes, applies union-find with path compression
//              and union-by-rank heuristics, and collects performance metrics.
// Target Users: Developers and researchers in parallel graph algorithms.
// ----------------------------------------------------------------------------

#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

// Constant for an uninitialized parent element in the union-find structure.
const int UNSET_ELEMENT = -1;
// Global variable to track cumulative inter-process communication time.
double commtime = 0;

/**
 * @brief Disjoint-set (union-find) structure for cycle detection.
 */
typedef struct Set {
    int elements;           ///< Total number of elements.
    int* canonicalElements; ///< Array holding parent pointers.
    int* rank;              ///< Array for union-by-rank heuristic.
} Set;

/**
 * @brief Structure representing a weighted graph.
 *
 * The graph is defined by its number of vertices, number of edges, and a sequential
 * array for edge triplets (source, destination, weight).
 */
typedef struct WeightedGraph {
    int edges;     ///< Total number of edges.
    int vertices;  ///< Total number of vertices.
    int* edgeList; ///< Edge data stored as consecutive triplets.
} WeightedGraph;

/**
 * @brief Allocates and initializes a weighted graph's edge list.
 *
 * @param graph    Pointer to the graph structure.
 * @param vertices Number of vertices.
 * @param edges    Number of edges.
 */
void newWeightedGraph(WeightedGraph* graph, const int vertices, const int edges) {
    graph->edges = edges;
    graph->vertices = vertices;
    // Allocate continuous memory for all edge triplets.
    graph->edgeList = (int*)calloc(edges * 3, sizeof(int));
}

/**
 * @brief Reads graph data from the specified file and populates the weighted graph.
 *
 * The file must include the number of vertices and edges on the first line, followed
 * by a list of edges formatted as "from to weight" for each subsequent line.
 *
 * @param graph         Pointer to the graph structure.
 * @param inputFileName Name of the file containing the graph.
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
    int vertices = 0;
    int edges = 0;
    // Read dimensions: vertices and edges count.
    fscanfResult = fscanf(inputFile, "%d %d", &vertices, &edges);
    newWeightedGraph(graph, vertices, edges);

    int from, to, weight;
    for (int i = 0; i < edges; i++) {
        fscanfResult = fscanf(inputFile, "%d %d %d", &from, &to, &weight);
        graph->edgeList[i * 3]     = from;
        graph->edgeList[i * 3 + 1] = to;
        graph->edgeList[i * 3 + 2] = weight;
        // Check if EOF reached unexpectedly.
        if (fscanfResult == EOF) {
            fprintf(stderr, "Error: Incomplete graph data in file, exiting!\n");
            fclose(inputFile);
            exit(EXIT_FAILURE);
        }
    }
    fclose(inputFile);
}

/**
 * @brief Prints the graph's edge list in "from, to, weight" format.
 *
 * Useful for debugging and verifying file reading.
 *
 * @param graph Pointer to the weighted graph structure.
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
 * @brief Initializes a disjoint-set with the specified number of elements.
 *
 * Sets all canonical pointers to UNSET_ELEMENT and allocates the rank array.
 *
 * @param set      Pointer to the Set structure.
 * @param elements Number of elements (vertices).
 */
void newSet(Set* set, const int elements) {
    set->elements = elements;
    set->canonicalElements = (int*)malloc(elements * sizeof(int));
    // Mark each element as uninitialized.
    memset(set->canonicalElements, UNSET_ELEMENT, elements * sizeof(int));
    set->rank = (int*)calloc(elements, sizeof(int));
}

/**
 * @brief Finds the representative of a given vertex with path compression.
 *
 * Reduces future search times by linking nodes directly to the representative.
 *
 * @param set    Pointer to the Set structure.
 * @param vertex The vertex for which to find the representative.
 * @return int   The representative element.
 */
int findSet(const Set* set, const int vertex) {
    if (set->canonicalElements[vertex] == UNSET_ELEMENT) {
        return vertex;
    } else {
        // Recursive call with path compression.
        set->canonicalElements[vertex] = findSet(set, set->canonicalElements[vertex]);
        return set->canonicalElements[vertex];
    }
}

/**
 * @brief Merges two subsets using union-by-rank.
 *
 * Reduces tree height by attaching the tree with smaller rank under the larger one.
 *
 * @param set     Pointer to the Set structure.
 * @param parent1 A vertex in the first subset.
 * @param parent2 A vertex in the second subset.
 */
void unionSet(Set* set, const int parent1, const int parent2) {
    int root1 = findSet(set, parent1);
    int root2 = findSet(set, parent2);

    if (root1 == root2) {
        return;  // Already in the same subset; no merge needed.
    }
    // Attach smaller tree under larger tree.
    else if (set->rank[root1] < set->rank[root2]) {
        set->canonicalElements[root1] = root2;
    } else if (set->rank[root1] > set->rank[root2]) {
        set->canonicalElements[root2] = root1;
    } else {
        // Equal ranks: choose one as new root and increment its rank.
        set->canonicalElements[root1] = root2;
        set->rank[root2] = set->rank[root1] + 1;
    }
}

/**
 * @brief Copies an edge (triplet: source, destination, weight) from one array to another.
 *
 * @param to   Destination array.
 * @param from Source array holding the edge.
 */
void copyEdge(int* to, int* from) {
    memcpy(to, from, 3 * sizeof(int));
}

/**
 * @brief Distributes the edge list across processes via MPI_Scatter.
 *
 * This function divides the edge list evenly among processes to parallelize
 * the search for the closest edge.
 *
 * @param edgeList      Full edge list in the graph.
 * @param edgeListPart  Buffer for storing the portion of the edge list.
 * @param elements      Total number of edges.
 * @param elementsPart  Pointer to number of edges per process.
 */
void scatterEdgeList(int* edgeList, int* edgeListPart, const int elements, int* elementsPart) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double mpi_start_time = MPI_Wtime();
    MPI_Scatter(edgeList, *elementsPart * 3, MPI_INT,
                edgeListPart, *elementsPart * 3, MPI_INT,
                0, MPI_COMM_WORLD);
    double mpi_end_time = MPI_Wtime();
    commtime += (mpi_end_time - mpi_start_time);

    // Adjust for remainder if total edges not divisible by process count.
    if (rank == size - 1 && elements % *elementsPart != 0) {
        *elementsPart = elements % *elementsPart;
    }

    // Check if the size and process configuration is unsupported.
    if (elements / 2 + 1 < size && elements != size) {
        if (rank == 0) {
            fprintf(stderr, "Error: Unsupported process/edge distribution, exiting!\n");
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}

/**
 * @brief Releases memory allocated for the disjoint-set.
 *
 * @param set Pointer to the Set structure.
 */
void deleteSet(Set* set) {
    free(set->canonicalElements);
    free(set->rank);
}

/**
 * @brief Frees the memory allocated for the weighted graph's edge list.
 *
 * @param graph Pointer to the graph structure.
 */
void deleteWeightedGraph(WeightedGraph* graph) {
    free(graph->edgeList);
}

/**
 * @brief Computes the MST using Boruvka's algorithm in a parallel MPI environment.
 *
 * The algorithm iteratively finds the lightest edge for each component and merges
 * connected components until the MST is formed. During this process, inter-process
 * communication is used to combine candidate edges.
 *
 * @param graph Pointer to the input graph.
 * @param mst   Pointer to the graph structure that will store the MST.
 */
void mstBoruvka(const WeightedGraph* graph, WeightedGraph* mst) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;
    bool parallel = (size != 1);

    // Process 0 broadcasts the number of edges and vertices to all processes.
    int edges, vertices;
    if (rank == 0) {
        edges = graph->edges;
        vertices = graph->vertices;
        double mpi_start_time = MPI_Wtime();
        MPI_Bcast(&edges, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&vertices, 1, MPI_INT, 0, MPI_COMM_WORLD);
        double mpi_end_time = MPI_Wtime();
        commtime += (mpi_end_time - mpi_start_time);
    } else {
        double mpi_start_time = MPI_Wtime();
        MPI_Bcast(&edges, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&vertices, 1, MPI_INT, 0, MPI_COMM_WORLD);
        double mpi_end_time = MPI_Wtime();
        commtime += (mpi_end_time - mpi_start_time);
    }

    // Distribute edge list among processes for parallel search.
    int edgesPart = (edges + size - 1) / size;
    int* edgeListPart = (int*)malloc(edgesPart * 3 * sizeof(int));
    if (parallel) {
        scatterEdgeList(graph->edgeList, edgeListPart, edges, &edgesPart);
    } else {
        edgeListPart = graph->edgeList;
    }

    // Initialize union–find structure for MST component management.
    Set* set = &(Set){ .elements = 0, .canonicalElements = NULL, .rank = NULL };
    newSet(set, vertices);

    int edgesMST = 0;
    // Array to store the best (minimum weight) edge for each component.
    int* closestEdge = (int*)malloc(vertices * 3 * sizeof(int));
    int* closestEdgeRecieved;
    if (parallel) {
        closestEdgeRecieved = (int*)malloc(vertices * 3 * sizeof(int));
    }

    // Iteratively merge components until MST is complete.
    for (int iter = 1; iter < vertices && edgesMST < vertices - 1; iter *= 2) {
        // Reset candidate edge for each vertex.
        for (int j = 0; j < vertices; j++) {
            closestEdge[j * 3 + 2] = INT_MAX;
        }

        // Process local edge list: determine candidate edge for merging.
        for (int j = 0; j < edgesPart; j++) {
            int* currentEdge = &edgeListPart[j * 3];
            int canonicalElements[2] = { findSet(set, currentEdge[0]),
                                           findSet(set, currentEdge[1]) };
            // Only consider edges connecting distinct components.
            if (canonicalElements[0] != canonicalElements[1]) {
                for (int k = 0; k < 2; k++) {
                    bool edgeNotSet = (closestEdge[canonicalElements[k] * 3 + 2] == INT_MAX);
                    bool weightSmaller = (currentEdge[2] < closestEdge[canonicalElements[k] * 3 + 2]);
                    if (edgeNotSet || weightSmaller) {
                        copyEdge(&closestEdge[canonicalElements[k] * 3], currentEdge);
                    }
                }
            }
        }

        // For parallel execution, merge candidate edge results from all processes.
        if (parallel) {
            int from, to;
            for (int step = 1; step < size; step *= 2) {
                if (rank % (2 * step) == 0) {
                    from = rank + step;
                    if (from < size) {
                        double mpi_start_time = MPI_Wtime();
                        MPI_Recv(closestEdgeRecieved, vertices * 3, MPI_INT, from, 0, MPI_COMM_WORLD, &status);
                        double mpi_end_time = MPI_Wtime();
                        commtime += (mpi_end_time - mpi_start_time);
                        // Combine received candidate edges with current results.
                        for (int i = 0; i < vertices; i++) {
                            int index = i * 3;
                            if (closestEdgeRecieved[index + 2] < closestEdge[index + 2]) {
                                copyEdge(&closestEdge[index], &closestEdgeRecieved[index]);
                            }
                        }
                    }
                } else if (rank % step == 0) {
                    to = rank - step;
                    double mpi_start_time = MPI_Wtime();
                    MPI_Send(closestEdge, vertices * 3, MPI_INT, to, 0, MPI_COMM_WORLD);
                    double mpi_end_time = MPI_Wtime();
                    commtime += (mpi_end_time - mpi_start_time);
                }
            }
            // Broadcast the merged candidate edges to all processes.
            double mpi_start_time = MPI_Wtime();
            MPI_Bcast(closestEdge, vertices * 3, MPI_INT, 0, MPI_COMM_WORLD);
            double mpi_end_time = MPI_Wtime();
            commtime += (mpi_end_time - mpi_start_time);
        }

        // Add the selected candidate edges to the MST and merge the corresponding sets.
        for (int j = 0; j < vertices; j++) {
            if (closestEdge[j * 3 + 2] != INT_MAX) {
                int from = closestEdge[j * 3];
                int to = closestEdge[j * 3 + 1];
                // Ensure the edge connects two distinct sets before adding.
                if (findSet(set, from) != findSet(set, to)) {
                    if (rank == 0) {
                        copyEdge(&mst->edgeList[edgesMST * 3], &closestEdge[j * 3]);
                    }
                    edgesMST++;
                    unionSet(set, from, to);
                }
            }
        }
    }

    // Clean up allocated memory.
    deleteSet(set);
    free(closestEdge);
    if (parallel) {
        free(closestEdgeRecieved);
        free(edgeListPart);
    }
}

/**
 * @brief Main entry point for the MPI-enabled Boruvka MST computation.
 *
 * Initializes MPI, reads graph data, runs Boruvka's algorithm, and prints
 * MST metrics including total weight, execution time, and communication overhead.
 *
 * @param argc Number of command line arguments.
 * @param argv Array of command line arguments (expects input file as argv[1]).
 * @return int EXIT_SUCCESS upon completion.
 */
int main(int argc, char* argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Initialize graph and MST containers.
    WeightedGraph* graph = &(WeightedGraph){ .edges = 0, .vertices = 0, .edgeList = NULL };
    WeightedGraph* mst = &(WeightedGraph){ .edges = 0, .vertices = 0, .edgeList = NULL };

    // Process 0 handles file I/O.
    if (rank == 0) {
        readGraphFile(graph, argv[1]);
        printf("Original Graph successfully loaded.\n");
        // Allocate MST container: MST has (vertices - 1) edges.
        newWeightedGraph(mst, graph->vertices, graph->vertices - 1);
    }

    double start = MPI_Wtime();
    mstBoruvka(graph, mst);

    if (rank == 0) {
        // Uncomment the next line to display full MST details.
        // printWeightedGraph(mst);
        unsigned long weightMST = 0;
        // Compute the total weight of the MST.
        for (int i = 0; i < mst->edges; i++) {
            weightMST += mst->edgeList[i * 3 + 2];
        }
        printf("Minimum Spanning Tree (Boruvka):\n");
        printf("MST weight: %lu\n", weightMST);
        printf("Execution Time elapsed: %f s\n", MPI_Wtime() - start);
        printf("Communication Time: %f s\n", commtime);
        // Cleanup graph memory.
        deleteWeightedGraph(graph);
        deleteWeightedGraph(mst);
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}