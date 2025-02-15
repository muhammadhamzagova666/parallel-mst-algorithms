// Language: C
// ----------------------------------------------------------------------------
// File: kruskal_mpi.c
// Description: Implements Kruskal's algorithm using MPI to compute the Minimum
//              Spanning Tree (MST) of a weighted graph in parallel. The algorithm
//              sorts edges using a custom merge sort and utilizes the union-find
//              data structure to avoid cycles. The program distributes edge data
//              among MPI processes and gathers sorted results for MST computation.
// Target Users: Developers and researchers interested in parallel graph algorithms.
// ----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <time.h>
#include <assert.h>

double total_comm_time = 0;  // Global to accumulate total communication time

/**
 * @brief Structure representing a weighted edge in the graph.
 */
typedef struct Edge {
    int u, v, w;  ///< u: source vertex, v: destination vertex, w: weight of the edge
} Edge;

int numVertices;             ///< Number of vertices in the graph
Edge* allEdges;              ///< Array containing all graph edges (only populated in rank 0)
Edge* chosenEdges;           ///< Array storing the MST edges chosen by Kruskal's algorithm
int* parent;                 ///< Array used for union-find: stores parent for each vertex
int numEdges;                ///< Total number of edges in the graph

/**
 * @brief Finds the representative element (set) for a given vertex using path compression.
 *
 * @param u The vertex whose set representative is to be found.
 * @return int The representative element of the set containing u.
 */
int findSet(int u) {
    if (parent[u] == u)
        return u;
    return parent[u] = findSet(parent[u]); // Recursively compress path for faster future queries.
}

/**
 * @brief Merges the sets containing vertices u and v.
 *
 * @param u First vertex.
 * @param v Second vertex.
 * @return int Returns 1 if the merge was successful (i.e., the nodes were in different sets),
 *             or 0 if they were already in the same set.
 */
int mergeSets(int u, int v) {
    int pu = findSet(u), pv = findSet(v);
    if (pu == pv)
        return 0;
    parent[pv] = pu;  // Merge the sets by linking pv's representative to pu.
    return 1;
}

/**
 * @brief Comparator function to compare two edges based on weight and vertex order.
 *
 * This comparator is used for sorting edges in ascending order. If two edges have the same
 * weight, the one with the smaller source vertex (u) is considered smaller; if those are equal,
 * the one with the smaller destination vertex (v) is considered smaller.
 *
 * @param x Pointer to the first edge.
 * @param y Pointer to the second edge.
 * @return int Non-zero if x should come before y; zero otherwise.
 */
int compareWeight(Edge* x, Edge* y) {
    if (x->w == y->w) {
        if (x->u == y->u)
            return x->v < y->v;
        return x->u < y->u;
    }
    return x->w < y->w;
}

/**
 * @brief Merges two sorted subarrays of edges into a single sorted array.
 *
 * The merging operation compares elements from left and right subarrays using the provided
 * comparison function.
 *
 * @param edges Destination array where the merged result will be stored.
 * @param larr Left subarray (sorted).
 * @param nl Number of elements in the left subarray.
 * @param rarr Right subarray (sorted).
 * @param nr Number of elements in the right subarray.
 * @param comparison Function pointer for comparing two Edge elements.
 * @param offset The starting index offset for the destination array.
 */
void merge(Edge edges[], Edge larr[], int nl, Edge rarr[], int nr, int (*comparison)(Edge*, Edge*), int offset) {
    int il = 0, ir = 0, j = offset;
    // Merge until one subarray is exhausted.
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
    // Copy any remaining elements from the left subarray.
    while (il < nl) {
        edges[j++] = larr[il++];
    }
    // Copy any remaining elements from the right subarray.
    while (ir < nr) {
        edges[j++] = rarr[ir++];
    }
}

/**
 * @brief Recursive merge sort implementation for an array of Edge elements.
 *
 * This function uses merge sort to sort the edge array according to the provided comparison 
 * function. It divides the array into halves, sorts each half recursively, and merges them.
 *
 * @param edges The array of Edge elements to sort.
 * @param n The number of elements in the array.
 * @param comparison Function pointer to compare two Edge elements.
 */
void mergeSort(Edge edges[], int n, int (*comparison)(Edge*, Edge*)) {
    if (n > 1) {
        int m = n / 2;
        Edge* larr = malloc(m * sizeof(Edge));
        Edge* rarr = malloc((n - m) * sizeof(Edge));
        memcpy(larr, edges, m * sizeof(Edge));
        memcpy(rarr, edges + m, (n - m) * sizeof(Edge));

        mergeSort(larr, m, comparison);
        mergeSort(rarr, n - m, comparison);

        merge(edges, larr, m, rarr, n - m, comparison, 0);
        free(larr); 
        free(rarr);
    }
}

/**
 * @brief Recursively merges gathered sorted subarrays from different MPI processes.
 *
 * The mergeGather function reduces multiple sorted segments into one fully sorted sequence.
 *
 * @param gather Array containing gathered edges from all processes.
 * @param sendcounts Array holding the number of elements each process contributed.
 * @param displ Array holding the displacement indices in the gathered array.
 * @param comparison Function pointer for Edge comparison.
 * @param l Left index of the current merge segment.
 * @param r Right index of the current merge segment.
 */
void mergeGather(Edge gather[], int sendcounts[], int displ[], int (*comparison)(Edge*, Edge*), int l, int r) {
    if (l >= r)
        return;
    int mid = (l + r) >> 1;
    mergeGather(gather, sendcounts, displ, comparison, l, mid);
    mergeGather(gather, sendcounts, displ, comparison, mid + 1, r);
    
    int nl = 0, nr = 0;
    for (int i = l; i <= r; i++) {
        if (i <= mid)
            nl += sendcounts[i];
        else
            nr += sendcounts[i];
    }
    Edge* larr = malloc(nl * sizeof(Edge));
    Edge* rarr = malloc(nr * sizeof(Edge));
    memcpy(larr, gather + displ[l], nl * sizeof(Edge));
    memcpy(rarr, gather + displ[l] + nl, nr * sizeof(Edge));
    merge(gather, larr, nl, rarr, nr, comparison, displ[l]);
    free(larr);
    free(rarr);
}

int main(int argc, char** argv) {
    int worldSize, worldRank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

    // Create a custom MPI datatype for the Edge structure.
    MPI_Datatype MPI_EDGE;
    int blockCount = 3;
    const int blockLengths[3] = { 1, 1, 1 };
    const MPI_Aint displacements[3] = { 0, sizeof(int), 2 * sizeof(int) };  // Memory layout is contiguous.
    MPI_Datatype blockTypes[3] = { MPI_INT, MPI_INT, MPI_INT };
    MPI_Type_create_struct(blockCount, blockLengths, displacements, blockTypes, &MPI_EDGE);
    MPI_Type_commit(&MPI_EDGE);

    // Start input processing and timing.
    clock_t t = clock();
    if (worldRank == 0) {
        printf("Enter number of nodes: ");
        scanf("%d", &numVertices);
        // Allocate maximum possible edges for a complete graph (upper bound).
        allEdges = (Edge*) malloc(numVertices * (numVertices + 1) / 2 * sizeof(Edge));
        for (int i = 0; i < numVertices; i++) {
            printf("{");
            for (int j = 0; j < numVertices; j++) {
                int x;
                // Generate a random weight for the edge.
                x = rand() % 99;
                printf("%d, ", x);
                // Skip self-loops and duplicate edges in an undirected graph.
                if (x == -1) continue;
                if (i >= j) continue;
                allEdges[numEdges].u = i;
                allEdges[numEdges].v = j;
                allEdges[numEdges].w = x;
                numEdges++;
            }
            printf("}\n");
        }
        // Ensure the graph is connected.
        assert(numEdges >= numVertices - 1);
    }
    
    // Synchronize processes before communication.
    MPI_Barrier(MPI_COMM_WORLD);

    // Broadcast number of edges to all processes.
    double mpi_start_time = MPI_Wtime(); 
    if (worldRank == 0) {
        for (int i = 1; i < worldSize; i++)
            MPI_Send(&numEdges, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    } else {
        MPI_Recv(&numEdges, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    double mpi_end_time = MPI_Wtime();
    printf("Process %d, MPI Send/Recv time: %f seconds\n", worldRank, mpi_end_time - mpi_start_time);
    
    // Process: Distribute edges among MPI processes for parallel sorting.
    {
        // Allocate arrays for gathered and scattered data.
        Edge* gathered = malloc(numEdges * sizeof(Edge));
        int scatteredSize = numEdges / worldSize;
        Edge* scattered = malloc((scatteredSize + 1) * sizeof(Edge));

        // Compute sendCounts and displacements for scatter/gather operations.
        int* sendCounts = malloc(worldSize * sizeof(int));
        int* displs = malloc(worldSize * sizeof(int));
        int rem = numEdges % worldSize;
        int sum = 0;
        for (int i = 0; i < worldSize; i++) {
            sendCounts[i] = numEdges / worldSize + (i < rem);
            displs[i] = sum;
            sum += sendCounts[i];
        }

        // Scatter edges to all processes.
        mpi_start_time = MPI_Wtime();
        MPI_Scatterv(allEdges, sendCounts, displs, MPI_EDGE, scattered, sendCounts[worldRank], MPI_EDGE, 0, MPI_COMM_WORLD);

        // Locally sort the edges using merge sort.
        mergeSort(scattered, sendCounts[worldRank], compareWeight);

        // Gather sorted edges back to root process.
        MPI_Gatherv(scattered, sendCounts[worldRank], MPI_EDGE, gathered, sendCounts, displs, MPI_EDGE, 0, MPI_COMM_WORLD);
        mpi_end_time = MPI_Wtime();
        printf("Process %d, MPI Scatterv/Gatherv time: %f seconds\n", worldRank, mpi_end_time - mpi_start_time);

        if (worldRank == 0) {
            // Copy gathered array back to allEdges and perform a final merge.
            for (int i = 0; i < numEdges; i++) {
                allEdges[i] = gathered[i];
            }
            mergeGather(gathered, sendCounts, displs, compareWeight, 0, worldSize - 1);
            for (int i = 0; i < numEdges; i++) {
                allEdges[i] = gathered[i];
            }
        }
        free(gathered); 
        free(scattered);
        free(sendCounts); 
        free(displs);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Kruskal's algorithm: select edges for the MST.
    int numChosen = 0;
    long long totalCost = 0;
    if (worldRank == 0) {
        // Initialize union-find structure.
        parent = (int*) malloc(numVertices * sizeof(int));
        for (int i = 0; i < numVertices; i++) {
            parent[i] = i;
        }
        chosenEdges = (Edge*) malloc(numEdges * sizeof(Edge));
        // Iterate over sorted edges and choose valid ones.
        for (int i = 0; i < numEdges; i++) {
            int u = allEdges[i].u;
            int v = allEdges[i].v;
            int w = allEdges[i].w;
            if (mergeSets(u, v)) {
                totalCost += w;
                chosenEdges[numChosen++] = allEdges[i];
                // Stop once the MST has (numVertices - 1) edges.
                if (numChosen == numVertices - 1)
                    break;
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Share the number of chosen MST edges across all processes.
    mpi_start_time = MPI_Wtime();
    if (worldRank == 0) {
        for (int i = 1; i < worldSize; i++)
            MPI_Send(&numChosen, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    } else {
        MPI_Recv(&numChosen, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    mpi_end_time = MPI_Wtime();
    printf("Process %d, 2nd MPI Send/Recv time: %f seconds\n", worldRank, mpi_end_time - mpi_start_time);

    // Output the MST results on the root process.
    if (worldRank == 0) {
        printf("MST Total Weight: %lld\n", totalCost);
        for (int i = 0; i < numChosen; i++) {
            printf("Edge: %d-%d\n", chosenEdges[i].u, chosenEdges[i].v);
        }
        double timeTaken = ((double)(clock() - t)) / CLOCKS_PER_SEC;
        printf("Execution Time: %f seconds\n", timeTaken);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Clean up MPI datatype and finalize the MPI environment.
    MPI_Type_free(&MPI_EDGE);
    MPI_Finalize();
    return 0;
}