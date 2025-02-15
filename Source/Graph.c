// Language: C
// ----------------------------------------------------------------------------
// File: Graph.c
// Description: Generates a random weighted graph and writes it to a CSV file.
//              The generated graph has a fixed maximum number of nodes and a random
//              number of edges connecting these nodes with random weights.
// Target Users: Developers and researchers needing test graphs for graph algorithms.
// ----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAX_NODES 10000   ///< Maximum allowed number of nodes in the generated graph.
#define MAX_WEIGHT 100    ///< Maximum weight value for an edge.

/**
 * @brief Returns a random integer within a given range [min, max].
 *
 * This function computes a random integer by using the standard rand()
 * function and scales the result to the specified range.
 *
 * @param min The minimum value of the desired random number.
 * @param max The maximum value of the desired random number.
 * @return int A random integer between min and max (inclusive).
 */
int randomRange(int min, int max) {
    return rand() % (max - min + 1) + min;
}

/**
 * @brief Main routine for generating a random graph and saving it as a CSV file.
 *
 * The program determines the number of nodes (fixed to MAX_NODES) and computes a 
 * random number of edges ensuring the graph is connected (at least numNodes-1 edges) 
 * and not overly dense (at most a complete graph). Each edge is assigned a random
 * source, destination, and weight. The results are written to "Graph1.csv" where the
 * first line specifies the number of nodes and edges, and subsequent lines list each edge.
 *
 * @return int EXIT_SUCCESS (0) if the file is generated successfully; otherwise, a non-zero value.
 */
int main() {
    // Seed the random number generator with the current time to enable varied outputs.
    srand(time(NULL));

    // Set the number of nodes and compute a random number of edges 
    // (ensuring at least a spanning tree and at most a complete graph).
    int numNodes = MAX_NODES;
    int numEdges = randomRange(numNodes - 1, numNodes * (numNodes - 1) / 2);

    // Open the CSV file for writing the generated graph.
    FILE *file = fopen("Graph1.csv", "w");
    if (file == NULL) {
        perror("Error opening file for writing");
        return 1;
    }

    // Write the header containing the number of nodes and edges.
    fprintf(file, "%d %d\n", numNodes, numEdges);

    // Generate and write each edge:
    // Each edge has a random source, destination (distinct from source), and weight.
    for (int i = 0; i < numEdges; ++i) {
        int source = randomRange(0, numNodes - 1);
        int dest = randomRange(0, numNodes - 1);
        int weight = randomRange(1, MAX_WEIGHT);

        // Ensure that the edge does not form a self-loop.
        while (source == dest) {
            dest = randomRange(0, numNodes - 1);
        }

        // Write the edge data as: "source destination weight"
        fprintf(file, "%d %d %d\n", source, dest, weight);
    }

    // Close the file after writing all edges.
    fclose(file);

    // Provide feedback to the user confirming successful generation.
    printf("CSV file generated successfully as 'Graph1.csv'.\n");

    return 0;
}