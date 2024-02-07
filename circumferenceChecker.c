/**
 * circumferenceChecker.c
 * 
 * Author: Jarne Renders (jarne.renders@kuleuven.be)
 * 
 * Circumference of graphs up to and including 6 vertices has been verified by
 * an independent program. As well as snarks up to and including order 36
 * whose circumference is 2 less than their order.
 * 
 * Longest induced path length verified up to and including 10 vertices by an
 * independent program. Counts of graphs without an induced path of length 4
 * correspond with http://oeis.org/A078564 checked up to and including 10
 * vertices. Note that their P5 is a path with 4 edges.
 * 
 * Longest induced cycle length verified up to and including 7 vertices by an
 * independent program. Counts of graphs without induced cycle of length 5
 * correspond with http://oeis.org/A078566 checked up to and including 10
 * vertices. 
 *
 */

#define USAGE "Usage: ./circumferenceChecker [-cf#|-pf#|-l] [-Cdo#] [-h]"

#define HELPTEXT \
"Count and or filter graphs depending on their circumference, length,\n\
induced cycles or paths.\n\
\n\
Graphs are read from stdin in graph6 format. Graphs are sent to stdout in\n\
graph6 format. If the input graph had a graph6 header, so will the\n\
output graph (if it passes through the filter).\n\
\n\
If no options are passed the program will compute the circumference of\n\
the input graphs.\n\
\n\
    -c, --induced-cycle\n\
            count the longest induced cycle of each graph and print in a table.\n\
    -C, --complement\n\
            print all the graphs that would not have been sent to stdout and\n\
            vice versa.\n\
    -d, --difference\n\
            count difference with order of the graph.\n\
    -f#, --forbidden=#\n\
            send all graphs to stdout that contain an induced path or cycle\n\
            of length # depending on the presence of -c or -p. Length of path\n\
            is the number of edges.\n\
    -h, --help\n\
            print help message.\n\
    -l, --length\n\
            find the length of each graph, i.e. the number of edges in a\n\
            longest path, and print in a table.\n\
    -o#, --output=#\n\
            send all graphs with value # in the table to stdout.\n\
    -p, --induced-path\n\
            count the longest induced path of each graph and print in a table.\n"


#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <getopt.h>
#include <time.h>
#include <string.h>
#include "libs/readGraph6.h"
#include "libs/bitset.h"
#include "libs/hamiltonicityMethods.h"

struct graph {
    bitset *adjacencyList;
    int nv;
};

void printGraph(struct graph *g) {
    for(int i = 0; i < g->nv; i++) {
        fprintf(stderr, "%d: ", i);
        forEach(nbr, g->adjacencyList[i]) {
            fprintf(stderr, "%d ", nbr);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

//******************************************************************************
//
//                   Methods for circumference checker
//
//******************************************************************************

bool canBeCycleOfLength(struct graph *g, bitset remainingVertices, int
lastElemOfPath, int firstElemOfPath, int cycleLength, int pathLength) {

    // Check if current path is a cycle of the required length.
    if((pathLength == cycleLength) &&
     contains(g->adjacencyList[firstElemOfPath], lastElemOfPath)) {
        return true;
    }

    // If start of path cannot be closed, path cannot become a cycle.
    if(isEmpty(intersection(g->adjacencyList[firstElemOfPath],
     remainingVertices))) { 
        return false;
    }

    bitset neighboursOfLastNotInPath = 
     intersection(g->adjacencyList[lastElemOfPath], remainingVertices);
    forEach(neighbour, neighboursOfLastNotInPath) {

        int oldElemOfPath = lastElemOfPath;

        removeElement(remainingVertices, neighbour);
        lastElemOfPath = neighbour; // Neighbour is the new last element.

        if (canBeCycleOfLength(g, remainingVertices, lastElemOfPath, 
         firstElemOfPath, cycleLength, pathLength + 1)) {
            return true;
        }

        add(remainingVertices, lastElemOfPath);
        lastElemOfPath = oldElemOfPath;
    }

    return false;
}


// Function to find an included vertex with the lowest degree.
int findLowestDegreeVertex(struct graph *g, bitset includedVertices) {

    // First included vertex.
    int start = next(includedVertices, -1);
    int lowestDegree = size(g->adjacencyList[start]);

    // Find an included vertex of lowest degree.
    forEachAfterIndex(u, includedVertices, start) {
        int degree = size(intersection(g->adjacencyList[u], includedVertices));

        if (lowestDegree > degree) {
            lowestDegree = degree;
            start = u;
        }
    }

    return start;
}

int getCircumference(struct graph *g, bitset excludedVertices) {

    // Check backwards from k = n to 3 if there is a cycle of length k.
    for(int i = g->nv; i > 2; i--) {

        bitset forbiddenVertices = excludedVertices;

        // Repeatedly check for a k-cycle in which the previous starting vertex
        // (v) is also forbidden.
        for(int j = g->nv - i; j >= 0; j--) {

            bitset includedVertices = 
             complement(forbiddenVertices, g->nv);

            if(isEmpty(includedVertices)) return false;

            int v = findLowestDegreeVertex(g, includedVertices); 

            // Loop over included neighbours of start and for each such
            // neighbour loop over the included neighbours of start that are of
            // higher index. 
            forEach(w, intersection(g->adjacencyList[v], includedVertices)) {
                forEachAfterIndex(u, 
                 intersection(g->adjacencyList[v], includedVertices), w) {

                    // Create path uvw. We have u > w, so that we eliminate the
                    // checking of paths which are mirrored.
                    bitset path = singleton(v);
                    add(path, u);
                    add(path, w);
                    bitset remainingVertices = 
                     difference(includedVertices, path);

                    if (canBeCycleOfLength(g, remainingVertices, u, w,
                     i - size(excludedVertices), 3)) {
                        return i;
                    }
                }
            }
            add(forbiddenVertices, v);
        }
    }
    return 0;
}


//******************************************************************************
//
//                          Methods for graph length
//
//******************************************************************************

// // Method 1: join K1 with graph and compute its circumference.

// // Free gNew->adjacencyList after!
// void joinWithK1(struct graph *g, struct graph *gNew) {
//     gNew->nv = g->nv + 1;
//     gNew->adjacencyList = malloc(gNew->nv * sizeof(bitset)); 

//     for(int i = 0; i < g->nv; i++) {
//         gNew->adjacencyList[i] = union(g->adjacencyList[i], singleton(g->nv));
//     }
//     gNew->adjacencyList[g->nv] = complement(EMPTY, g->nv);
// }

// int getLength(struct graph *g) {
//     struct graph g2;
//     joinWithK1(g, &g2); 

//     int length = getCircumference(&g2, EMPTY) - 2;
//     if(length < 0) length = 0;

//     free(g2.adjacencyList);
//     return length;
// }

// Method 2: directly search for a longest path. 

// Paths have an active end to which gets built, hence starting with uv will not
// yield the same paths as starting with vu
void searchLongestSuperPath(struct graph *g, bitset remainingVertices,
 int lastElemOfPath, int firstElemOfPath, int *orderOfLongestPath,
 int orderOfPath) {

    //  If found path of largest possible length, we are done.
    if(*orderOfLongestPath == g->nv) {
        return;
    }

    // Check whether we current path is longest.
    if(orderOfPath > *orderOfLongestPath) {
        *orderOfLongestPath = orderOfPath;
    }
    
    bitset neighboursOfLastNotInPath = 
     intersection(g->adjacencyList[lastElemOfPath], remainingVertices);

    forEach(neighbour, neighboursOfLastNotInPath) {

        int oldElemOfPath = lastElemOfPath;

        removeElement(remainingVertices, oldElemOfPath);

        // Neighbour is the new last element.
        lastElemOfPath = neighbour; 

        searchLongestSuperPath(g, remainingVertices, lastElemOfPath,
         firstElemOfPath, orderOfLongestPath, orderOfPath + 1);

        add(remainingVertices, oldElemOfPath);
        lastElemOfPath = oldElemOfPath;
    }
}

int getLength(struct graph *g) {

    int orderOfLongestPath = 0;

    // For each vertex find a longest path starting with v.
    for(int v = 0; v < g->nv; v++) {

        bitset remainingVertices = complement(singleton(v), g->nv);

        // The number of vertices gets stored in orderOfLongestPath
        forEach(w, g->adjacencyList[v]) {

            removeElement(remainingVertices, w);

            searchLongestSuperPath(g, remainingVertices, w, v,
             &orderOfLongestPath, 2);

            add(remainingVertices, w);
        }
    }

    //  Length of a path is number of edges in it.
    int pathLength = orderOfLongestPath - 1;
    if(pathLength < 0) pathLength = 0;

    return pathLength;
}

//******************************************************************************
//
//                      Longest induced cycle
//
//******************************************************************************

void lengthOfLongestInducedSuperCycle(struct graph *g, bitset remainingVertices,
 int lastElemOfPath, int firstElemOfPath, int *longestCycleLength,
 unsigned long long int numberOfLengths[], int pathLength) {

    // Check whether we have a longest induced path already and whether this
    // path is a cycle.
    if(contains(g->adjacencyList[firstElemOfPath], lastElemOfPath)) {
        if(pathLength > *longestCycleLength) {
            *longestCycleLength = pathLength;
        }

        // Counter number of lengths of induced cycles. Used for keeping track
        // of graphs with forbidden lengths of induced cycles.
        numberOfLengths[pathLength]++;

        //  Need to return because we already have a cycle. 
        return;
    }

    // Check if cycle can still be closed with remaining vertices.
    if(isEmpty(intersection(g->adjacencyList[firstElemOfPath], 
     remainingVertices))) { 
        return;
    }

    bitset neighboursOfLastNotInPath = 
     intersection(g->adjacencyList[lastElemOfPath], remainingVertices);

    forEach(neighbour, neighboursOfLastNotInPath) {

        int oldElemOfPath = lastElemOfPath;

        //  Delete neighbours of oldElemOfPath since no chord are allowed in the
        //  induced cycle.
        remainingVertices = difference(remainingVertices,
         g->adjacencyList[oldElemOfPath]);

        bitset deletedVertices = difference(g->adjacencyList[oldElemOfPath],
         remainingVertices);

        lastElemOfPath = neighbour;

        lengthOfLongestInducedSuperCycle(g, remainingVertices, lastElemOfPath,
         firstElemOfPath, longestCycleLength, numberOfLengths, pathLength + 1);

        //  Once oldElemOfPath is again an end of the path, the neighbours we
        //  removed should be considered again.
        remainingVertices = union(remainingVertices, deletedVertices);
        lastElemOfPath = oldElemOfPath;
    }
}

int getLongestInducedCycleLength(struct graph *g, 
 unsigned long long int numberOfLengths[]) {

    int longestInducedCycleLength = 0;

    // Loop over all vertices v and find the longed induced cycle going through
    // it.
    for(int v = 0; v < g->nv; v++) {

        bitset includedVertices = complement(EMPTY, g->nv);

        // Loop over included neighbours w of v and for each w loop over the
        // included neighbours u of v that are of higher index than w.
        forEachAfterIndex(w, intersection(g->adjacencyList[v], includedVertices),
         v) {
            forEachAfterIndex(u, 
             intersection(g->adjacencyList[v], includedVertices), w) {

                bitset remainingVertices = difference(includedVertices,
                 union(g->adjacencyList[v], singleton(v)));

                // Stores length of longest induced cycle containing uvw in
                // longestInducedCycleLength if it is the largest length
                // encountered.
                lengthOfLongestInducedSuperCycle(g, remainingVertices, u, w,
                 &longestInducedCycleLength, numberOfLengths, 3);
            }
        }
    }

    return longestInducedCycleLength;
}


//******************************************************************************
//
//                          Longest induced path
//
//******************************************************************************

void searchLongestInducedSuperPath(struct graph *g, bitset remainingVertices,
 int lastElemOfPath, int firstElemOfPath, int *orderOfLongestInducedPath,
 unsigned long long int numberOfLengths[], int orderOfPath) {

    // Count frequency of induced paths. Length of a path is number of edges,
    // not vertices. Used for keeping track of graphs with forbidden induced
    // path sizes.
    numberOfLengths[orderOfPath - 1]++;

    // Check whether we have a longest induced path already.
    if(orderOfPath > *orderOfLongestInducedPath) {
        *orderOfLongestInducedPath = orderOfPath;
    }
    
    bitset neighboursOfLastNotInPath = 
     intersection(g->adjacencyList[lastElemOfPath], remainingVertices);

    forEach(neighbour, neighboursOfLastNotInPath) {

        int oldElemOfPath = lastElemOfPath;

        remainingVertices = 
         difference(remainingVertices, g->adjacencyList[oldElemOfPath]);
        bitset deletedVertices = 
         difference(g->adjacencyList[oldElemOfPath], remainingVertices);

        // Neighbour is the new last element.
        lastElemOfPath = neighbour; 

        searchLongestInducedSuperPath(g, remainingVertices, lastElemOfPath,
         firstElemOfPath, orderOfLongestInducedPath, numberOfLengths,
         orderOfPath + 1);

        remainingVertices = union(remainingVertices, deletedVertices);
        lastElemOfPath = oldElemOfPath;
    }
}

int getLongestInducedPathLength(struct graph *g, 
 unsigned long long int numberOfLengths[]) {

    int orderOfLongestInducedPath = 0;

    // For each vertex find a longest induced path starting with v.
    for(int v = 0; v < g->nv; v++) {

        bitset remainingVertices = complement(g->adjacencyList[v], g->nv);
        removeElement(remainingVertices, v);

        // The number of vertices gets stored inorderOfLongestInducedPath 
        forEach(w, g->adjacencyList[v]) {
            searchLongestInducedSuperPath(g, remainingVertices, w, v,
             &orderOfLongestInducedPath, numberOfLengths, 2);
        }
    }

    //  Length of a path is number of edges in it.
    int pathLength = orderOfLongestInducedPath - 1;
    if(pathLength < 0) pathLength = 0;

    return pathLength;
}

//******************************************************************************
//
//                          Parsing flags
//
//******************************************************************************

struct options {
    bool cycleFlag;
    bool pathFlag;
    bool differenceFlag;
    bool complementFlag;
    bool lengthFlag;
    int forbiddenLength;
    int output;
};

bool shouldOutput(struct graph *g, int length, 
 unsigned long long int numberOfLengths[],  int optionsNumber, 
 struct options *options) {

    bool forbiddenLengthTooLong;

    //  ComplementTrue is true if complementFlag is false, false if
    //  complementFlag is true.
    bool complementTrue = (!true != !options->complementFlag);
    bool complementFalse = (!false != !options->complementFlag);

    switch(optionsNumber) {

        //  Not present: -d, -f#.
        case 0:
            if(length == options->output) {
                return complementTrue;
            }
            return complementFalse;

        //  Present: -d, not present: -f#.
        case 1:
            if(g->nv - length == options->output) {
                return complementTrue;
            }
            return complementFalse;

        //  Present: -f#, Not present: -d
        case 2:

            // Length of path can be at most n-1 (edges!) of cycles at most n. 
            forbiddenLengthTooLong =
             (options->pathFlag && options->forbiddenLength >= g->nv) ||
             (options->cycleFlag && options->forbiddenLength > g->nv);

            if(forbiddenLengthTooLong || 
             numberOfLengths[options->forbiddenLength] == 0) {
                return complementTrue;
            }
            return complementFalse;

        case 3:
            fprintf(stderr, "Do not use -d with -f#.\n");
            fprintf(stderr, "%s\n", USAGE);
            exit(1);

        default:
            fprintf(stderr, "Invalid combination of options.\n");
            fprintf(stderr, "%s\n", USAGE);
            exit(1);
    }
}

void printTable(struct options *options,
 unsigned long long int frequencies[], char *tableString) {

    for (int i = 0; i < BITSETSIZE; ++i) {
        if(frequencies[i] != 0) {
            fprintf(stderr, "\n \t%16lld graphs: %s%s = %d",
             frequencies[i],
             options->differenceFlag ? "order - " : "",
             tableString, i);
        }
    }
    fprintf(stderr, "\n");
}

void printNumberGraphsOutput(struct options *options,
 long long unsigned int passedGraphs, char *tableString) {

    if(options->output != -1 || options->forbiddenLength != -1) {
        if(passedGraphs > 0) {
            fprintf(stderr, "%lld graphs sent to stdout.\n", passedGraphs);
        }
        else {
            if(options->output != -1) {
                fprintf(stderr, "\nNo graphs found %s %s%s %d \n",
                 options->complementFlag ? "without" : "with",
                 options->differenceFlag ? "order - " : "",
                 tableString, options->output);
            }
            if(options->forbiddenLength != -1) {
                fprintf(stderr, 
                 "\nNo graphs found %s induced %s of forbidden length %d\n",
                 options->complementFlag ? "with" : "without",
                 options->pathFlag ? "path" : "cycle", options->forbiddenLength);
            }
        }
    }
}

int main(int argc, char ** argv) {

    struct options options = {0};
    options.forbiddenLength = -1;
    options.output = -1;
    char* tableString = "circumference";

    int opt;
    while (1) {
        int option_index = 0;
        static struct option long_options[] = 
        {   
            {"induced-cycle", no_argument, NULL, 'c'},
            {"complement", no_argument, NULL, 'C'},
            {"difference", no_argument, NULL, 'd'},
            {"forbidden", required_argument, NULL, 'f'},
            {"help", no_argument, NULL, 'h'},
            {"length", no_argument, NULL, 'l'},
            {"output", required_argument, NULL, 'o'},
            {"induced-path", no_argument, NULL, 'p'},
        };

        opt = getopt_long(argc, argv, "cCdf:hlo:p", long_options, &option_index);
        if (opt == -1) break;
        switch(opt) {
            case 'c':
                options.cycleFlag = true;
                tableString = "longest induced cycle";
                break;
            case 'C':
                options.complementFlag = true;
                break;
            case 'd':
                options.differenceFlag = true;
                break;
            case 'f':
                options.forbiddenLength = (int) strtol(optarg, (char **)NULL, 10);
                break;
            case 'h':
                fprintf(stderr, "%s\n", USAGE);
                fprintf(stderr, "%s", HELPTEXT);
                return 0;
            case 'l':
                options.lengthFlag = true;
                tableString = "graph length";
                break;
            case 'o':
                options.output = (int) strtol(optarg, (char **)NULL, 10);
                break;
            case 'p':
                options.pathFlag = true;
                tableString = "longest induced path";
                break;
            case '?':
                fprintf(stderr,"Error: Unknown option: %c\n", optopt);
                fprintf(stderr, "%s\n", USAGE);
                fprintf(stderr,
                 "Use ./circumferenceChecker --help"
                 " for more detailed instructions.\n");
                return 1;
        }
    }

    if(options.forbiddenLength != -1 && options.output != -1) {
        fprintf(stderr,
         "Error: -f# should not be used with -o#.\n");
        fprintf(stderr, "%s\n", USAGE);
        return 1;       
    }

    if(options.cycleFlag && options.pathFlag) {
        fprintf(stderr, "Invalid combination of options.\n");
        fprintf(stderr, "%s\n", USAGE);
        return 1;
    }
    if((options.cycleFlag || options.pathFlag) && options.lengthFlag) {
        fprintf(stderr, "Invalid combination of options.\n");
        fprintf(stderr, "%s\n", USAGE);
        return 1;
    }
    if(!options.cycleFlag && !options.pathFlag && 
     options.forbiddenLength != -1) {
        fprintf(stderr, "Use -f only with -c or -p.\n");
        fprintf(stderr, "%s\n", USAGE);
        return 1;
    }

    int optionsNumber = (options.differenceFlag ? 1 : 0) |
                        (options.forbiddenLength != -1 ? 2 : 0);

    unsigned long long int counter = 0;
    unsigned long long int skippedGraphs = 0;
    unsigned long long int passedGraphs = 0;
    unsigned long long int frequencies[BITSETSIZE] =
     { [ 0 ... BITSETSIZE-1 ] = 0 };

    clock_t start = clock();

    //  Start looping over lines of stdin.
    char * graphString = NULL;
    size_t size;
    while(getline(&graphString, &size, stdin) != -1) {

        struct graph g;
        g.nv = getNumberOfVertices(graphString);
        if(g.nv == -1 || g.nv > BITSETSIZE - 1) {
            fprintf(stderr, "Skipping invalid graph!\n");
            skippedGraphs++;
            continue;
        }
        bitset adjacencyList[g.nv];
        if(loadGraph(graphString, g.nv, adjacencyList) == -1) {
            fprintf(stderr, "Skipping invalid graph!\n");
            skippedGraphs++;
            continue;
        }
        g.adjacencyList = adjacencyList;
        counter++;

        // Length is largest length of (induced) cycle(or path). numberOfLenghts
        // keeps track of each length encountered for the induced paths or
        // cycles. Used for not counting forbidden induced cycle or path
        // lengths.
        int length;
        unsigned long long int numberOfLengths[BITSETSIZE] = 
         { [ 0 ... BITSETSIZE-1 ] = 0 };

        if(options.cycleFlag) {
            length = getLongestInducedCycleLength(&g, numberOfLengths);
        }
        else if(options.pathFlag) {
            length = getLongestInducedPathLength(&g, numberOfLengths);
        }
        else if(options.lengthFlag) {
            length = getLength(&g);
        }
        else {
            length = getCircumference(&g, EMPTY);
        }

        if(shouldOutput(&g, length, numberOfLengths,  optionsNumber, &options)) {
            passedGraphs++;
            printf("%s", graphString);
        }
        if(options.differenceFlag) {
            frequencies[g.nv - length]++;
        }
        else {
            frequencies[length]++;
        }
    }

    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    free(graphString);

    // Print data
    printTable(&options, frequencies, tableString);

    // Mention how many graphs were output
    printNumberGraphsOutput(&options, passedGraphs, tableString);

    // Mention how many graphs checked
    fprintf(stderr,"\rChecked %lld graphs in %f seconds.\n",
     counter, time_spent);

    return 0;
}