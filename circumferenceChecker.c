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

#define USAGE "Usage: ./circumferenceChecker [-cf#|-pf#] [-Cdo#] [-h]"

#define HELPTEXT \
"Count and or filter graphs depending on their circumference, induced cycles\n\
or paths.\n\
\n\
Graphs are read from stdin in graph6 format. Graphs are sent to stdout in\n\
graph6 format. If the input graph had a graph6 header, so will the\n\
output graph (if it passes through the filter).\n\
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
            print help message\n\
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

void printGraph(bitset adjacencyList[], int numberOfVertices) {
    for(int i = 0; i < numberOfVertices; i++) {
        fprintf(stderr, "%d: ", i);
        forEach(nbr, adjacencyList[i]) {
            fprintf(stderr, "%d ", nbr);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

bool canBeCycleOfLength(bitset adjacencyList[], bitset remainingVertices, int
lastElemOfPath, int firstElemOfPath, int cycleLength, int pathLength) {

    if((pathLength == cycleLength) && contains(adjacencyList[firstElemOfPath], lastElemOfPath)) {
        return true;
    }
    if(isEmpty(intersection(adjacencyList[firstElemOfPath],remainingVertices))) { 
        return false;
    }

    bitset neighboursOfLastNotInPath = 
     intersection(adjacencyList[lastElemOfPath], remainingVertices);
    forEach(neighbour, neighboursOfLastNotInPath) {

        int oldElemOfPath = lastElemOfPath;

        removeElement(remainingVertices, neighbour);
        lastElemOfPath = neighbour; // Neighbour is the new last element.

        if (canBeCycleOfLength(adjacencyList, remainingVertices, lastElemOfPath,
         firstElemOfPath, cycleLength, pathLength + 1)) {
            return true;
        }

        add(remainingVertices, lastElemOfPath);
        lastElemOfPath = oldElemOfPath;
    }

    return false;
}

int getCircumference(bitset adjacencyList[], int numberOfVertices, bitset excludedVertices) {
    for(int i = numberOfVertices; i > 2; i--) {
        bitset forbiddenVertices = excludedVertices;
        for(int j = numberOfVertices - i; j >= 0; j--) {

            bitset includedVertices = complement(forbiddenVertices, numberOfVertices);

            if(isEmpty(includedVertices)) return false;

            // First included vertex.
            int startingVertex = next(includedVertices,-1);
            int lowestDegree = size(adjacencyList[startingVertex]);

            //  Find an included vertex of lowest degree.
            forEachAfterIndex(includedVertex, includedVertices, startingVertex){
                int degreeOfNeighbour = 
                 size(intersection(adjacencyList[includedVertex], includedVertices));
                if(lowestDegree > degreeOfNeighbour) {
                    lowestDegree = degreeOfNeighbour;
                    startingVertex = includedVertex;
                }
            }

            // Loop over included neighbours of startingVertex and for each such
            // neighbour loop over the included neighbours of startingVertex that are
            // of higher index.
            forEach(secondElemOfPath, intersection(adjacencyList[startingVertex], includedVertices)) {
                forEachAfterIndex(lastElemOfPath, intersection(adjacencyList[startingVertex], includedVertices), secondElemOfPath) {

                    // Create path, lastElemOfPath, startingVertex, secondElemOfPath.
                    // We have lastElemOfPath > secondElemOfPath, so that we
                    // eliminate the checking of paths which are mirrored.
                    bitset path = singleton(startingVertex);
                    add(path, lastElemOfPath);
                    add(path, secondElemOfPath);
                    bitset remainingVertices = difference(includedVertices, path);

                    if (canBeCycleOfLength(adjacencyList, remainingVertices,
                     lastElemOfPath, secondElemOfPath, i - size(excludedVertices), 3)) {
                        return i;
                    }
                }
            }
            add(forbiddenVertices, startingVertex);
        }
    }
    return 0;
}

void lengthOfLongestInducedSuperCycle(bitset adjacencyList[], bitset remainingVertices, int
 lastElemOfPath, int firstElemOfPath, int *longestCycleLength, unsigned long long int numberOfLengths[], int pathLength) {

    // Check whether we have a longest induced path already and whether this path
    // is a cycle.
    if(contains(adjacencyList[firstElemOfPath], lastElemOfPath)) {
        if(pathLength > *longestCycleLength) {
            *longestCycleLength = pathLength;
        }
        numberOfLengths[pathLength]++;
        return; //  Need to return because we already have a cycle.
    }

    // Check if cycle can still be closed with remaining vertices.
    if(isEmpty(intersection(adjacencyList[firstElemOfPath],remainingVertices))) { 
        return;
    }

    bitset neighboursOfLastNotInPath = 
     intersection(adjacencyList[lastElemOfPath], remainingVertices);
    forEach(neighbour, neighboursOfLastNotInPath) {

        int oldElemOfPath = lastElemOfPath;

        //  Delete neighbours of oldElemOfPath since no chord are allowed in the induced cycle.
        remainingVertices = difference(remainingVertices, adjacencyList[oldElemOfPath]);
        bitset deletedVertices = difference(adjacencyList[oldElemOfPath], remainingVertices);
        lastElemOfPath = neighbour;

        lengthOfLongestInducedSuperCycle(adjacencyList, remainingVertices, lastElemOfPath,
         firstElemOfPath, longestCycleLength, numberOfLengths, pathLength + 1);

        //  Once oldElemOfPath is again an end of the path, the neighbours we
        //  removed should be considered again.
        remainingVertices = union(remainingVertices, deletedVertices);
        lastElemOfPath = oldElemOfPath;
    }
}

int getLongestInducedCycleLength(bitset adjacencyList[], int numberOfVertices, unsigned long long int numberOfLengths[]) {
    int longestInducedCycleLength = 0;
    for(int i = 0; i < numberOfVertices; i++) {
        bitset includedVertices = complement(EMPTY, numberOfVertices);

        // Loop over included neighbours of startingVertex and for each such
        // neighbour loop over the included neighbours of startingVertex that are
        // of higher index.
        forEachAfterIndex(secondElemOfPath, intersection(adjacencyList[i], includedVertices), i) {
            forEachAfterIndex(lastElemOfPath, intersection(adjacencyList[i], includedVertices), secondElemOfPath) {

                bitset remainingVertices = difference(includedVertices, union(adjacencyList[i], singleton(i)));

                lengthOfLongestInducedSuperCycle(adjacencyList, remainingVertices,
                 lastElemOfPath, secondElemOfPath, &longestInducedCycleLength, numberOfLengths, 3);
            }
        }
    }
    return longestInducedCycleLength;
}

void lengthOfLongestInducedSuperPath(bitset adjacencyList[], bitset remainingVertices, int
 lastElemOfPath, int firstElemOfPath, int *longestCycleLength, unsigned long long int numberOfLengths[], int pathLength) {

    // Count frequency of induced paths. length of a path is number of edges, not vertices. 
    numberOfLengths[pathLength - 1]++;

    // Check whether we have a longest induced path already.
    if(pathLength > *longestCycleLength) {
        *longestCycleLength = pathLength;
    }
    
    bitset neighboursOfLastNotInPath = 
     intersection(adjacencyList[lastElemOfPath], remainingVertices);
    forEach(neighbour, neighboursOfLastNotInPath) {

        int oldElemOfPath = lastElemOfPath;

        remainingVertices = difference(remainingVertices, adjacencyList[oldElemOfPath]);
        bitset deletedVertices = difference(adjacencyList[oldElemOfPath], remainingVertices);
        lastElemOfPath = neighbour; // Neighbour is the new last element.

        lengthOfLongestInducedSuperPath(adjacencyList, remainingVertices, lastElemOfPath,
         firstElemOfPath, longestCycleLength, numberOfLengths, pathLength + 1);

        remainingVertices = union(remainingVertices, deletedVertices);
        lastElemOfPath = oldElemOfPath;
    }
}

int getLongestInducedPathLength(bitset adjacencyList[], int numberOfVertices, unsigned long long int numberOfLengths[]) {
    int longestInducedPathLength = 0;

    for(int i = 0; i < numberOfVertices; i++) {
        bitset remainingVertices = complement(adjacencyList[i], numberOfVertices);
        removeElement(remainingVertices, i);
        forEach(neighbour, adjacencyList[i]) {
            lengthOfLongestInducedSuperPath(adjacencyList, remainingVertices,
                neighbour, i, &longestInducedPathLength, numberOfLengths, 2);
        }
    }

    //  Length of a path is number of edges in it.
    return longestInducedPathLength - 1;
}

bool shouldOutput(int numberOfVertices, int length, unsigned long long int numberOfLengths[], int output, int forbiddenLength, int optionsNumber, bool complementFlag, bool pathFlag) {
    //  True if complementFlag is false, false if complementFlag is true.
    bool complementTrue = (!true != !complementFlag);
    bool complementFalse = (!false != !complementFlag);
    switch(optionsNumber) {
        case 0:
            if(length == output) {
                return complementTrue;
            }
            return complementFalse;
        case 1:
            if(numberOfVertices - length == output) {
                return complementTrue;
            }
            return complementFalse;
        case 2:
            if((pathFlag ? forbiddenLength >= numberOfVertices : forbiddenLength > numberOfVertices) || numberOfLengths[forbiddenLength] == 0) {
                return complementTrue;
            }
            return complementFalse;
        default:
            fprintf(stderr, "Invalid combination of options.\n");
            fprintf(stderr, "%s\n", USAGE);
            exit(1);
    }
}


int main(int argc, char ** argv) {
    bool cycleFlag = false;
    bool pathFlag = false;
    bool differenceFlag = false;
    bool complementFlag = false;
    int forbiddenLength = -1;
    char* tableString = "circumference";
    int output = -1;
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
            {"output", required_argument, NULL, 'o'},
            {"induced-path", no_argument, NULL, 'p'},
        };

        opt = getopt_long(argc, argv, "cCdf:ho:p", long_options, &option_index);
        if (opt == -1) break;
        switch(opt) {
            case 'c':
                cycleFlag = true;
                tableString = "longest induced cycle";
                break;
            case 'C':
                complementFlag = true;
                break;
            case 'd':
                differenceFlag = true;
                tableString = "order - circumference";
                break;
            case 'f':
                forbiddenLength = (int) strtol(optarg, (char **)NULL, 10);
                break;
            case 'h':
                fprintf(stderr, "%s\n", USAGE);
                fprintf(stderr, "%s", HELPTEXT);
                return 0;
            case 'o':
                output = (int) strtol(optarg, (char **)NULL, 10);
                break;
            case 'p':
                pathFlag = true;
                tableString = "longest induced path";
                break;
            case '?':
                fprintf(stderr,"Error: Unknown option: %c\n", optopt);
                fprintf(stderr, "%s\n", USAGE);
                fprintf(stderr,
                 "Use ./circumferenceChecker --help for more detailed instructions.\n");
                return 1;
        }
    }

    if(forbiddenLength != -1 && output != -1) {
        fprintf(stderr,
         "Warning: -f need not be used with -o. This will probably not give the result you expect.\n");
        fprintf(stderr,
         "Use ./circumferenceChecker --help for more detailed instructions.\n");
    }

    if(cycleFlag && pathFlag) {
        fprintf(stderr, "Invalid combination of options.\n");
        fprintf(stderr, "%s\n", USAGE);
        return 1;
    }
    if(!cycleFlag && !pathFlag && forbiddenLength != -1) {
        fprintf(stderr, "Use -f only with -c or -p.\n");
        fprintf(stderr, "%s\n", USAGE);
        return 1;
    }

    int optionsNumber = (differenceFlag ? 1 : 0) |
                        (forbiddenLength != -1 ? 2 : 0);

    unsigned long long int counter = 0;
    unsigned long long int skippedGraphs = 0;
    unsigned long long int passedGraphs = 0;
    unsigned long long int frequencies[BITSETSIZE] = { [ 0 ... BITSETSIZE-1 ] = 0 };

    clock_t start = clock();

    //  Start looping over lines of stdin.
    char * graphString = NULL;
    size_t size;
    while(getline(&graphString, &size, stdin) != -1) {
        int nVertices = getNumberOfVertices(graphString);
        if(nVertices == -1 || nVertices > BITSETSIZE) {
            fprintf(stderr, "Skipping invalid graph!\n");
            skippedGraphs++;
            continue;
        }
        bitset adjacencyList[nVertices];
        if(loadGraph(graphString, nVertices, adjacencyList) == -1) {
            fprintf(stderr, "Skipping invalid graph!\n");
            skippedGraphs++;
            continue;
        }
        counter++;


        int length;
        unsigned long long int numberOfLengths[nVertices + 1];
        memset(numberOfLengths,0, (nVertices+1)*sizeof(unsigned long long int));
        if(cycleFlag) {
            length = getLongestInducedCycleLength(adjacencyList, nVertices, numberOfLengths);
        }
        else if(pathFlag) {
            length = getLongestInducedPathLength(adjacencyList, nVertices, numberOfLengths);
        }
        else {
            length = getCircumference(adjacencyList, nVertices, EMPTY);
        }

        if(shouldOutput(nVertices, length, numberOfLengths, output, forbiddenLength, optionsNumber, complementFlag, pathFlag)) {
            passedGraphs++;
            printf("%s", graphString);
        }
        if(differenceFlag) {
            frequencies[nVertices - length]++;
        }
        else {
            frequencies[length]++;
        }
    }
    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;

    for (int i = 0; i < BITSETSIZE; ++i) {
        if(frequencies[i] != 0) {
            fprintf(stderr, "\n \t%16lld graphs: %s = %d", frequencies[i],tableString, i);
        }
    }
    fprintf(stderr, "\n");
    if(passedGraphs > 0 || output != -1 || forbiddenLength != -1) {
        if(frequencies[output] > 0 || complementFlag) {
            fprintf(stderr, "%lld graphs sent to stdout.\n", passedGraphs);
        }
        else
            fprintf(stderr, "\nNo graphs found with %s %d \n", tableString, output);
    }

    fprintf(stderr,"\rChecked %lld graphs in %f seconds.\n", counter, time_spent);


    return 0;
}