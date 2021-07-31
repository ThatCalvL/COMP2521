///////////////////////////////////////////////////////////////////////

// COMP2521 Assignment 2 Part 2
// Written by Calvin Li z5242094
// 8 Nov 2020

///////////////////////////////Headers////////////////////////////////////////
#include "CentralityMeasures.h"
#include "Dijkstra.h"
#include "Graph.h"
#include "PQ.h"
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

#define INFINITY INT_MAX

////////////////////////////// Helpers///////////////////////////////////////////
static void unvisitAll(int visited[], int size);
//count # edges from src vertex.
static double countEdge(Graph g, Vertex src, ShortestPaths path);
//calculate sum of all paths.
static double sumPath(Graph g, ShortestPaths path);
//find Tthe amount of possible paths between src and dest.
static int cntPaths (PredNode **pred, Vertex src, Vertex dest, int visited[]);
// return the total number of paths i.e. the denominator.
static double numPath (ShortestPaths pathOne, Graph g, int visited[], int i, int j);
// count the #common paths of both part that passes through v;
// Inspired by Adrian in the tutorial assignment helping session.
static double numPassV (ShortestPaths pathOne, Graph g, int v, int src, int dest, int visited[]);


////////////////////////////// Task 1//////////////////////////////////////////
NodeValues closenessCentrality(Graph g){
    NodeValues c;
    c.numNodes = GraphNumVertices(g) - 1;
    c.values = malloc(GraphNumVertices(g) * sizeof(double));
    // run dijkstra on each node in graph, each time, add up the distance and
    // save in c.values array in respective index
    for (int src = 0; src < GraphNumVertices(g); src++) {
        ShortestPaths newpath = dijkstra(g, src);
        //count # vertices that i can reach
        double n = countEdge(g, src, newpath);

        //calculate sum of all paths
        double sum = 0.0;
        sum = sumPath(g, newpath);

        //get final result
        double N = (double)(c.numNodes);
        double numerator;
        double denominator;
        // check whether the demonimator is 0
        // if so, set the value to 0.0.
        if (N == 0.0 || sum == 0.0) {
            c.values[src] = 0.0;
        }
        else{
            numerator = (n - 1.0)*(n - 1.0);
            denominator = (N)*sum;
            c.values[src] = numerator/denominator;
        }
    }
     return c;
}

////////////////////////////// Task 2 //////////////////////////////////////////
NodeValues betweennessCentrality(Graph g){
    NodeValues b;
    // set variable and malloc double array
    b.numNodes = GraphNumVertices(g) - 1;
    b.values = malloc (GraphNumVertices(g) * sizeof(double));
    int a;
    // set all values to zero
    for (a = 0; a < GraphNumVertices(g); a++) {
        b.values[a] = 0;
    }
    // loop through all the source nodes
    int src; // source
    for (src = 0; src < GraphNumVertices(g); src++) {
        int dest;
        // loop through all the destination nodes
        for (dest = 0; dest < GraphNumVertices(g); dest++) {
            // skip when destination is source
            if (dest != src) {
                // run dijkstra on the source
                ShortestPaths path = dijkstra(g, src);
                // create visited array to be passed in to numPaths
                int visited[GraphNumVertices(g)];
                double total = numPath (path, g, visited, src, dest);
                // skip if total == 0;
                if (total != 0) {
                    int v;
                    // loop through all nodes in the graph to calculate
                    // betweenness
                    for (v = 0; v < GraphNumVertices(g); v++) {
                        // check if v is equal to either i or j
                        double result = numPassV (path, g, v, src, dest, visited);
                        // divide with total and add to the value in the
                        // respective index
                        double betweenness = result / total;
                        b.values[v] += betweenness;
                    }
                }
            }
        }
    }
    return b;
}


////////////////////////////// Task 3 //////////////////////////////////////////
NodeValues betweennessCentralityNormalised(Graph g){

    NodeValues bn  = betweennessCentrality(g);
    int i;
    double n = GraphNumVertices(g);
    for (i = 0; i < GraphNumVertices(g); i++) {
        double numer = bn.values[i];
        double denom = (n - 1.0) * (n - 2.0);
        double result = numer/denom;
        bn.values[i] = result;
    }
    return bn;
}


///////////////////////////////////////////////////////////////////////
void showNodeValues(NodeValues values){
    int i;
    for (i = 0; i <= values.numNodes; i++) {
        printf("%d: %.6lf\n", i, values.values[i]);
    }
}

void freeNodeValues(NodeValues values){
    free(values.values);
}
///////////////////////////////////////////////////////////////////////

////////////////////////////Helper Implementation//////////////////////

// marks all nodes in visited array to unvisited
static void unvisitAll(int visited[], int size) {
    int i;
    for (i = 0; i < size; i++) {
        visited[i] = -1;
    }
}

static double countEdge(Graph g, Vertex src, ShortestPaths path) {
    double count = 1.0;
    int j;
    for (j = 0; j < GraphNumVertices(g); j++) {
        if (path.pred[j] != NULL) {
            count++;
        }
    }
    return count;
}

static double sumPath(Graph g, ShortestPaths path) {
    double dist = 0.0;
    int k;
    double sum = 0.0;
    for (k = 0; k < GraphNumVertices(g); k++) {
        dist = (double)(path.dist[k]);
        sum = sum + dist;
    }
    return sum;
}

// Adapted from depth-first search in COMP2521 lec slides
// The amount of possible paths found between src and dest is returned
static int cntPaths (PredNode **pred, Vertex src, Vertex dest, int visited[]) {
    if (src == dest) {
        return 1;
    }
    int count = 0;
    PredNode *curr = pred[src];
    while (curr != NULL) {
        if (visited[curr->v] == -1) {
            visited[curr->v] = 1;
            if (curr->v == dest) {
                visited[curr->v] = -1;
            }
            count += cntPaths(pred, curr->v, dest, visited);
        }
        curr = curr->next;
    }
    visited[src] = -1;
    return count;
}

// return the total number of paths i.e. the denominator.
static double numPath (ShortestPaths pathOne, Graph g, int visited[], int src, int dest) {
    // set all nodes to be unvisited
    unvisitAll(visited, GraphNumVertices(g));
    // run numPaths to first get total number of paths for the
    // current src, dest pair
    double total = cntPaths(pathOne.pred, dest, src, visited);
    return total;
}

// count the #common paths of both part that passes through v;
//inpired by Adrian Martinez the tutor of this term.
// He suggests the way of computing the #paths through v from src to dest.
static double numPassV (ShortestPaths path, Graph g, int v, int src, int dest, int visited[]) {
    double result = 0.0;
    // check if v is equal to either i or j
    if (v != src && v != dest) {
        unvisitAll(visited, GraphNumVertices(g));
        // counts number of paths from src to i
        int partOne = cntPaths(path.pred, v, src, visited);
        unvisitAll(visited, GraphNumVertices(g));
        // counts number of paths from i to dest
        int partTwo = cntPaths(path.pred, dest, v, visited);
        // multiply to get total number of paths that i
        // appears in
        result = partOne * partTwo;
    }
    return result;
}
