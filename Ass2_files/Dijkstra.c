///////////////////////////////////////////////////////////////////////

// COMP2521 Assignment 2 Part 1
// Written by Calvin Li z5242094
// 4 Nov 2020

/////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include "Dijkstra.h"
#include "PQ.h"
#include "Graph.h"

// Task 1
#define INFINITY 9999999
/////////////////////////////////////////////////////////////////////////
/*---------------------------HELPER----------------------------*/
static void updatePred(PredNode **list, int v, int w) ;
static ShortestPaths newShortPath(Graph g, Vertex src);
static PredNode* newPredNode(Vertex v);

/////////////////////////////////////////////////////////////////////////

ShortestPaths dijkstra(Graph g, Vertex src) {
	//create a new shortestpath ADT && Initialize it with vertex src
	ShortestPaths new;
	new = newShortPath(g, src);

	//init dist and priority queue
	PQ q = PQNew();
	int i;
	for (i = 0; i < new.numNodes; i++) {
		new.dist[i] = INFINITY;
		new.pred[i] = NULL;
		PQInsert(q, i, new.dist[i]);
	}
	new.dist[src] = 0;
	new.pred[src] = NULL;

	PQInsert(q, src, new.dist[src]);

	while (!PQIsEmpty(q)) {
		Vertex s = PQDequeue(q);
		AdjList out = GraphOutIncident(g, s); // find the weight for each dequeued vertex
		AdjList curr = out;
		while (curr != NULL) {
            Vertex v = curr->v;
			int len = curr->weight;
			if ((new.dist[s] + len) < (new.dist[v])) {     
                new.dist[v] = new.dist[s] + len;
                PredNode *k = new.pred[v];
                while (k != NULL) {
                    PredNode *eraseNow = k;
                    k = k->next;
                    free(eraseNow);
                }
                PredNode *newNode = newPredNode(s);
                new.pred[v] = newNode;
                PQInsert(q, v, new.dist[v]);
            } else if ((new.dist[s] + len) == (new.dist[v])) {
				//update dist
				new.dist[curr->v] = (new.dist[s]) + len;
				//update queue
				PQUpdate(q, curr->v, new.dist[curr->v]);
				//update pred
				// if the vertex has beem visited but find a path whose distance is same as previous
                updatePred(new.pred, s, curr->v);
            }
			curr = curr->next;
		}
	}
	// checck every distance and set them to 0 if infi.
	int c;
	for (c = 0; c < GraphNumVertices(g); c++) {
		if (new.pred[c] == NULL) {
			new.dist[c] = 0;
		}
	}

	PQFree(q);
	return new;
}

/////////////////////////////////////////////////////////////////////////

void  showShortestPaths(ShortestPaths s) {
	int i = 0;
	printf("Node %d\n",s.src);
	printf("  Distance\n");
	for (i = 0; i < s.numNodes; i++) {
		if(i == s.src)
	    	printf("    %d : X\n",i);
		else
			printf("    %d : %d\n",i,s.dist[i]);
	}
	printf("  Preds\n");
	for (i = 0; i < s.numNodes; i++) {
		printf("    %d : ",i);
		PredNode* curr = s.pred[i];
		while(curr!=NULL) {
			printf("[%d]->",curr->v);
			curr = curr->next;
		}
		printf("NULL\n");
	}
}

void  freeShortestPaths(ShortestPaths s) {
	int i;
	//free distance
	free(s.dist);
	//free pred
	for (i = 0; i < s.numNodes; i++) {
		if (s.pred[i] == NULL) continue;
		PredNode *curr = s.pred[i];
		while (curr != NULL) {
			PredNode *temp = curr->next;
			free(curr);
			curr = temp;
		}
	}
	free(s.pred);
}

/////////////////////////////////////////////////////////////////////////
/*----------------HELPER IMPLEMENTATION--------------*/
static void updatePred(PredNode **list, int v, int w) {
	//create new node
	PredNode *new = malloc(sizeof(PredNode));
	assert(new != NULL);
	new->v = v;
	new->next = list[w];
	list[w] = new;
}


static ShortestPaths newShortPath(Graph g, Vertex src){
	ShortestPaths new;
	new.dist = malloc(sizeof(int) * (GraphNumVertices(g)));
	new.pred = malloc(sizeof(struct PredNode) * (GraphNumVertices(g)));
	assert((new.dist) != NULL);
	assert((new.pred) != NULL);
	new.src = src;
	new.numNodes = GraphNumVertices(g);
	return new;
}

// creates newPredNode wiith given vertex
static PredNode* newPredNode(Vertex v) {
    PredNode *p = malloc(sizeof(struct PredNode));
    p->v = v;
    p->next = NULL;
    return p;
}
