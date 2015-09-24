/*
Copyright (c) 2013-2015 Daniel Marbach

We release this software open source under an MIT license (see below). If this
software was useful for your scientific work, please cite our paper available at:
http://regulatorycircuits.org

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 */
package edu.mit.magnum.netprop;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import edu.uci.ics.jung.algorithms.importance.BetweennessCentrality;

import edu.mit.magnum.*;
import edu.mit.magnum.net.*;

/**
 * Compute topological properties of nodes
 * - Degree (indegree, outdegree)
 * - Betweenness
 * - Clustering coefficient
 */
public class BasicProperties extends NetworkProperties {
	
	/** Node degree */
	private Integer degree_[] = null;
	/** Node outdegree (for directed networks) */
	private Integer outdegree_[] = null;
	/** Node indegree (for directed networks) */
	private Integer indegree_[] = null;
	/** Betweenness centrality */
	private Double betweenness_[] = null;
	/** Clustering coefficient */
	private Double clusteringCoeff_[] = null;
	/** Mean network clustering coefficient */
	private double meanClusteringCoeff_ = -1;
	
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public BasicProperties(Magnum mag, Network network) {
		super(mag, network);
	}
	
	
	// ----------------------------------------------------------------------------

	/** Compute all metrics specified in Settings */
	public ArrayList<Double> run() {
		
		// Compute all metrics
		if (mag.set.computeDegree_)
			computeDegree();
		if (mag.set.computeBetweenness_)
			computeBetweenness();
		if (mag.set.computeClusteringCoefficient_)
			computeClusteringCoefficient();
		
		ArrayList<Double> networkMeans = new ArrayList<Double>();
		networkMeans.add(meanClusteringCoeff_);
		return networkMeans;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Add the computed node properties of this analyzer to the given map */
	public void addNodeProperties(LinkedHashMap<String,Number[]> map) {
		
		if (degree_ != null) 
			map.put("degree", degree_);
		if (outdegree_ != null)
			map.put("outdegree", outdegree_);
		if (indegree_ != null)
			map.put("indegree", indegree_);
		if (betweenness_ != null)
			map.put("betweennessCentrality", betweenness_);
		if (clusteringCoeff_ != null)
			map.put("clusteringCoeff", clusteringCoeff_);
	}

	
	// ----------------------------------------------------------------------------

	/** Compute the degree for every node */
	public void computeDegree() {
		
		mag.log.println("Computing node degrees...");
		
		// Degree
		degree_ = new Integer[numNodes_];
		nodeProperties_.put("degree", degree_);
		
		for (int i=0; i<numNodes_; i++) {
			Node node = network_.getNode(i);
			degree_[i] = graph_.degree(node);
		}
		
		// Indegree and outdegree
		if (isDirected_) {
			indegree_ = new Integer[numNodes_];
			outdegree_ = new Integer[numNodes_];
			nodeProperties_.put("indegree", indegree_);
			nodeProperties_.put("outdegree", outdegree_);

			for (int i=0; i<numNodes_; i++) {
				Node node = network_.getNode(i);
				indegree_[i] = graph_.inDegree(node);
				outdegree_[i] = graph_.outDegree(node);
			}
		}
	}

	
	// ----------------------------------------------------------------------------

	/** Compute the betweenness for every node */
	public void computeBetweenness() {
		
		mag.log.println("Computing betweenness centrality...");

		// Compute betweenness
		BetweennessCentrality<Node, Edge> betweennessScorer = new BetweennessCentrality<Node, Edge>(graph_);
		betweennessScorer.setRemoveRankScoresOnFinalize(false);
		betweennessScorer.evaluate();
				
		// Get result
		betweenness_ = new Double[numNodes_];
		nodeProperties_.put("betweenness", betweenness_);

		for (int i=0; i<numNodes_; i++) {
			Node node = network_.getNode(i);
			betweenness_[i] = betweennessScorer.getVertexRankScore(node);
		}
	}


	// ----------------------------------------------------------------------------

	/** Compute the clustering coefficient for every node */
	public void computeClusteringCoefficient() {
		
		mag.log.println("Computing clustering coefficient...");

		// Compute clustering coefficient
		//Map<Node, Double> clust = Metrics.clusteringCoefficients(graph_);
		
		clusteringCoeff_ = new Double[numNodes_];
		meanClusteringCoeff_ = 0;
		nodeProperties_.put("clusteringCoeff", clusteringCoeff_);
		
		for (int i=0; i<numNodes_; i++) {
			Node node = network_.getNode(i);
			ArrayList<Node> neighbors = new ArrayList<Node>(network_.getNeighborsNoSelf(node));
			int numNeighbors = neighbors.size();
			
			// If the node has 0 or 1 neighbors, the clustering coeff is 0 
			if (numNeighbors < 2) {
				clusteringCoeff_[i] = 0.0;
				continue;
			}
			
			// The number of edges among neighbors
			int numEdges = 0;
			
			// For each directed pair of neighbors, count if there is an edge
			// This gives the correct result both for directed and undirected networks
            for (int k=0; k<numNeighbors; k++)
            	for (int l=0; l<numNeighbors; l++)
            		if (l != k && graph_.isPredecessor(neighbors.get(k), neighbors.get(l)))
            			numEdges++;

            // Number of observed divided by number of possible edges between neighbors
            double numPossibleEdges = numNeighbors * (numNeighbors-1);
            clusteringCoeff_[i] = numEdges / numPossibleEdges;
            meanClusteringCoeff_ += clusteringCoeff_[i];
		}
		meanClusteringCoeff_ /= numNodes_;
	}

	
	// ============================================================================
	// PRIVATE METHODS


	
	// ============================================================================
	// SETTERS AND GETTERS

	public Integer[] getDegree() { return degree_; }
	public Integer[] getOutdegree() { return outdegree_; }
	public Integer[] getIndegree() { return indegree_; }
	public Double[] getBetweenness() { return betweenness_; }
	public Double[] getClusteringCoeff() { return clusteringCoeff_; }
	
}
