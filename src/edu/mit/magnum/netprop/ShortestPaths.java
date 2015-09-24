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

import java.util.Collection;
import java.util.HashSet;
import java.util.Map;

import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraDistance;
import edu.mit.magnum.Magnum;
import edu.mit.magnum.ProgressMonitor;
import edu.mit.magnum.net.*;


/**
 * Compute properties based on shortest paths between nodes
 */
public class ShortestPaths extends PairwiseProperties {
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public ShortestPaths(Magnum mag, Network network, boolean computeCentrality) {
		
		super(mag, network, "shortestPaths", "closenessCentrality", computeCentrality);
	}
	
	
	// ----------------------------------------------------------------------------

	/** Compute the degree for every node */
	public void computeK() {
		
		mag.log.println("Computing shortest paths...");
		
		// Dijkstra distance
		K_ = new DenseDoubleMatrix2D(numNodes_, numRefNodes_);
		
		// Initialize at -1
		for (int i=0; i<numNodes_; i++)
			for (int j=0; j<numRefNodes_; j++)
				K_.set(i, j, -2);
		
		//boolean averaging = true;
		//boolean ignoreMissing = true;
		//boolean ignoreSelfDistances = true;
		//DistanceCentralityScorer<Node, Edge> scorer = new DistanceCentralityScorer<Node, Edge>(graph_, averaging, ignoreMissing, ignoreSelfDistances);
		
		// The code below is inspired from DistanceCentralityScorer
		
		// Class to compute Dijkstra distances
		boolean cached = true;
		DijkstraDistance<Node,Edge> dijkstra = new DijkstraDistance<Node,Edge>(graph_, cached);
		
		// DijkstraShortestPath is a lightweight extension of DijkstraDistance, allows to reconstruct the shortest paths
		//DijkstraShortestPath<Node,Edge> distance = new DijkstraShortestPath<Node,Edge>(graph_, cached);

		// Print status after doing 'freq' regulators
		ProgressMonitor progress = new ProgressMonitor(mag, numNodes_);

		for (int i=0; i<numNodes_; i++) {
			// Print progress
			progress.iteration(i);

			Node source = network_.getNode(i);
			// The set of target nodes for which the distance needs to be computed
			Collection<Node> targets = null;
			
			// For undirected networks, do not compute distance between prior nodes twice
			if (isDirected_ || !network_.isRefNode(source)) {
				targets = network_.getRefNodes();
				
			} else {
				// The set of priors for which distance has not yet been computed
				targets = new HashSet<Node>();

				for (int j=0; j<numRefNodes_; j++) {
					Node priorNode = network_.getRefNode(j);
					
					int m = network_.getNodeIndex(priorNode);
					int n = network_.getRefNodeIndex(source);
					
					if (K_.get(m, n) != -2)
						K_.set(i, j, K_.get(m, n));
					else
						targets.add(priorNode);
				}
			}
			
			// Compute Dijkstra distance to all prior nodes
		    Map<Node,Number> v_distances = dijkstra.getDistanceMap(source, targets);
		    		    
		    // Store in distance matrix
			for (Node priorNode : targets) {
				Double d = (Double) v_distances.get(priorNode);
				int k = network_.getRefNodeIndex(priorNode);

				// d is null for non-reachable nodes and zero on the diagonal ("self distance")
				if (d == null)
					K_.set(i, k, -1);
				else
					K_.set(i, k, d);
			}		    
		}
		progress.done();
	}

	
	// ----------------------------------------------------------------------------

	/** Closeness centrality (the mean shortest path length of a node to all prior nodes) */
	public void computeCentrality() {
		
		// Initialize
		centrality_ = new Double[numNodes_];
		for (int i=0; i<numNodes_; i++)
			centrality_[i] = -1.0;
		
		for (int i=0; i<numNodes_; i++) {
			double sum = 0;
			int numReachable = 0;
			
			for (int j=0; j<numRefNodes_; j++) {
				// Distance is -1 if there is no path
				assert K_.get(i, j) == -1 || K_.get(i, j) >= 0;
				// Self distance was ignored, i.e., should be -1
				assert network_.getNode(i).equals(network_.getRefNode(j)) ? K_.get(i, j) == 0 : true;
				// Ignore missing (-1)
				if (K_.get(i, j) >= 0) {
					sum += K_.get(i, j);
					numReachable++;
				}
			}
					
			// The node itself is by definition reachable with distance 0
			// By definition, we give isolated nodes closeness centrality 0
			if (sum == 0) {
				if (numReachable == 0)
					centrality_[i] = 0.0;
				else
					centrality_[i] = 1e12; // positive infinity (let's not go all the way up for compatibility when reading the files)
				
			} else if (network_.getUseRefNodes()){
				// Includes the 0 self-distance if this node is a reference node
				centrality_[i] = numReachable / sum;
			
			} else {
				// Ignore the 0 self-distance => sum stays the same, numReachable decreases by one
				centrality_[i] = (numReachable-1.0) / sum;
			}
		}
	}

	
	// ============================================================================
	// PRIVATE METHODS
	
	// ============================================================================
	// SETTERS AND GETTERS

}
