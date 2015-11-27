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
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import ch.unil.gpsutils.ProgressMonitor;
import edu.mit.magnum.Magnum;
import edu.mit.magnum.net.*;


/**
 * Compute Tanimoto coefficient for:
 * (1) Regulator sets of target genes
 * (2) Target gene sets of regulators
 *     
 * The network has to be directed, self-loops are allowed.
 */
public class TanimotoCoefficient extends PairwiseProperties {
		
	/** Compute similarity between targets, otherwise compute similarity between regulators */
	private boolean computeTargetSimilarity_ = true;
	
	private ArrayList<Node> nodes_ = null;
	/** The nodes for which pairwise similarity is computed (the targets or the regulators) and their index */
	private HashMap<Node, Integer> nodeIndexes_ = null;
	
	/** The sum of squares of the edge weights for each node */
	private double[] sumOfSquares_ = null;
	/** The TFs / targets of the nodes */
	private ArrayList<HashSet<Node>> opposites_ = null;
	
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor (network must be undirected) */
	public TanimotoCoefficient(Magnum mag, Network network, boolean computeTargetSimilarity, boolean computeCentrality) {
		
		super(mag, network, "tmp", "tmp", computeCentrality);
		computeTargetSimilarity_ = computeTargetSimilarity;
		if (computeTargetSimilarity_) {
			name_ = "targetTanimoto";
			nameCentrality_ = "targetTanimoto";
		} else {
			name_ = "tfTanimoto";
			nameCentrality_ = "tfTanimoto";
		}
		
		if (!network_.getIsDirected())
			throw new RuntimeException("Tanimoto similarity only implemented for directed networks");
	}
	
	
	// ----------------------------------------------------------------------------

	/** Compute the pairwise tanimoto coefficient for all target genes or regulators */
	public void computeK() {

		mag.log.println("Computing pairwise Tanimoto coefficient for " + (computeTargetSimilarity_ ? "TARGETS" : "TFs") + "...");

		// Initialize nodes
		initialize();
		
		K_ = new DenseDoubleMatrix2D(numNodes_, numNodes_);
		ProgressMonitor progress = new ProgressMonitor(mag.log, numNodes_);
		
		for (int i=0; i<numNodes_; i++) {
			progress.iteration(i);
			K_.set(i, i, 1);
			
			for (int j=i+1; j<numNodes_; j++) {
				double x = computeTanimoto(i, j);
				K_.set(i, j, x);
				K_.set(j, i, x);
			}
		}
		progress.done();
	}

	
	// ============================================================================
	// PRIVATE METHODS

	/** Initialize */
	private void initialize() {

		// Get the nodes for which similarity is to be computed as a list
		if (computeTargetSimilarity_)
			nodes_ = new ArrayList<Node>(network_.getTargetNodes());
		else
			nodes_ = new ArrayList<Node>(network_.getRegulatorNodes());
		numNodes_ = nodes_.size();	
		numRefNodes_ = numNodes_;
		// Sort nodes alphabetically
		Collections.sort(nodes_, Node.getNodeComparator());
		
		// Create hashmap with node index
		nodeIndexes_ = new HashMap<Node, Integer>(numNodes_);
		for (int i=0; i<numNodes_; i++) 
			nodeIndexes_.put(nodes_.get(i), i);
		
		// Initialize sum of squares and opposites
		opposites_ = new ArrayList<HashSet<Node>>(numNodes_);
		for (int i=0; i<numNodes_; i++)
			opposites_.add(new HashSet<Node>());
		sumOfSquares_ = new double[numNodes_];
		
		for (int i=0; i<numNodes_; i++) {
			Node node_i = nodes_.get(i);
			HashSet<Node> opposites_i = opposites_.get(i);
			sumOfSquares_[i] = 0;
			
			// Get the in / out edges
			if (computeTargetSimilarity_) {
				for (Edge edge_k : graph_.getInEdges(node_i)) {
					opposites_i.add(graph_.getSource(edge_k));
					sumOfSquares_[i] += edge_k.w_ * edge_k.w_; 
				}		
			} else {
				for (Edge edge_k : graph_.getOutEdges(node_i)) {
					opposites_i.add(graph_.getDest(edge_k));
					sumOfSquares_[i] += edge_k.w_ * edge_k.w_; 
				}
			}
			assert opposites_i.size() > 0;
			opposites_.add(opposites_i);
		}
	}
	

	// ----------------------------------------------------------------------------

	/** Compute similarity for a keyed edge pair (given the two end nodes) */
	private double computeTanimoto(int i, int j) {
		
		Node n_i = nodes_.get(i);
		Node n_j = nodes_.get(j);
		
		HashSet<Node> opposites_i = opposites_.get(i);
		HashSet<Node> opposites_j = opposites_.get(j);
		
		if (n_i == n_j)
			return 1.0;
				
		// a_i * a_j
		double product = 0;

		// Compute product
		if (computeTargetSimilarity_) {
			for (Node n_ik : opposites_i) {
				if (opposites_j.contains(n_ik)) {
					Edge e_i = graph_.findEdge(n_ik, n_i);
					Edge e_j = graph_.findEdge(n_ik, n_j);
					product += e_i.w_ * e_j.w_;
				}
			}
		} else {
			for (Node n_ik : opposites_i) {
				if (opposites_j.contains(n_ik)) {
					Edge e_i = graph_.findEdge(n_i, n_ik);
					Edge e_j = graph_.findEdge(n_j, n_ik);
					product += e_i.w_ * e_j.w_;
				}
			}			
		}
			
		// Tanimoto coefficient
		return product / (sumOfSquares_[i] + sumOfSquares_[j] - product);
	}



	// ============================================================================
	// SETTERS AND GETTERS

	public ArrayList<Node> getNodes() { return nodes_; }
	
}
