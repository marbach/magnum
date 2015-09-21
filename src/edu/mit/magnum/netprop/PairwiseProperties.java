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

import cern.colt.matrix.DoubleMatrix2D;
import edu.mit.magnum.FileExport;
import edu.mit.magnum.Magnum;
import edu.mit.magnum.MagnumUtils;
import edu.mit.magnum.net.*;


/**
 * Analyzer for pairwise node properties (similarity / distance)
 */
public abstract class PairwiseProperties extends NetworkProperties {
		
	/** Pairwise similarity/distance matrix (NxN) */
	protected DoubleMatrix2D K_ = null;
	/** Node centrality (avg. similarity/distance) */
	protected Double[] centrality_ = null;
	
	/** Set true to compute node centralities */
	protected boolean computeCentrality_ = false;
	/** Flag indicating that K_ and centrality_ have been exported to a file */
	protected boolean saved_ = false;
	
	/** The name of the distance/similarity measure (e.g., 'kstepKernel', 'shortestPath', etc.) */
	protected String name_ = null;
	/** The name of the centrality measure (e.g., 'kstepKernelCentrality', 'closenessCentrality', etc.) */
	protected String nameCentrality_ = null;
	
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor (network must be undirected) */
	public PairwiseProperties(Network network, String name, String nameCentrality, boolean computeCentrality) {
		
		super(network);
		name_ = name;
		nameCentrality_ = nameCentrality;
		computeCentrality_ = computeCentrality;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Compute all metrics specified in Settings */
	public ArrayList<Double> run() {
		
		long t0 = System.currentTimeMillis();
		computeK();
		long t1 = System.currentTimeMillis();
		Magnum.log.println("Run time: " + MagnumUtils.chronometer(t1-t0));

		if (computeCentrality_)
			computeCentrality();
		return null;
	}

	
	// ----------------------------------------------------------------------------

	/** Compute all metrics specified in Settings */
	abstract public void computeK();

	
	// ----------------------------------------------------------------------------

	public void addNodeProperties(LinkedHashMap<String,Number[]> map) {
		
		if (centrality_ != null) 
			map.put(nameCentrality_, centrality_);
	}

	
	// ----------------------------------------------------------------------------

	/** Export distance / similiarity matrix K_ (if saved_ is not true) */
	public void saveK(String basicFilename) {
		
		if (saved_)
			return;
		
		// Ref nodes are only implemented for ShortestPaths
		if (network_.getUseRefNodes())
			throw new RuntimeException("Ref nodes implementation incomplete");
		
		// Add suffix to filename
		String filename = Magnum.set.outputDirectory_ + "/" + basicFilename + "_" + name_ + ".txt"; //+ (isDirected_? "_dir.txt" : "_undir.txt");

		// The file writer
		FileExport writer = new FileExport(filename, Magnum.set.compressFiles_);

		// Write the header
		//writer.print("node");
		for (int j=0; j<numRefNodes_; j++)
			writer.print("\t" + network_.getRefNode(j).getId());
		writer.print("\n");

		// For each node
		for (int i=0; i<numNodes_; i++) {
			// Node label
			writer.print(network_.getNode(i).getId());

			// Shortest paths
			for (int j=0; j<numRefNodes_; j++)
				writer.print("\t" + MagnumUtils.toStringScientific10(K_.get(i, j)));
			writer.print("\n");
		}
		// Close writer
		writer.close();
		saved_ = true;
	}

	
	// ----------------------------------------------------------------------------

	/** Centrality based on K_ (avg. value for each node, excluding the diagonal) */
	public void computeCentrality() {
		
		centrality_ = new Double[numNodes_];
		
		for (int i=0; i<numNodes_; i++) {
			double sum = 0;
			for (int j=0; j<numNodes_; j++)
				if (i != j)
					sum += K_.get(i, j);
			centrality_[i] = sum / (numNodes_-1);
		}
	}

	
	// ============================================================================
	// PRIVATE METHODS

	/** Normalize the given symmetric matrix (divide by the max, the max is computed on the lower triangular part) */
	protected void normalizeSym(DoubleMatrix2D X) {
		
		if (X.rows() != X.columns())
			throw new IllegalArgumentException("X must be a square, symmetric matrix");

		// Find the max value
		double max = X.get(0, 0);
		for (int i=0; i<X.rows(); i++)
			for (int j=0; j<=i; j++) // Only look at lower triangular part
				if (X.get(i, j) > max)
					max = X.get(i, j);
		
		// Normalize entire matrix
		for (int i=0; i<X.rows(); i++)
			for (int j=0; j<X.columns(); j++)
				X.set(i, j, X.get(i, j)/max);

	}

	// ============================================================================
	// SETTERS AND GETTERS

	public String getName() { return name_; }
	public DoubleMatrix2D getK() { return K_; }
	public Double[] getCentrality() { return centrality_; }
	
}
