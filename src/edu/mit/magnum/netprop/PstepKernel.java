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
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.Blas;
import cern.colt.matrix.linalg.SeqBlas;
import edu.mit.magnum.Magnum;
import edu.mit.magnum.MagnumUtils;
import edu.mit.magnum.net.*;


/**
 * Compute p-step random walk kernel (Smola & Kondor, 2003)
 * K = (a*I − L)^p ,  with a >= 2 
 * Where L is the normalized Laplacian: L = I - D^(-1/2)*A*D^(-1/2)
 * D is the nxn diagonal degree matrix with D_ii = Sum_j(A_ij)
 *     
 * Note that the network has to be undirected and self-loops are ignored
 */
public class PstepKernel extends PairwiseProperties {
		
	/** The alpha parameter (must be >= 2) */
	private double alpha_ = 2;
	/** Steps p of random walk kernel (ordered list of positive integers given in increasing order, the kernel for each listed p will be saved) */
	private ArrayList<Integer> p_ = null;
	/** The total number of steps (equal to the last element of p_) */
	private int numSteps_ = -1;
	/** Normalize the kernel matrix (divide by the max) */
	private boolean normalize_ = true;
	
	/** The normalized laplacian */
	private SparseDoubleMatrix2D normalizedLaplacian_ = null;

	/** Node centrality for each alpha / step */
	private LinkedHashMap<String,Double[]> pstepCentrality_ = null;

	/** Colt linear algebra methods */
	private final Algebra alg_ = new Algebra();
	/** Colt basic linear algebra system */
	private final Blas blas_ = SeqBlas.seqBlas;

	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor (network must be undirected) */
	public PstepKernel(Network network, double alpha, ArrayList<Integer> p, boolean normalize, boolean computeCentrality) {
		
		super(network, "pstepKernel", "pstepKernel", computeCentrality);
		alpha_ = alpha;
		p_ = p;
		numSteps_ = p_.get(p_.size()-1);
		normalize_ = true;
		pstepCentrality_ = new LinkedHashMap<String,Double[]>();
		
		if (alpha_ < 2)
			throw new IllegalArgumentException("Alpha must be greater or equal 2");
		
		if (!MagnumUtils.posIntIncreasing(p_))
			throw new RuntimeException("p must be and ordered list of positive integers, given in increasing order");
		
		if (isDirected_)
			throw new IllegalArgumentException("P-step kernels are not implemented for directed networks");
		if (numRefNodes_ != numNodes_)
			Magnum.log.warning("Specified reference nodes will be ignored by p-step kernel");
		// We now explicitly check for self-loops and abort if there are, because it screws up the degrees
		//if (!Settings.removeSelfLoops_)
			//Ngsea.warning("Self-loops will be ignored by p-step kernel");
	}
	
	
	/** Constructor using default alpha (network must be undirected) */
	public PstepKernel(Network network, ArrayList<Integer> p) {
		
		this(network, Magnum.set.pstepKernelAlpha_, p, Magnum.set.pstepKernelNormalize_, Magnum.set.exportNodeProperties_);
	}


	// ----------------------------------------------------------------------------

	/** Compute the p-step kernel matrix. Matrix multiplication could be done more efficiently by exploiting symmetry. */
	public void computeK() {

		Magnum.log.println("Computing normalized Laplacian...");		
		normalizedLaplacian_ = network_.computeNormalizedLaplacian();
		
		Magnum.log.println("Computing " + numSteps_ + "-step kernel with alpha=" + alpha_ + ":");
		Magnum.log.println("Step 1...");

		// K = (a*I - L)^p ,  with a >= 2
		// B := a*I - L
		SparseDoubleMatrix2D B = new SparseDoubleMatrix2D(numNodes_, numNodes_); // initializes at 0
		for (int i=0; i<numNodes_; i++)
			B.set(i, i, alpha_);
				
		// daxpy(double alpha, DoubleMatrix2D A, DoubleMatrix2D B)
        // Combined matrix scaling; B = B + alpha*A.
		blas_.daxpy(-1, normalizedLaplacian_, B);	
		
		// Using pow() or blas_.dggm() is much slower than the implementation below, actually
		// pow() uses blas. This is surprising, because pow() shows my implementation as the
		// naive approach and they do something much more sophisticated with blas.
		// Maybe the problem is that blas doesn't leverage SparseDoubleMatrix? Weird.
		//K_ = alg_.pow(B, p_);

		// K = B^1
		// Wow, the multiplication below is much faster if this is a dense matrix, probably access is faster
		K_ = new DenseDoubleMatrix2D(numNodes_, numNodes_);
		blas_.dcopy(B, K_);

		saved_ = false;
		saveStep(1, alpha_);

		// K = B^p
		for (int i=2; i<=p_.get(p_.size()-1); i++) {
			Magnum.log.println("Step " + i + "...");
			saved_ = false;

			//long t0 = System.currentTimeMillis();

			// mult() just calls DoubleMatrix2D.zMult(), which is implemented differently by Sparse and Dense matrices
			// The sparse and dense matrix implementations are optimized and too difficult to understand/modify
			// for symmetric matrixes.
			// Note, it seems runtime is fastest if one matrix is sparse and the other dense
			K_ = alg_.mult(B, K_);
			//K_ = B.zMult(K_, null);
			
			//long t1 = System.currentTimeMillis();
			//Ngsea.println("Run time: " + NgseaUtils.chronometer(t1-t0));

			// Save step, also computes centrality
			saveStep(i, alpha_);
		}		
		
		// Delete Laplacian
		normalizedLaplacian_ = null;
	}


	// ----------------------------------------------------------------------------

	public void addNodeProperties(LinkedHashMap<String,Number[]> map) {
		
		if (pstepCentrality_.size() != p_.size())
			throw new RuntimeException("Centrality has not been saved for every specified step");
		
		for (String centrality : pstepCentrality_.keySet())
			map.put(centrality, pstepCentrality_.get(centrality));
	}

	
	// ============================================================================
	// PRIVATE METHODS

	/** Save step i if it is specified in p_ (also normalizes K_ if normalize_ is set) */
	private void saveStep(int i, double alpha) {
		
		// Return if this step should not be saved
		if (!p_.contains(i))
			return;
		
		if (normalize_)
			normalizeSym(K_);
		
		String suffix = "_alpha" + alpha + (network_.getIsWeighted() ? "_weighted" : "");
		name_ = i + "stepKernel" + suffix;
		nameCentrality_ = i + "stepKernelCentrality" + suffix;
		
		if (Magnum.set.exportPairwiseNodeProperties_)
			saveK(MagnumUtils.extractBasicFilename(network_.getFile().getName(), false));
		
		if (computeCentrality_) {
			computeCentrality();
			pstepCentrality_.put(nameCentrality_, centrality_);
		}
	}


	// ============================================================================
	// SETTERS AND GETTERS

	public DoubleMatrix2D getK() { return K_; }
	public Double[] getPstepKernelCentrality() { return centrality_; }
	
}
