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
package edu.mit.magnum.enrich;

import edu.mit.magnum.Settings;
import edu.mit.magnum.gene.Gene;


/**
 * 
 */
public class EnrichmentPairwise extends Enrichment {

	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public EnrichmentPairwise(FunctionalData functData, GeneScoreList geneScores, LabelPermuter permuter) {

		super(functData, geneScores, permuter);
		
		// Check that functional data indexes of genes are within bounds of functional data
		// so that we can use getQuick() later on
		int N = functData_.rows();
		if (N != functData_.columns())
			throw new RuntimeException("Expected square matrix");
		for (Gene gene : geneScores_.getGenes())
			if (gene.getFunctDataIndex() >= N)
				throw new RuntimeException("Functional data index out of bounds");
	}


	// ============================================================================
	// PRIVATE METHODS

	/** Update the running sum with the given gene */
	protected void updateRunningSum() {

		int curGeneIndex = geneScores_.getGene(currentK_).getFunctDataIndex();
		
		// For all previous genes
		for (int i=0; i<currentK_; i++) {
			int prevGeneIndex = geneScores_.getGene(i).getFunctDataIndex();
			double w = functData_.getQuick(curGeneIndex, prevGeneIndex);
			
			if (!Double.isNaN(w)) {
				runningSum_ += w;
				runningCount_++;
			}
		}
	}
	
	
	// ----------------------------------------------------------------------------

	/** Compute connectivity at the current position (currentK_) */
	protected double computeConnectivity() {
		
		int N = currentK_ + 1;
		if (N <= 1)
			return 0;

		// 2 * ... because we only summed the upper triangular part
		//return 2 * runningSum_ / (N * (N-1));
		return runningSum_ / runningCount_;
	}

	
	// ----------------------------------------------------------------------------

	/** Compute connectivity at the current position (currentK_) */
	protected double computeSlidingWindowConnectivity() {
		
		int N = currentK_ + 1;
		if (N <= 1)
			return 0;

		int windowStart = Math.max(0, currentK_-Settings.slidingWindowSize_+1);
		double sum = 0;
		int count = 0;
		
		for (int i=windowStart; i<N; i++) {
			int gene_i = geneScores_.getGene(i).getFunctDataIndex();
			
			//for (int j=i+1; j<N; j++) {  // <-- connectivity within sliding window
			for (int j=0; j<i; j++) {    // <-- connectivity with all previous genes
				int gene_j = geneScores_.getGene(j).getFunctDataIndex();
				
				double w = functData_.getQuick(gene_i, gene_j);
				if (!Double.isNaN(w)) {
					sum += w;
					count++;
				}
			}
		}
		return sum / count;
	}

	
}
