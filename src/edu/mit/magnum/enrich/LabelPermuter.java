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

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import edu.mit.magnum.Magnum;
import edu.mit.magnum.gene.Gene;


/**
 * Maps genes to rows of the functional data matrix, provides functionality for label permutation
 */
public class LabelPermuter {

	/** The magnum instance */
	private Magnum mag;

	/** The functional / network data */
	protected FunctionalData functData_ = null;
	/** The genes sorted by centrality */
	private ArrayList<Gene> genes_ = null;
	
	/** The number of bins to use (1 bin = do not correct for centrality) */
	private int numBins_ = -1;
	/** The genes binned according to centrality (degree or average kernel similarity) */
	private ArrayList<ArrayList<Gene>> binnedGenes_ = null;
	/** The corresponding rows / indexes in the functional data matrix */
	private ArrayList<ArrayList<Integer>> binnedIndexes_ = null;
	/** The sum of the centralities of all genes */
	private double centralityVolume_ = -1;

	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor for pairwise funct data */
	public LabelPermuter(Magnum mag, FunctionalData functData, ArrayList<Gene> genes, int numBins) {
		
		this(mag, functData, genes, numBins, -1);
	}

	
	/** Constructor for per gene funct data, specifying which column should be used */
	public LabelPermuter(Magnum mag, FunctionalData functData, ArrayList<Gene> genes, int numBins, int functDataCol) {
		
		this.mag = mag;
		functData_ = functData;
		genes_ = new ArrayList<Gene>(genes);
		numBins_ = numBins;
		
		// Has to be done in this order
		initializeFunctDataIndexes();
		initializeCentrality(functDataCol);
		initializeBins();
	}

	
	// ----------------------------------------------------------------------------

	/** Permute labels of genes within the same bin */
	public void shuffle() {
		
		assert numBins_ == binnedGenes_.size();
		assert numBins_ == binnedIndexes_.size();
		
		// We could shuffle either the indexes or the genes, doesn't matter
		for (int i=0; i<numBins_; i++)
			Collections.shuffle(binnedIndexes_.get(i), mag.set.jdkRng_);
		
		// Assign the new indexes to the genes
		for (int i=0; i<numBins_; i++) {
			ArrayList<Gene> geneBin = binnedGenes_.get(i);
			ArrayList<Integer> indexBin = binnedIndexes_.get(i);
			assert geneBin.size() == indexBin.size();
			
			for (int j=0; j<geneBin.size(); j++)
				geneBin.get(j).setFunctDataIndex(indexBin.get(j));
		}
	}

	
	// ============================================================================
	// PRIVATE METHODS
		
	/** Initialize funct data indexes */
	private void initializeFunctDataIndexes() {
		
		for (Gene gene : genes_) {
			// Get the index of this gwas gene in the functional data matrix 
			Integer index = functData_.getGenes().get(gene.getId());
			if (index == null)
				throw new RuntimeException("Gene not found");
			gene.setFunctDataIndex(index);
		}		
	}

	
	// ----------------------------------------------------------------------------

	/** Initialize centrality and sort genes accordingly */
	private void initializeCentrality(int functDataCol) {
		
		// If it's pairwise data, functDataCol should be -1
		assert (functData_.getIsPairwiseData() || functDataCol >= 0);
		
		centralityVolume_ = 0;

		if (functData_.getIsPairwiseData()) {
			// Avg. of row i
			for (Gene gene_i : genes_) {
				double sum = 0;
				int count = 0;
				for (Gene gene_j : genes_) {
					double w = functData_.get(gene_i.getFunctDataIndex(), gene_j.getFunctDataIndex());
					// Exclude self and NaNs (gene pairs to be excluded)
					if (gene_i != gene_j && !Double.isNaN(w)) {
						sum += w;
						count++;
					}
				}
				// We define the centrality as the sum, not the mean. However, since some elements may have been
				// excluded, we have to compute the mean first and then multiply by the number of genes (-1 to exclude self)
				// Note, centrality is defined in this way because that's what we need in the formulat to compute modularity later on.
				double c = (sum/count) * (genes_.size() - 1);
				if (count == 0)
					c = 0;
				gene_i.setCentrality(c);
			}		
			
		} else {
			for (Gene gene : genes_)
				gene.setCentrality(functData_.get(gene.getFunctDataIndex(), functDataCol));
		}
		
		// Comparator to sort genes by centrality
		final class GeneComparator implements Comparator<Gene> {
			public int compare(Gene g1, Gene g2) {
				return -Double.compare(g1.getCentrality(), g2.getCentrality());
			}
		}
		// Sort the genes
		Collections.sort(genes_, new GeneComparator());
		
		// Compute the sum
		for (Gene gene : genes_)
			centralityVolume_ += gene.getCentrality();
	}

	
	// ----------------------------------------------------------------------------

	/** Bin the genes according to their current order (sort first) */
	private void initializeBins() {
		
		binnedGenes_ = new ArrayList<ArrayList<Gene>>(numBins_);
		binnedIndexes_ = new ArrayList<ArrayList<Integer>>(numBins_);
		
		// Average number of genes per bin
		double numGenesPerBin = genes_.size() / (double) numBins_;
		// Initialize the array lists
		int initCapacity = (int)(numGenesPerBin+1);
		for (int i=0; i<numBins_; i++) {
			binnedGenes_.add(new ArrayList<Gene>(initCapacity));
			binnedIndexes_.add(new ArrayList<Integer>(initCapacity));
		}

		// Bin the gwas genes and initialize with correct funct data index
		for (int i=0; i<genes_.size(); i++) {
			// The current bin
			int bin = (int) Math.floor(i / numGenesPerBin);
			binnedGenes_.get(bin).add(genes_.get(i));
			binnedIndexes_.get(bin).add(genes_.get(i).getFunctDataIndex());
		}		

	}


	// ============================================================================
	// GETTERS AND SETTERS

	public double getCentralityVolume() { return centralityVolume_; }
	
}
