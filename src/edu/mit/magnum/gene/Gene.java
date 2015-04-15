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
package edu.mit.magnum.gene;


/**
 * A gene and it's associated SNPs.
 * Note, equals() and hashcode() are based on Node.id_.
 */
public class Gene extends GenomicElement {

	/** The gene name / symbol (use GenomicElement.id_ for the ensembl/entrez ID) */
	public String symbol_ = null;
	
	/** The index of this gene in the gene functional data matrix */
	private int functDataIndex_ = -1;
	/** The centrality of the gene defined by the functional data (e.g., avg. similarity) */
	private double centrality_ = 0;
	/** The association scores / p-values */
	private double[] score_ = null;

	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public Gene(String id) {

		super(id);
	}

	
	/** Constructor */
	public Gene(String id, String symbol) {
		
		super(id);
		symbol_ = symbol;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Get the TSS (start if on + strand, end if on - strand) */
	public int getTss() {
		
		if (posStrand_)
			return start_;
		else
			return end_;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Get the the symbol or "NA" if the symbol is NULL */
	public String getSymbolOrNA() {
		
		return (symbol_ == null) ? "NA" : symbol_;
	}

		
	// ============================================================================
	// GETTERS AND SETTERS
		
	public void setFunctDataIndex(int i) { functDataIndex_ = i; }
	public int getFunctDataIndex() { return functDataIndex_; }
	
	public void setCentrality(double c) { centrality_ = c; }
	public double getCentrality() { return centrality_; }
	
	public void setScore(double[] x) { score_ = x; }
	public double[] getScore() { return score_; }
	public double getScore(int i) { return score_[i]; }
	
}
