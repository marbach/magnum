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

import edu.mit.magnum.Magnum;


/**
 * 
 */
public class EmpiricalPvals {

	/** The indexes for which enrichment is computed */
	private ArrayList<Integer> k_ = null;
	/** The number of permutations that were done */
	private int numPermut_ = -1;
	/** The number of points of each curve */
	private int numPoints_ = -1;
	
	/** Used for pval computation at every point of a curve */
	private ArrayList<Percentile> percentiles_ = null;

	/** The significance values for which boundary curves should be drawn */
	private ArrayList<Double> significanceLevels_ = null;
	/** The boundary curves at the given significance levels */
	private ArrayList<Curve> curvesSignificance_ = null;
	/** The median curve */
	private Curve curveMedian_ = null;

	/** The empirical p-value of curveObs_ at every point */
	private Curve curvePval_ = null;
	/** The most significant p-value */
	private double minPval_ = -1;
	/** The cutoff (k_) at which the most significant p-value occurs */
	private int minPvalK_ = -1;
	/** Indicates whether it's enrichment or depletion */
	private boolean minPvalIsEnrichment_ = false;
	
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	EmpiricalPvals(ArrayList<Curve> curvesPermut, ArrayList<Integer> k) {
		
		k_ = k;
		initialize(curvesPermut);
		//initializeFdr();
	}


	// ----------------------------------------------------------------------------

	/** Compute curves for the given significance levels */
	public void computeCurvesSignificance() {
		
		// The median curve
		curveMedian_ = new Curve(numPoints_);
		// The boundary curves for the given pvals
		curvesSignificance_ = new ArrayList<Curve>();
		for (int i=0; i<2*significanceLevels_.size(); i++)
			curvesSignificance_.add(new Curve(numPoints_));
		
		// Compute the indexes corresponding to the given percentiles
		ArrayList<Integer> indexes = new ArrayList<Integer>();
		for (int i=0; i<significanceLevels_.size(); i++)
			indexes.add((int) Math.round(significanceLevels_.get(i)*numPermut_));

		// Add the corresponding points to the boundary curves
		for (int p=0; p<numPoints_; p++) {
			ArrayList<Double> pctiles = percentiles_.get(p).getValues(indexes);
			curveMedian_.addPoint(pctiles.get(0));
			for (int i=1; i<pctiles.size(); i++)
				curvesSignificance_.get(i-1).addPoint(pctiles.get(i));
		}
	}
	
	
	// ----------------------------------------------------------------------------

	/** Compute pval curve and min pval */
	public void computePvalCurve(Curve curveObs, boolean printInfo) {
		
		curvePval_ = new Curve(numPoints_);
		
		// Compute p-value for observed
		for (int p=0; p<numPoints_; p++) {
			double obs = curveObs.getValue(p);
			double pval = percentiles_.get(p).pValue(obs);
			curvePval_.addPoint(pval);
		}
		
		// Compute most significant p-value, the sign indicates if it's higher/lower than expected
		minPval_ = 1;
		minPvalK_ = -1;
		minPvalIsEnrichment_ = false;
		
		for (int p=0; p<curvePval_.getNumPoints(); p++) {
			if (k_.get(p) < Magnum.set.AUCStart_)
				continue;
			
			double pval = curvePval_.getValue(p);
			if (pval > 0.5)
				pval = 1 - pval;
			
			if (pval < minPval_) {
				minPval_ = pval;
				minPvalK_ = k_.get(p);
				minPvalIsEnrichment_ = curvePval_.getValue(p) > 0.5;
			}
		}
		
		if (printInfo) {
			Magnum.log.println("Most significant point:");
			Magnum.log.println("- " + (minPvalIsEnrichment_ ? "ENRICHMENT" : "DEPLETION"));
			Magnum.log.println("- k    = " + minPvalK_);
			Magnum.log.println("- pval = " + minPval_);
			Magnum.log.println();
		}
	}

	
		
	// ============================================================================
	// PRIVATE METHODS

	/** Initialize */
	private void initialize(ArrayList<Curve> curvesPermut) {
		
		// numPermut_, numPoints_
		numPermut_ = curvesPermut.size();
		numPoints_ = k_.size();
		
		// significanceLevels_ (copy because we modify below)
		significanceLevels_ = new ArrayList<Double>(Magnum.set.pval_);
		// Check that significance levels are below 0.5
		for (int i=0; i<significanceLevels_.size(); i++)
			if (significanceLevels_.get(i) >= 0.5)
				throw new IllegalArgumentException("Significance levels must not be <0.5 (found: " + significanceLevels_.get(i) + ")");

		// Divide by two if two-sided test
		if (Magnum.set.twoSidedTest_)
			for (int i=0; i<significanceLevels_.size(); i++)
				significanceLevels_.set(i, significanceLevels_.get(i)/2.0);

		// percentiles_
		percentiles_ = new ArrayList<Percentile>(numPoints_);
		for (int p=0; p<numPoints_; p++) {
			// Create the percentiles for this position in the list
			Percentile pctile = new Percentile();
			percentiles_.add(pctile);
			
			// Add all points/curves at this position in the list
			for (int i=0; i<numPermut_; i++)
				pctile.addPoint(curvesPermut.get(i).getPoint(p));
		}
	}
	
	
	// ----------------------------------------------------------------------------

//	/** Initialize fdr_ */
//	private void initializeFdr() {
//		
//		HashSet<Curve> falsePositives = new HashSet<Curve>();
//		fdr_ = new double[numPermut_/2];
//		
//		for (int i=0; i<numPermut_/2; i++) {
//			for (int p=0; p<numPoints_; p++) {
//				percentiles_.get(index)
//			}
//		}
//	}

	// ============================================================================
	// GETTERS AND SETTERS

	public Curve getCurveMedian() { return curveMedian_; }
	public ArrayList<Double> getSignificanceLevels() { return significanceLevels_; }
	public ArrayList<Curve> getCurvesSignificance() { return curvesSignificance_; }
	public Curve getCurvePval() { return curvePval_; }
	public double getMinPval() { return minPval_; }
	public int getMinPvalK() { return minPvalK_; }
	public boolean getMinPvalIsEnrichment() { return minPvalIsEnrichment_; }

}
