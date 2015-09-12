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

import cern.colt.matrix.DoubleMatrix2D;
import edu.mit.magnum.FileExport;
import edu.mit.magnum.Magnum;
import edu.mit.magnum.MagnumUtils;
import edu.mit.magnum.ProgressMonitor;
import edu.mit.magnum.Settings;
import edu.mit.magnum.gene.Gene;



/**
 * 
 */
abstract public class Enrichment {

	/** The functional data */
	protected DoubleMatrix2D functData_ = null;
	/** Indicates if the data is pairwise */
	protected boolean isPairwiseData_ = false;
	/** The gene scores */
	protected GeneScoreList geneScores_ = null;
	/** The number of genes */
	protected int numGenes_ = -1;
	/** Maps genes to rows of the genePropertyMatrix_, provides functionality for label permutation */
	protected LabelPermuter permuter_ = null;	

	/** Current position in the list of genes */
	protected int currentK_ = -1;
	/** The value of the previous point/sum in the curve before division */
	protected double runningSum_ = -1;
	protected int runningCount_ = -1;
	/** The value of the expected previous point/sum in the curve before division */
	protected double runningSumExpected_ = -1;
	protected int runningCountExpected_ = -1;
	
	/** The indexes for which enrichment is computed */
	protected ArrayList<Integer> k_ = null;

	/** The observed enrichment curve for the original / unpermuted gene list */
	protected Curve curveObs_ = null;	
	/** The observed enrichment curve for the original / unpermuted gene list using sliding windows */
	protected Curve curveObsSlidingWindow_ = null;	
	/** The curve with the gene scores for the original / unpermuted gene list */
	protected Curve curveObsGeneScores_ = null;
	/** The median of the permuted curves at every point */
	protected Curve curveMedian_ = null;
	/** The median of the permuted curves at every point for sliding windows */
	protected Curve curveMedianSlidingWindow_ = null;
	/** The random enrichment curves obtained from permuted gene lists */
	protected ArrayList<Curve> curvesPermut_ = null;
	/** The random enrichment curves obtained from permuted gene lists using sliding windows */
	protected ArrayList<Curve> curvesPermutSlidingWindow_ = null;
	
	/** Used to compute confidence intervals */
	protected EmpiricalPvals empiricalPvals_ = null;
	/** Used to compute confidence intervals for sliding window enrichment */
	protected EmpiricalPvals empiricalPvalsSlidingWindow_ = null;

	/** The area under the curves for the observed and permuted curves */
	protected ArrayList<double[]> AUCs_ = null;
	/** The corresponding empirical p-values */
	double[] pvals_ = null;

	/** The number of random samples / curves */
	protected int numPermutations_ = -1;
	/** Number of random permutations for which enrichment curves are exported (smaller or equal numPermutations) */
	protected int numPermutationsExport_ = -1;
	
	/** The curve that is currently being computed */
	protected Curve curCurve_ = null;
	/** The sliding window curve that is currently being computed */
	protected Curve curCurveSlidingWindow_ = null;

	
	// ============================================================================
	// PUBLIC METHODS

	/** Constructor */
	public Enrichment(FunctionalData functData, GeneScoreList geneScores, LabelPermuter permuter) {

		numPermutations_ = Settings.numPermutations_;
		numPermutationsExport_ = Settings.numPermutationsExport_;
		if (numPermutationsExport_ == -1)
			numPermutationsExport_ = numPermutations_;
		functData_ = functData.getData();
		isPairwiseData_ = functData.getIsPairwiseData();
		geneScores_ = geneScores;
		numGenes_ = geneScores.getNumGenes();
		// Not if mapping to entrez (see Handler() constructor, we print a warning)
		//assert geneScores.getNumGenes() == functData.getNumGenes();
		permuter_ = permuter;

		// k_: the points to be computed / plotted
		initializeK();
				
		// Display info
		if (Settings.curveCutoff_ < 1)
			Magnum.log.println("- " + 100*Settings.curveCutoff_ + "% top genes used");
		Magnum.log.print("- " + k_.size() + " points per curve ");
		if (Settings.varCurveResolution_ > 0) 
			Magnum.log.print("(variable resolution, delta=" + Settings.varCurveResolution_ + ")\n");
		else
			Magnum.log.print("(fixed resolution, delta=" + Settings.constCurveResolution_ + ")\n");
		Magnum.log.println("- " + numPermutations_ + " permutations");
		Magnum.log.println("- " + Settings.numBins_ + " bins for within-degree permutation");
		if (Settings.excludedGenesDistance_ > 0)
			Magnum.log.println("- Excluding gene pairs with respective windows <" + Settings.excludedGenesDistance_ + "mb apart");
		Magnum.log.println("- Sliding window size: " + Settings.slidingWindowSize_);
		Magnum.log.println();

	}

	
	// ----------------------------------------------------------------------------

	/** Run enrichment analysis */
	public void run() {

		// Compute expected
		//expected_ = computeSum();
		
		// Compute enrichment curve for original / unpermuted case (also computes curveExpected_)
		Magnum.log.println("Computing enrichment curve for unpermuted list:");
		long t0 = System.currentTimeMillis();
		
		computeCurve(true);
		curveObs_ = curCurve_;
		curveObsSlidingWindow_ = curCurveSlidingWindow_;
		
		long t1 = System.currentTimeMillis();
		Magnum.log.println("Estimated runtime for " + numPermutations_ + " random permutations: " + MagnumUtils.chronometer(numPermutations_*(t1-t0)));
				
		// Do random permutations
		computePermutCurves();
		
		// Compute empirical p-values based on random permutations
		empiricalPvals_ = new EmpiricalPvals(curvesPermut_, k_);
		empiricalPvals_.computeCurvesSignificance();
		empiricalPvals_.computePvalCurve(curveObs_, false); //true);
		curveMedian_ = empiricalPvals_.getCurveMedian();

		// Compute empirical p-values based on random permutations
		if (Settings.slidingWindowSize_ > 0) {
			empiricalPvalsSlidingWindow_ = new EmpiricalPvals(curvesPermutSlidingWindow_, k_);
			empiricalPvalsSlidingWindow_.computeCurvesSignificance();
			empiricalPvalsSlidingWindow_.computePvalCurve(curveObsSlidingWindow_, false);
			curveMedianSlidingWindow_ = empiricalPvalsSlidingWindow_.getCurveMedian();
		}
		
		// Compute AUC for expected and permut curves
		computeAUC();		
		// Compute p-values
		computePvals();
	}
	
	
	// ----------------------------------------------------------------------------

	/** Reset tagged genes and variables used for curve computation */
	public void reset() {
		
		runningSum_ = 0;
		runningCount_ = 0;
		runningSumExpected_ = 0;
		runningCountExpected_ = 0;
		currentK_ = 0;
	}


	// ----------------------------------------------------------------------------

	/** Save enrichment curves for observed and permuted lists */
	public void save(String filename) {
		
		saveCurves(filename + ".txt");
		saveAUCs(filename + "_AUC.txt");
	}

	
	// ============================================================================
	// ABSTRACT METHODS

	/** Update the running sum with the next gene (currentK_) */
	abstract protected void updateRunningSum();
	
//	/** Update the expected running sum with the next gene (currentK_) */
//	abstract protected void updateExpectedRunningSum();
	
	/** Compute enrichment at the current position (currentK_) */
	abstract protected double computeConnectivity();
	/** Compute enrichment at the current position (currentK_) using sliding window */
	abstract protected double computeSlidingWindowConnectivity();

//	/** Compute the expected enrichment at the current position (currentK_) */
//	abstract protected double computeExpectedEnrichment();
	
	
	// ============================================================================
	// PRIVATE METHODS

	/** Compute the curve for the given ranked list */
	private void computeCurve(boolean isObs) {
		
		// Initialize
		reset();
		curCurve_ = new Curve(k_.size());
		curCurveSlidingWindow_ = new Curve(k_.size());

		ProgressMonitor progress = null;
		if (isObs) {
			curveObsGeneScores_ = new Curve(k_.size());
			progress = new ProgressMonitor(k_.size());
		}
		
		int kIndex = 0; // index pointing to the next k_

		// Walk down the gene list
		for (currentK_=0; currentK_<numGenes_; currentK_++) {			
			// Update the sum
			updateRunningSum();
						
			// Add current point to curve
			if (currentK_ == k_.get(kIndex)) {
				// Overall enrichment
				double enrich = computeConnectivity();
				curCurve_.addPoint(enrich);

				// Window enrichment
				if (Settings.slidingWindowSize_ > 0) {
					double enrichWindow = computeSlidingWindowConnectivity();
					// For the first points, sliding window and overall enrichment is the same (only when doing within window enrichment)
					//assert (currentK_ > Settings.slidingWindowSize_) || 
					//		(Double.isNaN(enrichWindow) && Double.isNaN(enrich)) || 
					//		(Math.abs(enrich-enrichWindow) < 1e-12);
					curCurveSlidingWindow_.addPoint(enrichWindow);
				}					
				
				if (isObs) {
					Gene curGene = geneScores_.getGene(currentK_);
					curveObsGeneScores_.addPoint(curGene.getScore(0));
					progress.iteration(kIndex);
				}
				
				// Check if we're done
				kIndex++;
				if (kIndex == k_.size()) {
					currentK_++; // not sure if needed
					break;
				}
			}
		}
		if (isObs)
			progress.done();
	}

	
	// ----------------------------------------------------------------------------

	/** Compute enrichment for permuted lists */
	private void computePermutCurves() {
		
		curvesPermut_ = new ArrayList<Curve>(numPermutations_);
		curvesPermutSlidingWindow_ = new ArrayList<Curve>(numPermutations_);
		
		//Ngsea.println("Computing enrichment for " + numPermutations_ + " random permutations");
		ProgressMonitor progress = new ProgressMonitor(numPermutations_);

		for (int i=0; i<numPermutations_; i++) {
			// Print progress
			progress.iteration(i);
			// Shuffle and compute curve
			permuter_.shuffle();
			computeCurve(false);
			curvesPermut_.add(curCurve_);
			if (Settings.slidingWindowSize_ > 0)
				curvesPermutSlidingWindow_.add(curCurveSlidingWindow_);
		}
		progress.done();
	}

	
	// ----------------------------------------------------------------------------

	/** Check that curves have consistent number of points */
	private void checkNumPoints() {

		for (int i = 0; i < curvesPermut_.size(); i++)
			if (curvesPermut_.get(i).getNumPoints() != k_.size())
				throw new RuntimeException("Enrichment curve of permuted list "
						+ i + " has inconsistent number of points");

		if (curveObs_.getNumPoints() != k_.size())
			throw new RuntimeException("Enrichment curve of observed list has inconsistent number of points");
		
		// Sliding window enrichment
		if (Settings.slidingWindowSize_ > 0) {
			for (int i = 0; i < curvesPermutSlidingWindow_.size(); i++)
				if (curvesPermutSlidingWindow_.get(i).getNumPoints() != k_.size())
					throw new RuntimeException("Sliding window enrichment curve of permuted list "
							+ i + " has inconsistent number of points");
			
			if (curveObsSlidingWindow_.getNumPoints() != k_.size())
				throw new RuntimeException("Sliding window enrichment curve of observed list has inconsistent number of points");
		}
	}
	
	
	// ----------------------------------------------------------------------------

	/** Initialize the points to be computed / plotted */
	private void initializeK() {

		int numGenes = (int) Math.round(Settings.curveCutoff_ * geneScores_.getNumGenes());
		
		// vary resolution
		if (Settings.varCurveResolution_ > 0) {
			k_ = new ArrayList<Integer>();

			int delta = Settings.varCurveResolution_;
			int nextK = delta-1;
			
			for (int i=0; i<numGenes; i++) {
				if (i == nextK) {
					k_.add(i);
					delta += Settings.varCurveResolution_;
					nextK = i + delta;
				}
			}
			
		// equidistant points based on curvesResolution_
		} else {
			k_ = new ArrayList<Integer>(numGenes/Settings.constCurveResolution_);
			for (int i = 0; i < numGenes; i++) {
				if ((i+1) % Settings.constCurveResolution_ == 0)
					k_.add(i);
			}
		}
		
		// Always add the last point, unless the curve is cut off
		if (k_.get(k_.size()-1) != numGenes-1)
			k_.add(numGenes-1);	
	}

	
	// ----------------------------------------------------------------------------

	/** Compute AUC for expected and permut curves */
	private void computeAUC() {
		
		AUCs_ = new ArrayList<double[]>();
		// Observed
		AUCs_.add(computeAUC(curveObs_));
		// Permutations
		for (int i=0; i<numPermutations_; i++)
			AUCs_.add(computeAUC(curvesPermut_.get(i)));
	}

	
	// ----------------------------------------------------------------------------

	/** Compute empirical p-values */
	private void computePvals() {
		
		// The observed AUCs
		double[] observed = AUCs_.get(0);
		
		// The number of permutations with greater AUCs
		int[] numGreater = new int[observed.length];
		for (int i=0; i<numGreater.length; i++)
			numGreater[i] = 0;

		assert AUCs_.size() == numPermutations_ + 1;
		for (int i=1; i<AUCs_.size(); i++) {
			double[] permut = AUCs_.get(i);
			for (int k=0; k<permut.length; k++)
				if (permut[k] > observed[k])
					numGreater[k]++;
		}
		
		pvals_ = new double[observed.length];
		for (int k=0; k<observed.length; k++)
			pvals_[k] = (numGreater[k] + 1.0) / (numPermutations_ + 1.0);
	}

	
	// ----------------------------------------------------------------------------

	/** Print the empirical p-values */
	public void printPvals() {

		Magnum.log.println("Ranked gene list\tP-value");
		Magnum.log.println(Math.round(100*Settings.curveCutoff_/4.0) + "%\t" + MagnumUtils.toStringScientific10(pvals_[4]));
		Magnum.log.println(Math.round(100*Settings.curveCutoff_/2.0) + "%\t" + MagnumUtils.toStringScientific10(pvals_[5]));
		Magnum.log.println(Math.round(100*3*Settings.curveCutoff_/4.0) + "%\t" + MagnumUtils.toStringScientific10(pvals_[6]));
		Magnum.log.println(Math.round(100*Settings.curveCutoff_) + "%\t" + MagnumUtils.toStringScientific10(pvals_[7]));
		Magnum.log.println();
	}

	

	// ----------------------------------------------------------------------------

	/** Compute AUC for the given curve with respect to the given reference curve */
	private double[] computeAUC(Curve curve) {

		// The first one is linear scale (no enrichment is 1), the second one is log2 scale (no enrichment is 0)
		double[] auc = new double[10];
		for (int i=0; i<10; i++)
			auc[i] = 0;
		
		final double ln_2 = Math.log(2);
		int numGenomeWideSignificant = geneScores_.getNumGenomeWideSignificant();
		
		double[] sumx = new double[5];
		for (int i=0; i<5; i++)
			sumx[i] = 0;
		
		// trapezoidal rule
		for (int i=1; i<k_.size(); i++) {
			if (k_.get(i) < Settings.AUCStart_)
				continue;
			
			double delta_x = k_.get(i) - k_.get(i-1);
			double y1 = curve.getValue(i-1)/curveMedian_.getValue(i-1);
			double y2 = curve.getValue(i)/curveMedian_.getValue(i);
			if (Double.isNaN(y1) || Double.isNaN(y2) || y1 == 0 || y2 ==0 || Double.isInfinite(y1) || Double.isInfinite(y2))
				continue;

			double area = delta_x * (y1 + y2 - 2); // -2 because we integrate with respect to y=1
			double areaLog = delta_x * (Math.log(y1)/ln_2 + Math.log(y2)/ln_2);
			
			int m = 4 - (i / (k_.size()/4));
			for (int j=0; j<m; j++) {
				auc[j] += area;
				auc[4+j] += areaLog;
				sumx[j] += delta_x; 
			}
			if (k_.get(i) <= numGenomeWideSignificant) {
				auc[8] += area;
				auc[9] += areaLog;
				sumx[4] += delta_x;
			}
		}
		for (int i=0; i<4; i++) {
			auc[i] /= 2.0 * sumx[i];
			auc[4+i] /= 2.0 * sumx[i];
		}
		auc[8] /= 2.0 * sumx[4];
		auc[9] /= 2.0 * sumx[4];
		
		return auc;
	}

//	/** Compute AUC for the given curve */
//	private double computeAUC(Curve curve) {
//		
//		double auc = 0;
//		// trapezoidal rule
//		for (int i=1; i<k_.size(); i++) {
//			double delta_x = k_.get(i) - k_.get(i-1);
//			auc += delta_x * ((curve.getValue(i-1) - curveExpected_.getValue(i-1)) + (curve.getValue(i) - curveExpected_.getValue(i))) / 2.0;
//		}
//		return auc;
//	}

	
	// ----------------------------------------------------------------------------

	/** Write observed, expected, FDR and permuted curves */
	private void saveCurves(String filename) {
		
		// Check that curves have consistent number of points
		checkNumPoints();
		// Create output file
		FileExport writer = new FileExport(filename, Settings.compressFiles_);

		// Write header
		writer.print("k");
		writer.print("\tobserved");
		writer.print("\tmedian");
		writer.print("\tp-value");
		
		// Significance curves
		ArrayList<Double> significanceLevels_ = empiricalPvals_.getSignificanceLevels();
		for (int i=0; i<significanceLevels_.size(); i++)
			writer.print("\t" + "p_" + significanceLevels_.get(i));
		for (int i=significanceLevels_.size()-1; i>=0; i--)
			writer.print("\t" + "p_" + (1-significanceLevels_.get(i)));

		writer.print("\tgeneScore");

		// Sliding window
		if (Settings.slidingWindowSize_ > 0) {
			writer.print("\tobservedWindow");
			writer.print("\tmedianWindow");
			writer.print("\tp-valueWindow");
		
			// Significance curves
			ArrayList<Double> significanceLevelsWindow_ = empiricalPvalsSlidingWindow_.getSignificanceLevels();
			for (int i=0; i<significanceLevelsWindow_.size(); i++)
				writer.print("\t" + "pWindow_" + significanceLevelsWindow_.get(i));
			for (int i=significanceLevelsWindow_.size()-1; i>=0; i--)
				writer.print("\t" + "pWindow_" + (1-significanceLevelsWindow_.get(i)));
		}
		
		// Curves for permuted lists
		for (int i=0; i<numPermutationsExport_; i++)
			writer.print("\t" + "permut" + (i+1));			
		writer.print("\n");		
		
		for (int p=0; p<k_.size(); p++) {
			// k
			writer.print(Integer.toString(k_.get(p)+1));

			// Enrichment curves
			writer.print("\t" + MagnumUtils.toStringScientific10(curveObs_.getValue(p)));
			writer.print("\t" + MagnumUtils.toStringScientific10(curveMedian_.getValue(p)));
			writer.print("\t" + MagnumUtils.toStringScientific10(empiricalPvals_.getCurvePval().getValue(p)));			
			for (int i=0; i<empiricalPvals_.getCurvesSignificance().size(); i++)	
				writer.print("\t" + MagnumUtils.toStringScientific10(empiricalPvals_.getCurvesSignificance().get(i).getValue(p)));

			// Gene score
			writer.print("\t" + MagnumUtils.toStringScientific10(curveObsGeneScores_.getValue(p)));
			
			// Sliding window enrichment curves
			if (Settings.slidingWindowSize_ > 0) {
				writer.print("\t" + MagnumUtils.toStringScientific10(curveObsSlidingWindow_.getValue(p)));
				writer.print("\t" + MagnumUtils.toStringScientific10(curveMedianSlidingWindow_.getValue(p)));
				writer.print("\t" + MagnumUtils.toStringScientific10(empiricalPvalsSlidingWindow_.getCurvePval().getValue(p)));			
				for (int i=0; i<empiricalPvalsSlidingWindow_.getCurvesSignificance().size(); i++)	
					writer.print("\t" + MagnumUtils.toStringScientific10(empiricalPvalsSlidingWindow_.getCurvesSignificance().get(i).getValue(p)));
			}
			
			// Curves for permuted lists
			for (int i=0; i<numPermutationsExport_; i++)
				writer.print("\t" + MagnumUtils.toStringScientific10(curvesPermut_.get(i).getValue(p)));

			writer.print("\n");
		}
		writer.close();
	}

	
	// ----------------------------------------------------------------------------

	/** Write the AUCs */
	private void saveAUCs(String filename) {
		
		FileExport writer = new FileExport(filename, Settings.compressFiles_);
		writer.println("AUC_1\tAUC_2\tAUC_3\tAUC_4\tAUClog2_1\tAUClog2_2\tAUClog2_3\tAUClog2_4\tAUC_GWS\tAUClog2_GWS");
		for (int i=0; i<AUCs_.size(); i++) {
			String line = MagnumUtils.toStringScientific10(AUCs_.get(i)[0]);
			for (int j=1; j<10; j++)
				line += "\t" + MagnumUtils.toStringScientific10(AUCs_.get(i)[j]);
			writer.println(line);
		}
		writer.close();
	}

		
	// ============================================================================
	// GETTERS AND SETTERS
	
	public Curve getCurveObs() { return curveObs_; }
//	public Curve getCurveExpected() { return curveExpected_; }
	//public double getExpected() { return expected_; }
	public ArrayList<double[]> getAUCs() { return AUCs_; }
}
