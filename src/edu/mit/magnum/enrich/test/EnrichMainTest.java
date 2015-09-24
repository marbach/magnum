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
package edu.mit.magnum.enrich.test;

import java.io.File;
import java.util.ArrayList;

import static org.junit.Assert.*;

import org.junit.*;

import edu.mit.magnum.Magnum;
import edu.mit.magnum.enrich.*;
import edu.mit.magnum.gene.*;

/**
 * Unit tests for EnrichmentTest
 */
public class EnrichMainTest {
	
	/** The magnum instance */
	private static Magnum mag = new Magnum();

	// ============================================================================
	// SETUP
	
	@BeforeClass
	public static void testSetup() {
		
		mag.set.resetToDefaults();
		mag.set.geneScoreFile_ = new File("src/edu/mit/magnum/enrich/test/simpleNet_genescores.txt");
		mag.set.constCurveResolution_ = 1;
		mag.set.varCurveResolution_ = -1;
		mag.set.curveCutoff_ = 1;
		mag.set.functionalDataCols_ = null;
		mag.set.numPermutations_ = 10;
		mag.set.numPermutationsExport_ = 10;
		mag.set.numBins_ = 1;
		mag.set.setRandomSeed(1);
		mag.set.pval_ = new ArrayList<Double>();
		mag.set.pval_.add(0.2);
		mag.set.geneCoordFile_ = new File("src/edu/mit/magnum/enrich/test/simpleNet_geneCoords.bed");
		mag.set.idTypeFunctionalData_ = "custom";
		mag.set.idTypeGeneScores_ = "custom";
		mag.set.excludeGenomeWideSignificantGenes_ = false;
	}

	@AfterClass
	public static void testCleanup() {
	}
	  
	// ============================================================================
	// TESTS

	/** Loading excludedGenes and excludedGeneParis files */
	@Test
	public void testPairwiseEnrichment_loadExcludedFromFiles() {

		// Settings
		mag.set.functionalDataFile_ = new File("src/edu/mit/magnum/enrich/test/simpleNet_testKernel.txt");
		mag.set.excludedGenesFile_ = new File("src/edu/mit/magnum/enrich/test/simpleNet_excludedGenes.txt");
		mag.set.excludedGenePairsFile_ = new File("src/edu/mit/magnum/enrich/test/simpleNet_excludedGenePairs.txt");
		mag.set.excludeXYChromosomes_ = false;
		mag.set.excludedGenesDistance_ = -1;
		
		// Run enrichment analysis
		EnrichMain enrichMain = new EnrichMain(mag);
		enrichMain.run();

		// TEST INITIALIZATION OF CENTRALITIES BY PERMUTER
		double[] c = { -1, 0.2650766, 0.3608349, 0.3998752, 0.3803405, 0.2583860, 0.3477147 };
	
		ArrayList<Gene> genes = enrichMain.getGeneScores().getGenes();
		double epsilon = 1e-6;
		assertEquals(genes.size(), 6);
		
		assertEquals(genes.get(0).id_, "6");
		assertEquals(genes.get(0).getCentrality(), 5*c[6], epsilon);
		
		assertEquals(genes.get(1).id_, "5");
		assertEquals(genes.get(1).getCentrality(), 5*c[5], epsilon);

		assertEquals(genes.get(2).id_, "4");
		assertEquals(genes.get(2).getCentrality(), 5*c[4], epsilon);

		assertEquals(genes.get(3).id_, "2");
		assertEquals(genes.get(3).getCentrality(), 5*c[2], epsilon);

		assertEquals(genes.get(4).id_, "3");
		assertEquals(genes.get(4).getCentrality(), 5*c[3], epsilon);

		assertEquals(genes.get(5).id_, "1");
		assertEquals(genes.get(5).getCentrality(), 5*c[1], epsilon);
 
		//assertEquals(genes.get(6).id_, "7");
		//assertEquals(genes.get(6).getCentrality(), 5*c[7], epsilon);
		
//		// TEST EXPECTED
//		double[] e = { 0.0000000, 0.4564072, 0.4823739, 0.4021919, 0.3942546, 0.3418840 };
//		Curve exp = handler.getEnrichment().getCurveExpected();
//		assertEquals(exp.getNumPoints(), 6);
//		for (int i=0; i<e.length; i++)
//			assertEquals(e[i], exp.getValue(i), epsilon);
//		//assertEquals(e[e.length-1], handler.getEnrichment().getExpected(), epsilon);
		
		// TEST OBSERVED
		Curve obs = enrichMain.getEnrichment().getCurveObs();
		assertEquals(obs.getNumPoints(), 6);
		// Expected distances
		double[] o = { 0.0000000, 0.7167235, 0.7108765, 0.3574122, 0.3762026, 0.3365148 };
		// Check distances
		for (int i=0; i<o.length; i++)
			assertEquals(o[i], obs.getValue(i), epsilon);
		
//		// TEST AUC
//		ArrayList<Double> auc = handler.getEnrichment().getAUCs();
//		assertEquals(0.4233026, auc.get(0), epsilon);
	}

	
	/** Excluding allosomes and gene pairs based on loaded coords */
	@Test
	public void testPairwiseEnrichment_excludeBasedOnCoords() {

		// Settings
		mag.set.functionalDataFile_ = new File("src/edu/mit/magnum/enrich/test/simpleNet_testKernel.txt");
		mag.set.excludedGenesFile_ = null;
		mag.set.excludedGenePairsFile_ = null;
		mag.set.excludeXYChromosomes_ = true;
		mag.set.excludedGenesDistance_ = 1;
		
		// Run enrichment analysis
		EnrichMain enrichMain = new EnrichMain(mag);
		enrichMain.run();

		// TEST INITIALIZATION OF CENTRALITIES BY PERMUTER
		double[] c = { -1, 0.2650766, 0.3608349, 0.3998752, 0.3803405, 0.2583860, 0.3477147 };
	
		ArrayList<Gene> genes = enrichMain.getGeneScores().getGenes();
		double epsilon = 1e-6;
		assertEquals(genes.size(), 6);
		
		assertEquals(genes.get(0).id_, "6");
		assertEquals(genes.get(0).getCentrality(), 5*c[6], epsilon);
		
		assertEquals(genes.get(1).id_, "5");
		assertEquals(genes.get(1).getCentrality(), 5*c[5], epsilon);

		assertEquals(genes.get(2).id_, "4");
		assertEquals(genes.get(2).getCentrality(), 5*c[4], epsilon);

		assertEquals(genes.get(3).id_, "2");
		assertEquals(genes.get(3).getCentrality(), 5*c[2], epsilon);

		assertEquals(genes.get(4).id_, "3");
		assertEquals(genes.get(4).getCentrality(), 5*c[3], epsilon);

		assertEquals(genes.get(5).id_, "1");
		assertEquals(genes.get(5).getCentrality(), 5*c[1], epsilon);
 		
		// TEST OBSERVED
		Curve obs = enrichMain.getEnrichment().getCurveObs();
		assertEquals(obs.getNumPoints(), 6);
		// Expected distances
		double[] o = { 0.0000000, 0.7167235, 0.7108765, 0.3574122, 0.3762026, 0.3365148 };
		// Check distances
		for (int i=0; i<o.length; i++)
			assertEquals(o[i], obs.getValue(i), epsilon);
	}

	
	/** Compute kernel instead of loading it from file */
	@Test
	public void testPairwiseEnrichment_computeKernel() {

		// Settings
		mag.set.functionalDataFile_ = null;
		mag.set.networkDir_ = new File(".");
		mag.set.networkFile_ = new File("src/edu/mit/magnum/netprop/test/simpleNet.txt");
		mag.set.isDirected_ = false;
		mag.set.isWeighted_ = false;
		mag.set.excludedGenesFile_ = null;
		mag.set.excludedGenePairsFile_ = null;
		mag.set.excludeXYChromosomes_ = true;
		mag.set.excludedGenesDistance_ = 1;
		mag.set.computePstepKernel_ = true;
		mag.set.pstepKernelP_ = new ArrayList<Integer>();
		mag.set.pstepKernelP_.add(4);
		mag.set.pstepKernelAlpha_ = 2;

		// Run enrichment analysis
		EnrichMain enrichMain = new EnrichMain(mag);
		enrichMain.run();

		// TEST INITIALIZATION OF CENTRALITIES BY PERMUTER
		double[] c = { -1, 0.2650766, 0.3608349, 0.3998752, 0.3803405, 0.2583860, 0.3477147 };
	
		ArrayList<Gene> genes = enrichMain.getGeneScores().getGenes();
		double epsilon = 1e-6;
		assertEquals(genes.size(), 6);
		
		assertEquals(genes.get(0).id_, "6");
		assertEquals(genes.get(0).getCentrality(), 5*c[6], epsilon);
		
		assertEquals(genes.get(1).id_, "5");
		assertEquals(genes.get(1).getCentrality(), 5*c[5], epsilon);

		assertEquals(genes.get(2).id_, "4");
		assertEquals(genes.get(2).getCentrality(), 5*c[4], epsilon);

		assertEquals(genes.get(3).id_, "2");
		assertEquals(genes.get(3).getCentrality(), 5*c[2], epsilon);

		assertEquals(genes.get(4).id_, "3");
		assertEquals(genes.get(4).getCentrality(), 5*c[3], epsilon);

		assertEquals(genes.get(5).id_, "1");
		assertEquals(genes.get(5).getCentrality(), 5*c[1], epsilon);
 		
		// TEST OBSERVED
		Curve obs = enrichMain.getEnrichment().getCurveObs();
		assertEquals(obs.getNumPoints(), 6);
		// Expected distances
		double[] o = { 0.0000000, 0.7167235, 0.7108765, 0.3574122, 0.3762026, 0.3365148 };
		// Check distances
		for (int i=0; i<o.length; i++)
			assertEquals(o[i], obs.getValue(i), epsilon);
	}

	
	// ----------------------------------------------------------------------------

//	/** Test enrichment for per gene average kernel similarity */
//	@Test
//	public void testIndividualEnrichment() {
//
//		// Settings
//		Settings.functionalDataFile_ = "data/test/simpleNet_testNodeProperties.txt";
//		
//		// Run enrichment analysis
//		Handler handler = new Handler();
//		handler.run();
//		
//		// TEST INITIALIZATION OF CENTRALITIES BY PERMUTER
//		HashMap<String, Gene> genes = handler.getGwas().getGenes();
//		double epsilon = 1e-6;
//		assertEquals(genes.size(), 6);
//		assertEquals(genes.get("1").getCentrality(), 0.26507663571757084, epsilon);
//		assertEquals(genes.get("2").getCentrality(), 0.36083486759002886, epsilon);
//		assertEquals(genes.get("3").getCentrality(), 0.39987517924678223, epsilon);
//		assertEquals(genes.get("4").getCentrality(), 0.44527829169338745, epsilon);
//		assertEquals(genes.get("5").getCentrality(), 0.3477147181566287, epsilon);
//		assertEquals(genes.get("6").getCentrality(), 0.3477147181566287, epsilon);
//
//		// TEST EXPECTED
//		Curve exp = handler.getEnrichment().getCurveExpected();
//		for (int i=0; i<exp.getNumPoints(); i++)
//			assertEquals(0.3610824, exp.getValue(i), epsilon);
//
//		Curve obs = handler.getEnrichment().getCurveObs();
//		
//		// Expected distances
//		double[] o = { 0.3477147, 0.3477147, 0.3802359, 0.3802359, 0.3753856, 0.3753856, 0.3753856, 0.3802836, 0.3610824 };
//		
//		// Check distances
//		for (int i=0; i<o.length; i++)
//			assertEquals(o[i], obs.getValue(i), epsilon);
//	}
//

	// ============================================================================
	// PRIVATE METHODS

}
