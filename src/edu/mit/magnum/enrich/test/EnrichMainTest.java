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

import java.util.ArrayList;

import static org.junit.Assert.*;

import org.junit.*;

import edu.mit.magnum.Settings;
import edu.mit.magnum.enrich.*;
import edu.mit.magnum.gene.*;

/**
 * Unit tests for EnrichmentTest
 */
public class EnrichMainTest {
	
	
	// ============================================================================
	// SETUP
	
	@BeforeClass
	public static void testSetup() {
		
		Settings.loadSettings();		
		Settings.geneScoreFile_ = "src/edu/mit/magnum/enrich/test/simpleNet_genescores.txt";
		Settings.constCurveResolution_ = 1;
		Settings.varCurveResolution_ = -1;
		Settings.curveCutoff_ = 1;
		Settings.functionalDataCols_ = null;
		Settings.numPermutations_ = 10;
		Settings.numPermutationsExport_ = 10;
		Settings.numBins_ = 1;
		Settings.randomSeed_ = 1;
		Settings.pval_ = new ArrayList<Double>();
		Settings.pval_.add(0.2);
		Settings.geneCoordFile_ = "src/edu/mit/magnum/enrich/test/simpleNet_geneCoords.bed";
		Settings.idTypeFunctionalData_ = "custom";
		Settings.idTypeGeneScores_ = "custom";
		Settings.excludeGenomeWideSignificantGenes_ = false;
		Settings.initialize();
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
		Settings.functionalDataFile_ = "src/edu/mit/magnum/enrich/test/simpleNet_testKernel.txt";
		Settings.excludedGenesFile_ = "src/edu/mit/magnum/enrich/test/simpleNet_excludedGenes.txt";
		Settings.excludedGenePairsFile_ = "src/edu/mit/magnum/enrich/test/simpleNet_excludedGenePairs.txt";
		Settings.ignoreAllosomes_ = false;
		Settings.excludedGenesDistance_ = -1;
		
		// Run enrichment analysis
		EnrichMain enrichMain = new EnrichMain();
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
		Settings.functionalDataFile_ = "src/edu/mit/magnum/enrich/test/simpleNet_testKernel.txt";
		Settings.excludedGenesFile_ = "";
		Settings.excludedGenePairsFile_ = "";
		Settings.ignoreAllosomes_ = true;
		Settings.excludedGenesDistance_ = 1;
		
		// Run enrichment analysis
		EnrichMain enrichMain = new EnrichMain();
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
		Settings.functionalDataFile_ = "";
		Settings.networkDir_ = ".";
		Settings.networkFile_ = "src/edu/mit/magnum/netprop/test/simpleNet.txt";
		Settings.isDirected_ = false;
		Settings.isWeighted_ = false;
		Settings.excludedGenesFile_ = "";
		Settings.excludedGenePairsFile_ = "";
		Settings.ignoreAllosomes_ = true;
		Settings.excludedGenesDistance_ = 1;
		Settings.computePstepKernel_ = true;
		Settings.pstepKernelP_ = new ArrayList<Integer>();
		Settings.pstepKernelP_.add(4);
		Settings.pstepKernelAlpha_ = new ArrayList<Double>();
		Settings.pstepKernelAlpha_.add(2.0);

		// Run enrichment analysis
		EnrichMain enrichMain = new EnrichMain();
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
