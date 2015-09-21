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
package edu.mit.magnum.test;

import static org.junit.Assert.*;

import org.junit.*;

import edu.mit.magnum.Magnum;


/**
 * Unit tests for NetworkAnalyzerBasicProperties
 */
public class MagnumOptionParserTest {
	

	// ============================================================================
	// SETUP
	
	@BeforeClass
	public static void testSetup() {
		Magnum.set.resetToDefaults();
	}

	@AfterClass
	public static void testCleanup() {
	}
	  
	// ============================================================================
	// TESTS

	/** Test parsing of command-line arguments */
	@Test
	public void testParse() {

		// Check the default values
		testDefaults();
		
		String[] args = {"--mode", "22", 
				"--seed", "23", 
				"--verbose", 
				"--outdir", "myOutdir",
				"--netdir", "myNetdir",
				"--net", "myNet",
				"--directed",
				"--weighted",
				"--cutoff", "0.42",
				"--pstep",
				"--nsteps", "5",
				"--degree",
				"--betweenness",
				"--clustcoeff",
				"--shortestpath",
				"--union",
				"--genes", "myGenes",
				"--scores", "myScores",
				"--cmatrix", "myCmatrix",
				"--excl", "myExcl",
				"--neighbors", "5",
				"--bins", "6",
				"--permut", "7",
				"--curve", "0.8"
				}; 
		
		Magnum.set.parse(args);
		
		assertEquals(22, Magnum.set.mode_);
		assertEquals(23, Magnum.set.getRandomSeed());
		assertEquals(true, Magnum.set.verbose_);
		assertEquals("myOutdir", Magnum.set.outputDirectory_);

		assertEquals("myNetdir", Magnum.set.networkDir_);
		assertEquals("myNet", Magnum.set.networkFile_);
		assertEquals(true, Magnum.set.isDirected_);
		assertEquals(true, Magnum.set.isWeighted_);
		assertEquals(0.42, Magnum.set.threshold_, 1e-12);
		
		assertEquals(true, Magnum.set.computePstepKernel_);
		assertEquals(5, (int)Magnum.set.pstepKernelP_.get(0));
		assertEquals(true, Magnum.set.computeDegree_);
		assertEquals(true, Magnum.set.computeBetweenness_);
		assertEquals(true, Magnum.set.computeClusteringCoefficient_);
		assertEquals(true, Magnum.set.computeShortestPathLengths_);
		assertEquals(true, Magnum.set.computeUnion_);

		assertEquals("myGenes", Magnum.set.geneCoordFile_);
		assertEquals("myScores", Magnum.set.geneScoreFile_);
		assertEquals("myCmatrix", Magnum.set.functionalDataFile_);
		assertEquals("myExcl", Magnum.set.excludedGenesFile_);
		assertEquals(5.0, Magnum.set.excludedGenesDistance_, 1e-12);
		assertEquals(6, Magnum.set.numBins_);
		assertEquals(7, Magnum.set.numPermutations_);
		assertEquals(0.8, Magnum.set.curveCutoff_, 1e-12);
		
	}

	
	/** 
	 * Test default values of settings
	 * THIS IS NOT A TEST ON PURPOSE
	 * It's called by testParse(), because it has to be done before options are parsed 
	 */
	private void testDefaults() {
		
		double eps = 1e-12; 
		
		assertEquals(0, Magnum.set.mode_);
		assertEquals(42, Magnum.set.getRandomSeed());
		assertEquals(false, Magnum.set.verbose_);
		assertEquals(".", Magnum.set.outputDirectory_);
		assertEquals("", Magnum.set.outputFilename_);
		assertEquals(true, Magnum.set.compressFiles_);

		assertEquals(null, Magnum.set.networkDir_);
		assertEquals(null, Magnum.set.networkFile_);
		assertEquals("TAB", Magnum.set.networkFileDelim_);
		assertEquals(true, Magnum.set.isDirected_);
		assertEquals(true, Magnum.set.removeSelfLoops_);
		assertEquals(true, Magnum.set.isWeighted_);
		assertEquals(0.0, Magnum.set.threshold_, eps);
		assertEquals(0.0, Magnum.set.superHubThreshold_, eps);
		assertEquals(null, Magnum.set.refNodesFile_);

		assertEquals(false, Magnum.set.computeUnion_);
		assertEquals(null, Magnum.set.networkGroupFile_);
		assertEquals("", Magnum.set.networkFilePrefix_);
		assertEquals(false, Magnum.set.computePairwiseSum_);
		assertEquals(null, Magnum.set.networkDir2_);

		assertEquals(false, Magnum.set.computeDegree_);
		assertEquals(false, Magnum.set.computeBetweenness_);
		assertEquals(false, Magnum.set.computeClusteringCoefficient_);
		assertEquals(false, Magnum.set.computeShortestPathLengths_);

		assertEquals(false, Magnum.set.computePstepKernel_);
		assertEquals(2.0, Magnum.set.pstepKernelAlpha_, eps);
		assertEquals(4, (int) Magnum.set.pstepKernelP_.get(0));
		assertEquals(true, Magnum.set.pstepKernelNormalize_);

		assertEquals(false, Magnum.set.computeTargetTanimoto_);
		assertEquals(false, Magnum.set.computeTfTanimoto_);

		assertEquals("", Magnum.set.outputSuffix_);
		assertEquals(true, Magnum.set.exportPairwiseNodeProperties_);
		assertEquals(true, Magnum.set.exportNodeProperties_);

		assertEquals(null, Magnum.set.genesToBeLoadedFile_);

		assertEquals(null, Magnum.set.chromosome_);
		assertEquals(true, Magnum.set.ignoreAllosomes_);

		assertEquals(true, Magnum.set.loadOnlyProteinCodingGenes_);

		assertEquals(null, Magnum.set.geneCoordFile_);
		assertEquals(1e-6, Magnum.set.genomeWideSignificanceThreshold_, eps);
		assertEquals(false, Magnum.set.excludeGenomeWideSignificantGenes_);
		assertEquals(null, Magnum.set.functionalDataFile_);
		//assertEquals(0, Magnum.set.functionalDataCols_.size());
		assertEquals(null, Magnum.set.functionalDataCols_);

		assertEquals(null, Magnum.set.excludedGenesFile_);
		assertEquals(null, Magnum.set.excludedGenePairsFile_);
		assertEquals(1.0, Magnum.set.excludedGenesDistance_, eps);

		assertEquals("custom", Magnum.set.idTypeGeneScores_);
		assertEquals("custom", Magnum.set.idTypeFunctionalData_);

		assertEquals(10000, Magnum.set.numPermutations_);
		assertEquals(100, Magnum.set.numBins_);
		assertEquals(false, Magnum.set.scaleKernel_);

		assertEquals(10, Magnum.set.constCurveResolution_);
		assertEquals(-1, Magnum.set.varCurveResolution_);
		assertEquals(0.2, Magnum.set.curveCutoff_, eps);
		assertEquals(-1, Magnum.set.slidingWindowSize_);

		assertEquals(0.01, Magnum.set.pval_.get(0), eps);
		assertEquals(0.05, Magnum.set.pval_.get(1), eps);
		assertEquals(false, Magnum.set.twoSidedTest_);
		assertEquals(false, Magnum.set.controlFDR_);
		assertEquals(100, Magnum.set.FDRStart_);
		assertEquals(10, Magnum.set.AUCStart_);
		assertEquals(0, Magnum.set.geneScoreIndexStart_);
		assertEquals(0, Magnum.set.geneScoreIndexEnd_);

		assertEquals(0, Magnum.set.numPermutationsExport_);
	}
	
	


}
