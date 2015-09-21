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

import edu.mit.magnum.MagnumOptionParser;
import edu.mit.magnum.MagnumSettings;


/**
 * Unit tests for NetworkAnalyzerBasicProperties
 */
public class MagnumOptionParserTest {
	

	// ============================================================================
	// SETUP
	
	@BeforeClass
	public static void testSetup() {
		MagnumSettings.loadSettings();
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
		
		MagnumOptionParser parser =  new MagnumOptionParser();
		parser.parse(args);
		
		assertEquals(22, MagnumSettings.mode_);
		assertEquals(23, MagnumSettings.randomSeed_);
		assertEquals(true, MagnumSettings.verbose_);
		assertEquals("myOutdir", MagnumSettings.outputDirectory_);

		assertEquals("myNetdir", MagnumSettings.networkDir_);
		assertEquals("myNet", MagnumSettings.networkFile_);
		assertEquals(true, MagnumSettings.isDirected_);
		assertEquals(true, MagnumSettings.isWeighted_);
		assertEquals(0.42, MagnumSettings.threshold_, 1e-12);
		
		assertEquals(true, MagnumSettings.computePstepKernel_);
		assertEquals(5, (int)MagnumSettings.pstepKernelP_.get(0));
		assertEquals(true, MagnumSettings.computeDegree_);
		assertEquals(true, MagnumSettings.computeBetweenness_);
		assertEquals(true, MagnumSettings.computeClusteringCoefficient_);
		assertEquals(true, MagnumSettings.computeShortestPathLengths_);
		assertEquals(true, MagnumSettings.computeUnion_);

		assertEquals("myGenes", MagnumSettings.geneCoordFile_);
		assertEquals("myScores", MagnumSettings.geneScoreFile_);
		assertEquals("myCmatrix", MagnumSettings.functionalDataFile_);
		assertEquals("myExcl", MagnumSettings.excludedGenesFile_);
		assertEquals(5.0, MagnumSettings.excludedGenesDistance_, 1e-12);
		assertEquals(6, MagnumSettings.numBins_);
		assertEquals(7, MagnumSettings.numPermutations_);
		assertEquals(0.8, MagnumSettings.curveCutoff_, 1e-12);
		
	}

	
	/** 
	 * Test default values of settings
	 * THIS IS NOT A TEST ON PURPOSE
	 * It's called by testParse(), because it has to be done before options are parsed 
	 */
	private void testDefaults() {
		
		double eps = 1e-12; 
		
		assertEquals(0, MagnumSettings.mode_);
		assertEquals(42, MagnumSettings.randomSeed_);
		assertEquals(false, MagnumSettings.verbose_);
		assertEquals(".", MagnumSettings.outputDirectory_);
		assertEquals("", MagnumSettings.outputFilename_);
		assertEquals(true, MagnumSettings.compressFiles_);

		assertEquals(".", MagnumSettings.networkDir_);
		assertEquals("", MagnumSettings.networkFile_);
		assertEquals("TAB", MagnumSettings.networkFileDelim_);
		assertEquals(false, MagnumSettings.isDirected_);
		assertEquals(false, MagnumSettings.removeSelfLoops_);
		assertEquals(false, MagnumSettings.isWeighted_);
		assertEquals(0.0, MagnumSettings.threshold_, eps);
		assertEquals(0.0, MagnumSettings.superHubThreshold_, eps);
		assertEquals("", MagnumSettings.refNodesFile_);

		assertEquals(false, MagnumSettings.computeUnion_);
		assertEquals("", MagnumSettings.networkGroupFile_);
		assertEquals("", MagnumSettings.networkFilePrefix_);
		assertEquals(false, MagnumSettings.computePairwiseSum_);
		assertEquals("", MagnumSettings.networkDir2_);

		assertEquals(false, MagnumSettings.computeDegree_);
		assertEquals(false, MagnumSettings.computeBetweenness_);
		assertEquals(false, MagnumSettings.computeClusteringCoefficient_);
		assertEquals(false, MagnumSettings.computeShortestPathLengths_);

		assertEquals(false, MagnumSettings.computePstepKernel_);
		assertEquals(2.0, MagnumSettings.pstepKernelAlpha_.get(0), eps);
		assertEquals(4, (int) MagnumSettings.pstepKernelP_.get(0));
		assertEquals(true, MagnumSettings.pstepKernelNormalize_);

		assertEquals(false, MagnumSettings.computeTargetTanimoto_);
		assertEquals(false, MagnumSettings.computeTfTanimoto_);

		assertEquals("", MagnumSettings.outputSuffix_);
		assertEquals(true, MagnumSettings.exportPairwiseNodeProperties_);
		assertEquals(true, MagnumSettings.exportNodeProperties_);

		assertEquals("", MagnumSettings.genesToBeLoadedFile_);

		assertEquals("", MagnumSettings.chromosome_);
		assertEquals(true, MagnumSettings.ignoreAllosomes_);

		assertEquals(true, MagnumSettings.loadOnlyProteinCodingGenes_);

		assertEquals("", MagnumSettings.geneCoordFile_);
		assertEquals("", MagnumSettings.geneCoordFile_);
		assertEquals(1e-6, MagnumSettings.genomeWideSignificanceThreshold_, eps);
		assertEquals(false, MagnumSettings.excludeGenomeWideSignificantGenes_);
		assertEquals("", MagnumSettings.functionalDataFile_);
		assertEquals(0, MagnumSettings.functionalDataCols_.size());

		assertEquals("", MagnumSettings.excludedGenesFile_);
		assertEquals("", MagnumSettings.excludedGenePairsFile_);
		assertEquals(1.0, MagnumSettings.excludedGenesDistance_, eps);

		assertEquals("custom", MagnumSettings.idTypeGeneScores_);
		assertEquals("custom", MagnumSettings.idTypeFunctionalData_);

		assertEquals(10000, MagnumSettings.numPermutations_);
		assertEquals(100, MagnumSettings.numBins_);
		assertEquals(false, MagnumSettings.scaleKernel_);

		assertEquals(10, MagnumSettings.constCurveResolution_);
		assertEquals(-1, MagnumSettings.varCurveResolution_);
		assertEquals(0.2, MagnumSettings.curveCutoff_, eps);
		assertEquals(-1, MagnumSettings.slidingWindowSize_);

		assertEquals(0.01, MagnumSettings.pval_.get(0), eps);
		assertEquals(0.05, MagnumSettings.pval_.get(1), eps);
		assertEquals(false, MagnumSettings.twoSidedTest_);
		assertEquals(false, MagnumSettings.controlFDR_);
		assertEquals(100, MagnumSettings.FDRStart_);
		assertEquals(10, MagnumSettings.AUCStart_);
		assertEquals(0, MagnumSettings.geneScoreIndexStart_);
		assertEquals(0, MagnumSettings.geneScoreIndexEnd_);

		assertEquals(0, MagnumSettings.numPermutationsExport_);
		assertEquals("", MagnumSettings.outputPrefix_);
	}
	
	


}
