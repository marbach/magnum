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
import edu.mit.magnum.Settings;


/**
 * Unit tests for NetworkAnalyzerBasicProperties
 */
public class MagnumOptionParserTest {
	

	// ============================================================================
	// SETUP
	
	@BeforeClass
	public static void testSetup() {
		Settings.loadSettings();
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
		
		assertEquals(22, Settings.mode_);
		assertEquals(23, Settings.randomSeed_);
		assertEquals(true, Settings.verbose_);
		assertEquals("myOutdir", Settings.outputDirectory_);

		assertEquals("myNetdir", Settings.networkDir_);
		assertEquals("myNet", Settings.networkFile_);
		assertEquals(true, Settings.isDirected_);
		assertEquals(true, Settings.isWeighted_);
		assertEquals(0.42, Settings.threshold_, 1e-12);
		
		assertEquals(true, Settings.computePstepKernel_);
		assertEquals(5, (int)Settings.pstepKernelP_.get(0));
		assertEquals(true, Settings.computeDegree_);
		assertEquals(true, Settings.computeBetweenness_);
		assertEquals(true, Settings.computeClusteringCoefficient_);
		assertEquals(true, Settings.computeShortestPathLengths_);
		assertEquals(true, Settings.computeUnion_);

		assertEquals("myGenes", Settings.geneCoordFile_);
		assertEquals("myScores", Settings.geneScoreFile_);
		assertEquals("myCmatrix", Settings.functionalDataFile_);
		assertEquals("myExcl", Settings.excludedGenesFile_);
		assertEquals(5.0, Settings.excludedGenesDistance_, 1e-12);
		assertEquals(6, Settings.numBins_);
		assertEquals(7, Settings.numPermutations_);
		assertEquals(0.8, Settings.curveCutoff_, 1e-12);
		
	}

	
	/** 
	 * Test default values of settings
	 * THIS IS NOT A TEST ON PURPOSE
	 * It's called by testParse(), because it has to be done before options are parsed 
	 */
	private void testDefaults() {
		
		double eps = 1e-12; 
		
		assertEquals(0, Settings.mode_);
		assertEquals(42, Settings.randomSeed_);
		assertEquals(false, Settings.verbose_);
		assertEquals(".", Settings.outputDirectory_);
		assertEquals("", Settings.outputFilename_);
		assertEquals(true, Settings.compressFiles_);

		assertEquals(".", Settings.networkDir_);
		assertEquals("", Settings.networkFile_);
		assertEquals("TAB", Settings.networkFileDelim_);
		assertEquals(false, Settings.isDirected_);
		assertEquals(false, Settings.removeSelfLoops_);
		assertEquals(false, Settings.isWeighted_);
		assertEquals(0.0, Settings.threshold_, eps);
		assertEquals(0.0, Settings.superHubThreshold_, eps);
		assertEquals("", Settings.refNodesFile_);

		assertEquals(false, Settings.computeUnion_);
		assertEquals("", Settings.networkGroupFile_);
		assertEquals("", Settings.networkFilePrefix_);
		assertEquals(false, Settings.computePairwiseSum_);
		assertEquals("", Settings.networkDir2_);

		assertEquals(false, Settings.computeDegree_);
		assertEquals(false, Settings.computeBetweenness_);
		assertEquals(false, Settings.computeClusteringCoefficient_);
		assertEquals(false, Settings.computeShortestPathLengths_);

		assertEquals(false, Settings.computePstepKernel_);
		assertEquals(2.0, Settings.pstepKernelAlpha_.get(0), eps);
		assertEquals(4, (int) Settings.pstepKernelP_.get(0));
		assertEquals(true, Settings.pstepKernelNormalize_);

		assertEquals(false, Settings.computeTargetTanimoto_);
		assertEquals(false, Settings.computeTfTanimoto_);

		assertEquals("", Settings.outputSuffix_);
		assertEquals(true, Settings.exportPairwiseNodeProperties_);
		assertEquals(true, Settings.exportNodeProperties_);

		assertEquals("", Settings.genesToBeLoadedFile_);

		assertEquals("", Settings.chromosome_);
		assertEquals(true, Settings.ignoreAllosomes_);

		assertEquals(true, Settings.loadOnlyProteinCodingGenes_);

		assertEquals("", Settings.geneCoordFile_);
		assertEquals("", Settings.geneCoordFile_);
		assertEquals(1e-6, Settings.genomeWideSignificanceThreshold_, eps);
		assertEquals(false, Settings.excludeGenomeWideSignificantGenes_);
		assertEquals("", Settings.functionalDataFile_);
		assertEquals(0, Settings.functionalDataCols_.size());

		assertEquals("", Settings.excludedGenesFile_);
		assertEquals("", Settings.excludedGenePairsFile_);
		assertEquals(1.0, Settings.excludedGenesDistance_, eps);

		assertEquals("custom", Settings.idTypeGeneScores_);
		assertEquals("custom", Settings.idTypeFunctionalData_);

		assertEquals(10000, Settings.numPermutations_);
		assertEquals(100, Settings.numBins_);
		assertEquals(false, Settings.scaleKernel_);

		assertEquals(10, Settings.constCurveResolution_);
		assertEquals(-1, Settings.varCurveResolution_);
		assertEquals(0.2, Settings.curveCutoff_, eps);
		assertEquals(-1, Settings.slidingWindowSize_);

		assertEquals(0.01, Settings.pval_.get(0), eps);
		assertEquals(0.05, Settings.pval_.get(1), eps);
		assertEquals(false, Settings.twoSidedTest_);
		assertEquals(false, Settings.controlFDR_);
		assertEquals(100, Settings.FDRStart_);
		assertEquals(10, Settings.AUCStart_);
		assertEquals(0, Settings.geneScoreIndexStart_);
		assertEquals(0, Settings.geneScoreIndexEnd_);

		assertEquals(0, Settings.numPermutationsExport_);
		assertEquals("", Settings.outputPrefix_);
	}
	
	


}
