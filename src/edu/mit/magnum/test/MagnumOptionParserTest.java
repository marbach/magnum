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

import java.io.File;

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

		String[] args = {"--mode", "22", 
				"--seed", "23", 
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
		assertEquals(new File("myOutdir"), Magnum.set.outputDirectory_);

		assertEquals(new File("myNetdir"), Magnum.set.networkDir_);
		assertEquals(new File("myNet"), Magnum.set.networkFile_);
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

		assertEquals(new File("myGenes"), Magnum.set.geneCoordFile_);
		assertEquals(new File("myScores"), Magnum.set.geneScoreFile_);
		assertEquals(new File("myCmatrix"), Magnum.set.functionalDataFile_);
		assertEquals(new File("myExcl"), Magnum.set.excludedGenesFile_);
		assertEquals(5.0, Magnum.set.excludedGenesDistance_, 1e-12);
		assertEquals(6, Magnum.set.numBins_);
		assertEquals(7, Magnum.set.numPermutations_);
		assertEquals(0.8, Magnum.set.curveCutoff_, 1e-12);
		
	}

}
