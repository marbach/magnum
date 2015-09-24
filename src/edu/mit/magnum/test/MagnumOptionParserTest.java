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
	
	/** The magnum instance */
	private static Magnum mag = new Magnum();

	// ============================================================================
	// SETUP
	
	@BeforeClass
	public static void testSetup() {
		mag.set.resetToDefaults();
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
		
		mag.set.parse(args);
		
		assertEquals(22, mag.set.mode_);
		assertEquals(23, mag.set.getRandomSeed());
		assertEquals(new File("myOutdir"), mag.set.outputDirectory_);

		assertEquals(new File("myNetdir"), mag.set.networkDir_);
		assertEquals(new File("myNet"), mag.set.networkFile_);
		assertEquals(true, mag.set.isDirected_);
		assertEquals(true, mag.set.isWeighted_);
		assertEquals(0.42, mag.set.threshold_, 1e-12);
		
		assertEquals(true, mag.set.computePstepKernel_);
		assertEquals(5, (int)mag.set.pstepKernelP_.get(0));
		assertEquals(true, mag.set.computeDegree_);
		assertEquals(true, mag.set.computeBetweenness_);
		assertEquals(true, mag.set.computeClusteringCoefficient_);
		assertEquals(true, mag.set.computeShortestPathLengths_);
		assertEquals(true, mag.set.computeUnion_);

		assertEquals(new File("myGenes"), mag.set.geneCoordFile_);
		assertEquals(new File("myScores"), mag.set.geneScoreFile_);
		assertEquals(new File("myCmatrix"), mag.set.functionalDataFile_);
		assertEquals(new File("myExcl"), mag.set.excludedGenesFile_);
		assertEquals(5.0, mag.set.excludedGenesDistance_, 1e-12);
		assertEquals(6, mag.set.numBins_);
		assertEquals(7, mag.set.numPermutations_);
		assertEquals(0.8, mag.set.curveCutoff_, 1e-12);
		
	}

}
