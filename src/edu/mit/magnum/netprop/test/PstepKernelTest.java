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
package edu.mit.magnum.netprop.test;

import static org.junit.Assert.*;
import java.util.ArrayList;
import org.junit.*;

import cern.colt.matrix.DoubleMatrix2D;
import edu.mit.magnum.Settings;
import edu.mit.magnum.net.*;
import edu.mit.magnum.netprop.*;


/**
 * Unit tests for NetworkAnalyzerShortestPaths
 */
public class PstepKernelTest {
	
	
	// ============================================================================
	// SETUP
	
	@BeforeClass
	public static void testSetup() {
		Settings.loadSettings();
		Settings.superHubThreshold_ = 0;
		Settings.computePstepKernel_ = true;
		Settings.exportNodeProperties_ = true;
		Settings.pstepKernelNormalize_ = true;
	}

	@AfterClass
	public static void testCleanup() {
	}
	  
	// ============================================================================
	// TESTS

	/** Pstep kernel computation */
	@Test
	public void testComputeK() {

		// Load undirected network without self-loops
		Network testNet = new Network("src/edu/mit/magnum/netprop/test/simpleNet.txt", false, false);
		
		ArrayList<Integer> numSteps = new ArrayList<Integer>();
		numSteps.add(1);
		numSteps.add(3);
		numSteps.add(4); 
		ArrayList<Double> alpha = new ArrayList<Double>();
		alpha.add(2.0);
		PstepKernel test = new PstepKernel(testNet, alpha, numSteps, true, false);
		
		// Compute pstep kernel and corresponding node centrality
		test.run();
		DoubleMatrix2D K = test.getK();
		//Double[] centrality = test.getCentrality();
		
		// Expected distances
		double[][] K_expected = {
				{ 0.716723549, 0.81088013, 0.4006134, 0.09458298, 0.009653335, 0.009653335 },
				{ 0.810880131, 1.00000000, 0.6279863, 0.24244096, 0.061433447, 0.061433447 },
				{ 0.400613398, 0.62798635, 0.6313993, 0.47930857, 0.245733788, 0.245733788 },
				{ 0.094582979, 0.24244096, 0.4793086, 0.82593857, 0.705029471, 0.705029471 },
				{ 0.009653335, 0.06143345, 0.2457338, 0.70502947, 0.726962457, 0.716723549 },
				{ 0.009653335, 0.06143345, 0.2457338, 0.70502947, 0.716723549, 0.726962457 }};
		
		// Check distances
		double epsilon = 1e-6;
		for (int i=0; i<6; i++)
			for (int j=0; j<6; j++)
				assertEquals(K_expected[i][j], K.get(i, j), epsilon);
	}

	


	// ============================================================================
	// PRIVATE METHODS

}
