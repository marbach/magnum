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
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import edu.mit.magnum.Settings;
import edu.mit.magnum.net.*;
import edu.mit.magnum.netprop.*;


/**
 * Unit tests for NetworkAnalyzerShortestPaths
 */
public class TanimotoCoefficientTest {
	
	
	// ============================================================================
	// SETUP
	
	@BeforeClass
	public static void testSetup() {
		Settings.loadSettings();
		Settings.superHubThreshold_ = 0;
	}

	@AfterClass
	public static void testCleanup() {
	}
	  
	// ============================================================================
	// TESTS

	/** Tanimoto computation for targets */
	@Test
	public void testTargetTanimoto() {

		// Load directed, weighted network with self-loops
		Network testNet = new Network("src/edu/mit/magnum/netprop/test/tanimotoTestNet.txt", true, false, true, 0);
		// Initialize
		TanimotoCoefficient tanimoto = new TanimotoCoefficient(testNet, true, false);
		// Compute tanimoto similarity matrix and corresponding node centrality
		tanimoto.run();

		DoubleMatrix2D K = tanimoto.getK();
		ArrayList<Node> nodes = tanimoto.getNodes();
		
		assertEquals(7, nodes.size());
		//Double[] centrality = test.getCentrality();
		
		// Expected Tanimoto coefficients
		SparseDoubleMatrix2D K_expected = new SparseDoubleMatrix2D(7, 7);
		K_expected.set(1, 2, 0.9285714);
		K_expected.set(2, 1, 0.9285714);
		
		K_expected.set(4, 5, 0.9863014);
		K_expected.set(5, 4, 0.9863014);

		K_expected.set(4, 6, 0.1584158);
		K_expected.set(6, 4, 0.1584158);

		K_expected.set(5, 6, 0.1551724);
		K_expected.set(6, 5, 0.1551724);

		K_expected.set(3, 6, 0.893617);
		K_expected.set(6, 3, 0.893617);

		for (int i=0; i<7; i++)
			K_expected.set(i, i, 1);

		// Check distances
		double epsilon = 1e-6;
		for (int i=0; i<7; i++)
			for (int j=0; j<7; j++)
				assertEquals(K_expected.get(i, j), K.get(i, j), epsilon);
	}

	
	/** Tanimoto computation for targets */
	@Test
	public void testRegulatorTanimoto() {

		// Load directed, weighted network with self-loops
		Network testNet = new Network("src/edu/mit/magnum/netprop/test/tanimotoTestNet.txt", true, false, true, 0);
		// Initialize
		TanimotoCoefficient tanimoto = new TanimotoCoefficient(testNet, false, false);
		// Compute tanimoto similarity matrix and corresponding node centrality
		tanimoto.run();

		DoubleMatrix2D K = tanimoto.getK();
		ArrayList<Node> nodes = tanimoto.getNodes();
		
		assertEquals(5, nodes.size());
		//Double[] centrality = test.getCentrality();
		
		// Expected Tanimoto coefficients
		SparseDoubleMatrix2D K_expected = new SparseDoubleMatrix2D(5, 5);
		K_expected.set(1, 2, 0.7419355);
		K_expected.set(2, 1, 0.7419355);
		
		K_expected.set(3, 4, 0.06363636);
		K_expected.set(4, 3, 0.06363636);

		for (int i=0; i<5; i++)
			K_expected.set(i, i, 1);

		// Check distances
		double epsilon = 1e-6;
		for (int i=0; i<5; i++)
			for (int j=0; j<5; j++)
				assertEquals(K_expected.get(i, j), K.get(i, j), epsilon);
	}
	


	// ============================================================================
	// PRIVATE METHODS

}
