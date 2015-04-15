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
public class ShortestPathsTest {
	
	
	// ============================================================================
	// SETUP
	
	@BeforeClass
	public static void testSetup() {
		Settings.loadSettings();
		Settings.superHubThreshold_ = 0;
		Settings.computeShortestPathLengths_ = true;
		Settings.exportNodeProperties_ = true;
	}

	@AfterClass
	public static void testCleanup() {
	}
	  
	// ============================================================================
	// TESTS

	/** Test undirected distances and closeness */
	@Test
	public void testUndirected() {

		// Load undirected network without self-loops
		Network testNet = new Network("src/edu/mit/magnum/netprop/test/degreeTestNet.txt", false, false);
		
		// Test that prior nodes work by setting them in arbitrary order
		ArrayList<Node> priorNodes = new ArrayList<Node>();
		priorNodes.add(testNet.getNode(6));
		priorNodes.add(testNet.getNode(5));
		priorNodes.add(testNet.getNode(1));
		priorNodes.add(testNet.getNode(2));
		priorNodes.add(testNet.getNode(0));
		priorNodes.add(testNet.getNode(4));
		priorNodes.add(testNet.getNode(3));
		testNet.setRefNodes(priorNodes);
		
		ShortestPaths test = new ShortestPaths(testNet, true);
		
		// Compute distance matrix and closeness centrality
		test.run();
		DoubleMatrix2D distances = test.getK();
		Double[] closeness = test.getCentrality();
		
		// Node indexes
		int[] node = new int[7];
		for (int i=0; i<7; i++)
			node[i] = testNet.getNodeIndex(Integer.toString(i+1));

		// Prior node indexes
		int[] priorNode = new int[7];
		for (int i=0; i<7; i++)
			priorNode[i] = testNet.getRefNodeIndex(new Node(Integer.toString(i+1)));

		// Expected distances
		int[][] expectedDistance = {
				{ 0, 1, 1, 1, 1, -1, -1 },
				{ 1, 0, 1, 2, 2, -1, -1 },
				{ 1, 1, 0, 1, 2, -1, -1 },
				{ 1, 2, 1, 0, 2, -1, -1 },
				{ 1, 2, 2, 2, 0, -1, -1 }, 
				{ -1, -1, -1, -1, -1, 0, 1 },
				{ -1, -1, -1, -1, -1, 1, 0 }};
		
		// Check distances
		double epsilon = 1e-6;
		for (int i=0; i<7; i++)
			for (int j=0; j<7; j++)
				assertEquals(expectedDistance[i][j], distances.get(node[i], priorNode[j]), epsilon);
		
		// Check closeness
		assertEquals(closeness[node[0]], 1/(4/4.0), epsilon);
		assertEquals(closeness[node[1]], 1/(6/4.0), epsilon);
		assertEquals(closeness[node[2]], 1/(5/4.0), epsilon);
		assertEquals(closeness[node[3]], 1/(6/4.0), epsilon);
		assertEquals(closeness[node[4]], 1/(7/4.0), epsilon);
		assertEquals(closeness[node[5]], 1/(1/1.0), epsilon);
		assertEquals(closeness[node[6]], 1/(1/1.0), epsilon);
	}

	
	// ----------------------------------------------------------------------------
	
	/** Test undirected distances and closeness with priors */
	@Test
	public void testUndirectedWithPriors() {

		// Load undirected network without self-loops
		Network testNet = new Network("src/edu/mit/magnum/netprop/test/degreeTestNet.txt", false, false);
		
		// Test that prior nodes work by setting them in arbitrary order
		ArrayList<Node> priors = new ArrayList<Node>();
		priors.add(testNet.getNode("4"));
		priors.add(testNet.getNode("2"));
		priors.add(testNet.getNode("7"));
		testNet.setRefNodes(priors);
		
		ShortestPaths test = new ShortestPaths(testNet, true);
		
		// Compute distance matrix and closeness centrality
		test.run();
		DoubleMatrix2D distances = test.getK();
		Double[] closeness = test.getCentrality();
		
		// Node indexes
		int[] node = new int[7];
		for (int i=0; i<7; i++)
			node[i] = testNet.getNodeIndex(Integer.toString(i+1));

		// Prior node indexes
		int[] priorNode = new int[3];
		priorNode[0] = testNet.getRefNodeIndex(priors.get(0));
		priorNode[1] = testNet.getRefNodeIndex(priors.get(1));
		priorNode[2] = testNet.getRefNodeIndex(priors.get(2));

		// Expected distances
		int[][] expectedDistance = {
				{ 1, 1, -1 },
				{ 2, 0, -1 },
				{ 1, 1, -1 },
				{ 0, 2, -1 },
				{ 2, 2, -1 }, 
				{ -1, -1, 1 },
				{ -1, -1, 0 }};
		
		// Check distances
		double epsilon = 1e-6;
		for (int i=0; i<7; i++)
			for (int j=0; j<3; j++)
				assertEquals(expectedDistance[i][j], distances.get(node[i], priorNode[j]), epsilon);
		
		// Check closeness
		assertEquals(closeness[node[0]], 1/(2/2.0), epsilon);
		assertEquals(closeness[node[1]], 1/(2/2.0), epsilon);
		assertEquals(closeness[node[2]], 1/(2/2.0), epsilon);
		assertEquals(closeness[node[3]], 1/(2/2.0), epsilon);
		assertEquals(closeness[node[4]], 1/(4/2.0), epsilon);
		assertEquals(closeness[node[5]], 1/(1/1), epsilon);
		assertEquals(closeness[node[6]], 1e12, epsilon);
	}


	// ----------------------------------------------------------------------------

	/** Test directed distances and closeness */
	@Test
	public void testDirected() {

		// Load undirected network without self-loops
		Network testNet = new Network("src/edu/mit/magnum/netprop/test/degreeTestNet.txt", true, false);
		
		// Test that prior nodes work by setting them in arbitrary order
		ArrayList<Node> priorNodes = new ArrayList<Node>();
		priorNodes.add(testNet.getNode(6));
		priorNodes.add(testNet.getNode(5));
		priorNodes.add(testNet.getNode(1));
		priorNodes.add(testNet.getNode(2));
		priorNodes.add(testNet.getNode(0));
		priorNodes.add(testNet.getNode(4));
		priorNodes.add(testNet.getNode(3));
		testNet.setRefNodes(priorNodes);
		
		ShortestPaths test = new ShortestPaths(testNet, true);
		
		// Compute distance matrix and closeness centrality
		test.run();
		DoubleMatrix2D distances = test.getK();
		Double[] closeness = test.getCentrality();
		
		// Node indexes
		int[] node = new int[7];
		for (int i=0; i<7; i++)
			node[i] = testNet.getNodeIndex(Integer.toString(i+1));

		// Prior node indexes
		int[] priorNode = new int[7];
		for (int i=0; i<7; i++)
			priorNode[i] = testNet.getRefNodeIndex(new Node(Integer.toString(i+1)));

		// Expected distances
		int[][] expectedDistance = {
				{ 0, 1, 1, 1, 1, -1, -1 },
				{ -1, 0, 1, -1, -1, -1, -1 },
				{ -1, -1, 0, -1, -1, -1, -1 },
				{ 1, 2, 1, 0, 2, -1, -1 },
				{ -1, -1, -1, -1, 0, -1, -1 }, 
				{ -1, -1, -1, -1, -1, 0, 1 },
				{ -1, -1, -1, -1, -1, -1, 0 }};
		
		// Check distances
		double epsilon = 1e-6;
		for (int i=0; i<7; i++)
			for (int j=0; j<7; j++)
				assertEquals(expectedDistance[i][j], distances.get(node[i], priorNode[j]), epsilon);
		
		// Check closeness
		assertEquals(closeness[node[0]], 1/(4/4.0), epsilon);
		assertEquals(closeness[node[1]], 1/(1/1.0), epsilon);
		assertEquals(closeness[node[2]], 1e12, epsilon);
		assertEquals(closeness[node[3]], 1/(6/4.0), epsilon);
		assertEquals(closeness[node[4]], 1e12, epsilon);
		assertEquals(closeness[node[5]], 1/(1/1.0), epsilon);
		assertEquals(closeness[node[6]], 1e12, epsilon);
	}

	
	// ----------------------------------------------------------------------------
	
	/** Test undirected distances and closeness with priors */
	@Test
	public void testDirectedWithPriors() {

		// Load undirected network without self-loops
		Network testNet = new Network("src/edu/mit/magnum/netprop/test/degreeTestNet.txt", true, false);
		
		// Test that prior nodes work by setting them in arbitrary order
		ArrayList<Node> priors = new ArrayList<Node>();
		priors.add(testNet.getNode("4"));
		priors.add(testNet.getNode("2"));
		priors.add(testNet.getNode("7"));
		testNet.setRefNodes(priors);
		
		ShortestPaths test = new ShortestPaths(testNet, true);
		
		// Compute distance matrix and closeness centrality
		test.run();
		DoubleMatrix2D distances = test.getK();
		Double[] closeness = test.getCentrality();
		
		// Node indexes
		int[] node = new int[7];
		for (int i=0; i<7; i++)
			node[i] = testNet.getNodeIndex(Integer.toString(i+1));

		// Prior node indexes
		int[] priorNode = new int[3];
		priorNode[0] = testNet.getRefNodeIndex(priors.get(0));
		priorNode[1] = testNet.getRefNodeIndex(priors.get(1));
		priorNode[2] = testNet.getRefNodeIndex(priors.get(2));

		// Expected distances
		int[][] expectedDistance = {
				{ 1, 1, -1 },
				{ -1, 0, -1 },
				{ -1, -1, -1 },
				{ 0, 2, -1 },
				{ -1, -1, -1 }, 
				{ -1, -1, 1 },
				{ -1, -1, 0 }};
		
		// Check distances
		double epsilon = 1e-6;
		for (int i=0; i<7; i++)
			for (int j=0; j<3; j++)
				assertEquals(expectedDistance[i][j], distances.get(node[i], priorNode[j]), epsilon);
		
		// Check closeness
		assertEquals(closeness[node[0]], 1/(2/2.0), epsilon);
		assertEquals(closeness[node[1]], 1e12, epsilon);
		assertEquals(closeness[node[2]], 0, epsilon);
		assertEquals(closeness[node[3]], 1/(2/2.0), epsilon);
		assertEquals(closeness[node[4]], 0, epsilon);
		assertEquals(closeness[node[5]], 1/(1/1), epsilon);
		assertEquals(closeness[node[6]], 1e12, epsilon);
	}



	// ============================================================================
	// PRIVATE METHODS

}
