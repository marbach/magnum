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

import org.junit.*;

import edu.mit.magnum.MagnumSettings;
import edu.mit.magnum.net.*;
import edu.mit.magnum.netprop.*;


/**
 * Unit tests for NetworkAnalyzerBasicProperties
 */
public class BasicPropertiesTest {
	
	/** Test network */
	private String ppiNetFile_ = "src/edu/mit/magnum/netprop/test/ppi.txt";

	// ============================================================================
	// SETUP
	
	@BeforeClass
	public static void testSetup() {
		MagnumSettings.loadSettings();
		MagnumSettings.superHubThreshold_ = 0;
	}

	@AfterClass
	public static void testCleanup() {
	}
	  
	// ============================================================================
	// TESTS

	/** Node degrees */
	@Test
	public void testDegree() {

		// Load undirected network without self-loops
		Network testNet = new Network("src/edu/mit/magnum/netprop/test/degreeTestNet.txt", false, true);
		BasicProperties test = new BasicProperties(testNet);

		// Check number of nodes and edges
		assertEquals(testNet.getNumNodes(), 7);
		assertEquals(testNet.getNumEdges(), 7);
		
		// Compute degrees
		test.computeDegree();
		Integer[] degree = test.getDegree();
		
		// Check degrees
		assertEquals((int)degree[testNet.getNodeIndex("1")], 4);
		assertEquals((int)degree[testNet.getNodeIndex("2")], 2);
		assertEquals((int)degree[testNet.getNodeIndex("3")], 3);
		assertEquals((int)degree[testNet.getNodeIndex("4")], 2);
		assertEquals((int)degree[testNet.getNodeIndex("5")], 1);
		assertEquals((int)degree[testNet.getNodeIndex("6")], 1);
		assertEquals((int)degree[testNet.getNodeIndex("7")], 1);
		
		// Load the ppi and test that the degrees sum up to the number of edges
		testNet = new Network(ppiNetFile_, false, true);
		test = new BasicProperties(testNet);
		
		assertEquals(testNet.getNumEdges(), 96967);

		// Compute node degrees 
		test.computeDegree();
		degree = test.getDegree();
		
		// Sum node degrees
		int sum = 0;
		for (int i=0; i<degree.length; i++)
			sum += degree[i];
		
		assertEquals(sum, 2*testNet.getNumEdges());
	}

	
	// ----------------------------------------------------------------------------

	/** Degrees for directed networks */
	@Test
	public void testDegreesDirected() {
		
		// Load directed network with self-loops
		Network testNet = new Network("src/edu/mit/magnum/netprop/test/degreeTestNet.txt", true, false);
		BasicProperties test = new BasicProperties(testNet);

		// Check number of nodes and edges
		assertEquals(testNet.getNumNodes(), 7);
		assertEquals(testNet.getNumEdges(), 9);
		
		// Compute degrees
		test.computeDegree();
		Integer[] degree = test.getDegree();
		Integer[] outdegree = test.getOutdegree();
		Integer[] indegree = test.getIndegree();
		
		// Check degrees
		assertEquals((int)degree[testNet.getNodeIndex("1")], 5);
		assertEquals((int)degree[testNet.getNodeIndex("2")], 2);
		assertEquals((int)degree[testNet.getNodeIndex("3")], 3);
		assertEquals((int)degree[testNet.getNodeIndex("4")], 3);
		assertEquals((int)degree[testNet.getNodeIndex("5")], 1);
		assertEquals((int)degree[testNet.getNodeIndex("6")], 1);
		// self-loops are only counted once
		assertEquals((int)degree[testNet.getNodeIndex("7")], 2);

		assertEquals((int)outdegree[testNet.getNodeIndex("1")], 4);
		assertEquals((int)outdegree[testNet.getNodeIndex("2")], 1);
		assertEquals((int)outdegree[testNet.getNodeIndex("3")], 0);
		assertEquals((int)outdegree[testNet.getNodeIndex("4")], 2);
		assertEquals((int)outdegree[testNet.getNodeIndex("5")], 0);
		assertEquals((int)outdegree[testNet.getNodeIndex("6")], 1);
		assertEquals((int)outdegree[testNet.getNodeIndex("7")], 1);

		assertEquals((int)indegree[testNet.getNodeIndex("1")], 1);
		assertEquals((int)indegree[testNet.getNodeIndex("2")], 1);
		assertEquals((int)indegree[testNet.getNodeIndex("3")], 3);
		assertEquals((int)indegree[testNet.getNodeIndex("4")], 1);
		assertEquals((int)indegree[testNet.getNodeIndex("5")], 1);
		assertEquals((int)indegree[testNet.getNodeIndex("6")], 0);
		assertEquals((int)indegree[testNet.getNodeIndex("7")], 2);
	}

	
	// ----------------------------------------------------------------------------

	/** Betweenness centrality */
	@Test
	public void testBetweenness() {
		
		// Load undirected network without self-loops
		Network testNet = new Network("src/edu/mit/magnum/netprop/test/degreeTestNet.txt", false, false);
		BasicProperties test = new BasicProperties(testNet);

		// Compute betweenness
		test.computeBetweenness();
		Double[] betweenness = test.getBetweenness();
		
		double epsilon = 1e-6;
		assertEquals(betweenness[testNet.getNodeIndex("1")], 3.5, epsilon);
		assertEquals(betweenness[testNet.getNodeIndex("2")], 0, epsilon);
		assertEquals(betweenness[testNet.getNodeIndex("3")], 0.5, epsilon);
		assertEquals(betweenness[testNet.getNodeIndex("4")], 0, epsilon);
		assertEquals(betweenness[testNet.getNodeIndex("5")], 0, epsilon);
		assertEquals(betweenness[testNet.getNodeIndex("6")], 0, epsilon);
		assertEquals(betweenness[testNet.getNodeIndex("7")], 0, epsilon);
	}


	// ----------------------------------------------------------------------------

	/** Betweenness centrality */
	@Test
	public void testBetweennessDirected() {
		
		// Load directed network with self-loops
		Network testNet = new Network("src/edu/mit/magnum/netprop/test/degreeTestNet.txt", true, false);
		BasicProperties test = new BasicProperties(testNet);

		// Compute betweenness
		test.computeBetweenness();
		Double[] betweenness = test.getBetweenness();
		
		double epsilon = 1e-6;
		assertEquals(betweenness[testNet.getNodeIndex("1")], 2, epsilon);
		assertEquals(betweenness[testNet.getNodeIndex("2")], 0, epsilon);
		assertEquals(betweenness[testNet.getNodeIndex("3")], 0, epsilon);
		assertEquals(betweenness[testNet.getNodeIndex("4")], 0, epsilon);
		assertEquals(betweenness[testNet.getNodeIndex("5")], 0, epsilon);
		assertEquals(betweenness[testNet.getNodeIndex("6")], 0, epsilon);
		assertEquals(betweenness[testNet.getNodeIndex("7")], 0, epsilon);
	}
	
	
	// ----------------------------------------------------------------------------

	/** Clustering coefficient */
	@Test
	public void testClusteringCoeff() {
		
		// Load undirected network with self-loops
		Network testNet = new Network("src/edu/mit/magnum/netprop/test/degreeTestNet.txt", false, false);
		BasicProperties test = new BasicProperties(testNet);

		// Compute betweenness
		test.computeClusteringCoefficient();
		Double[] clust = test.getClusteringCoeff();
		
		double epsilon = 1e-6;
		assertEquals(clust[testNet.getNodeIndex("1")], clusteringCoeffUndir(2*2, 4), epsilon);
		assertEquals(clust[testNet.getNodeIndex("2")], clusteringCoeffUndir(2*1, 2), epsilon);
		assertEquals(clust[testNet.getNodeIndex("3")], clusteringCoeffUndir(2*2, 3), epsilon);
		assertEquals(clust[testNet.getNodeIndex("4")], clusteringCoeffUndir(2*1, 2), epsilon);
		assertEquals(clust[testNet.getNodeIndex("5")], 0, epsilon);
		assertEquals(clust[testNet.getNodeIndex("6")], 0, epsilon);
		assertEquals(clust[testNet.getNodeIndex("7")], 0, epsilon);
	}

	
	// ----------------------------------------------------------------------------

	/** Directed clustering coefficient */
	@Test
	public void testClusteringCoeffDirected() {
		
		// Load undirected network with self-loops
		Network testNet = new Network("src/edu/mit/magnum/netprop/test/degreeTestNet.txt", true, false);
		BasicProperties test = new BasicProperties(testNet);

		// Compute betweenness
		test.computeClusteringCoefficient();
		Double[] clust = test.getClusteringCoeff();
		
		double epsilon = 1e-6;
		assertEquals(clust[testNet.getNodeIndex("1")], clusteringCoeffUndir(2, 4), epsilon);
		assertEquals(clust[testNet.getNodeIndex("2")], clusteringCoeffUndir(1, 2), epsilon);
		assertEquals(clust[testNet.getNodeIndex("3")], clusteringCoeffUndir(3, 3), epsilon);
		assertEquals(clust[testNet.getNodeIndex("4")], clusteringCoeffUndir(1, 2), epsilon);
		assertEquals(clust[testNet.getNodeIndex("5")], 0, epsilon);
		assertEquals(clust[testNet.getNodeIndex("6")], 0, epsilon);
		assertEquals(clust[testNet.getNodeIndex("7")], 0, epsilon);
	}


	// ============================================================================
	// PRIVATE METHODS

	/** Compute clustering coeff from number of edges between neighbors and number of neighbors */
	private double clusteringCoeffUndir(int numEdgesNeighbors, int numNeighbors) {
		
		return numEdgesNeighbors / (double)(numNeighbors*(numNeighbors-1));
	}

}
