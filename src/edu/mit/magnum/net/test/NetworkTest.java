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
package edu.mit.magnum.net.test;

import static org.junit.Assert.*;

import org.junit.*;

import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import edu.mit.magnum.MagnumSettings;
import edu.mit.magnum.net.*;


/**
 * Unit tests for Network.
 * 
 * Test network: data/test/ppi.txt
 * Nodes: 13977
 * Edges: 195444
 * Multi-edges: 96967
 * Single-edges: 98477
 * Self-loops: 1510
 * Single-edge, no self loops: 96967
 * Nodes with at least one edge after removing self-loops: 13942
 */
public class NetworkTest {
	
	/** Test network */
	private String ppiNetFile_ = "src/edu/mit/magnum/netprop/test/ppi.txt";
	
	/** Network properties of test network */
	private int ppiNumNodes_ = 13942;
	private int ppiNumUndirEdges_ = 96967;
	private int ppiNumNodesAllowSelfLoops_ = 13977;
	private int ppiNumUndirEdgesAllowSelfLoops_ = 98477;
	
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

	/** Undirected load */
	@Test
	public void testUndirectedLoad() {

		// Create undirected network without self-loops
		Network ppi = new Network(ppiNetFile_, false, true);
		assertEquals(ppi.getNumNodes(), ppiNumNodes_);
		assertEquals(ppi.getNumEdges(), ppiNumUndirEdges_);
		
		// Create undirected network with self-loops
		ppi = new Network(ppiNetFile_, false, false);
		assertEquals(ppi.getNumNodes(), ppiNumNodesAllowSelfLoops_);
		assertEquals(ppi.getNumEdges(), ppiNumUndirEdgesAllowSelfLoops_);
	}

	
	// ----------------------------------------------------------------------------

	/** Weighted load with threshold */
	@Test
	public void testWeightedLoad() {

		// Create undirected network without self-loops
		Network net = new Network("src/edu/mit/magnum/netprop/test/hierarchicalScaleFreeLevel1.txt", false, true, true, 0.5);		
		assertEquals(net.getNumNodes(), 4);
		assertEquals(net.getNumEdges(), 3);	
		
		double epsilon = 1e-12;
		assertEquals(0.8494727767538279, net.getEdge("0", "1").w_, epsilon);
		assertEquals(0.95, net.getEdge("0", "2").w_, epsilon);
		assertEquals(0.96, net.getEdge("0", "3").w_, epsilon);
		assertEquals(0.8494727767538279, net.getEdge("1", "0").w_, epsilon);
		assertEquals(0.95, net.getEdge("2", "0").w_, epsilon);
		assertEquals(0.96, net.getEdge("3", "0").w_, epsilon);

		// Create directed network without self-loops
		net = new Network("src/edu/mit/magnum/netprop/test/hierarchicalScaleFreeLevel1.txt", true, true, true, 0.5);		
		assertEquals(net.getNumNodes(), 4);
		assertEquals(net.getNumEdges(), 4);	
		
		assertEquals(0.8494727767538279, net.getEdge("0", "1").w_, epsilon);
		assertEquals(0.95, net.getEdge("0", "2").w_, epsilon);
		assertEquals(0.9002343488391489, net.getEdge("0", "3").w_, epsilon);
		assertNull(net.getEdge("1", "0"));
		assertNull(net.getEdge("2", "0"));
		assertEquals(0.96, net.getEdge("3", "0").w_, epsilon);
	}

	
	// ----------------------------------------------------------------------------

	/** Laplacian */
	@Test
	public void testComputeLaplacian() {
		
		Network net = new Network("src/edu/mit/magnum/netprop/test/simpleNet.txt", false, true);
		SparseDoubleMatrix2D L = net.computeLaplacian();
		
		double[][] L_expected = {
				{1, -1, 0, 0, 0, 0},
				{-1, 2, -1, 0, 0, 0},
				{0, -1, 2, -1, 0, 0},
				{0, 0, -1, 3, -1, -1},
				{0, 0, 0, -1, 2, -1},
				{0, 0, 0, -1, -1, 2}};
		
		double epsilon = 1e-12;
		for (int i=0; i<net.getNumNodes(); i++)
			for (int j=0; j<net.getNumNodes(); j++)
				assertEquals(L_expected[i][j], L.get(i, j), epsilon);
	}

	
	// ----------------------------------------------------------------------------

	/** Normalized Laplacian */
	@Test
	public void testComputeNormalizedLaplacian() {
		
		Network net = new Network("src/edu/mit/magnum/netprop/test/simpleNet.txt", false, true);
		SparseDoubleMatrix2D L = net.computeNormalizedLaplacian();
		
		double[][] L_expected = {
				{1.0000000, -0.7071068,  0.0000000,  0.0000000,  0.0000000,  0.0000000},
				{-0.7071068,  1.0000000, -0.5000000,  0.0000000,  0.0000000,  0.0000000},
				{0.0000000, -0.5000000,  1.0000000, -0.4082483,  0.0000000,  0.0000000},
				{0.0000000,  0.0000000, -0.4082483,  1.0000000, -0.4082483, -0.4082483},
				{0.0000000,  0.0000000,  0.0000000, -0.4082483,  1.0000000, -0.5000000},
				{0.0000000,  0.0000000,  0.0000000, -0.4082483, -0.5000000,  1.0000000}};
				
		double epsilon = 1e-6;
		for (int i=0; i<net.getNumNodes(); i++)
			for (int j=0; j<net.getNumNodes(); j++)
				assertEquals(L_expected[i][j], L.get(i, j), epsilon);
	}

	
}
