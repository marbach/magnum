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
package edu.mit.magnum.netops.test;

import static org.junit.Assert.*;

import java.io.File;

import org.junit.*;

import edu.mit.magnum.Magnum;
import edu.mit.magnum.net.Network;
import edu.mit.magnum.netops.Union;


/**
 * Unit tests for NetworkUnion
 */
public class UnionTest {
	
	/** The magnum instance */
	private static Magnum mag = new Magnum();

	// ============================================================================
	// SETUP
	
	@BeforeClass
	public static void testSetup() {
		mag.set.resetToDefaults();
		mag.set.superHubThreshold_ = 0;
		mag.set.computeUnion_ = true;
		mag.set.isDirected_ = true;
		mag.set.isWeighted_ = true;
	}

	@AfterClass
	public static void testCleanup() {
	}
	  
	// ============================================================================
	// TESTS

	/** Node degrees */
	@Test
	public void testUnionMax() {

		// Initialize netops
		Union union = new Union(mag, new File("src/edu/mit/magnum/netops/test/net.e"));
		// Compute union
		Network net = union.run();
		assertEquals(8, net.getNumEdges());
		
		double eps = 1e-12;
		assertEquals(0.2, net.getEdge("r1", "r2").w_, eps); 
		assertEquals(0.1, net.getEdge("r1", "g1").w_, eps); 
		assertEquals(0.3, net.getEdge("r2", "r1").w_, eps); 
		assertEquals(0.4, net.getEdge("r2", "g1").w_, eps); 
		assertEquals(0.8, net.getEdge("r2", "g2").w_, eps); 
		assertEquals(1, net.getEdge("r2", "g3").w_, eps); 
		assertEquals(0.6, net.getEdge("r3", "g1").w_, eps); 
		assertEquals(0.7, net.getEdge("r3", "g4").w_, eps); 
	}

	
	// ----------------------------------------------------------------------------

	// ============================================================================
	// PRIVATE METHODS

}
