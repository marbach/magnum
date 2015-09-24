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
package edu.mit.magnum.netprop;

import java.util.ArrayList;
import java.util.LinkedHashMap;

import edu.uci.ics.jung.graph.AbstractTypedGraph;
import edu.mit.magnum.Magnum;
import edu.mit.magnum.net.*;


/**
 * Abstract class, extended by classes that compute node / network properties:
 * - NodeProperties
 */
public abstract class NetworkProperties {

	/** The magnum instance */
	protected Magnum mag;

	/** The network that is being analyzed */
	protected Network network_ = null;
	/** The JUNG graph of the network */
	protected AbstractTypedGraph<Node, Edge> graph_ = null;
	
	/** Defines whether the network is directed or undirected */
	protected boolean isDirected_ = true;
	/** The number of nodes in the network */
	protected int numNodes_ = -1;
	/** The number of reference nodes in the network */
	protected int numRefNodes_ = -1;
	
	/** Computed node properties such as degree, betweenness, etc. */
	protected LinkedHashMap<String, Number[]> nodeProperties_ = null;
		
	
	// ============================================================================
	// ABSTRACT METHODS

	/** Compute all metrics specified in Settings */
	public abstract ArrayList<Double> run();
		
	/** Add the computed node properties of this analyzer to the given map */
	public void addNodeProperties(LinkedHashMap<String,Number[]> map) { }

	/** Save all pairwise node properties of this subclass */
	public void saveK() { }

	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public NetworkProperties(Magnum mag, Network network) {
		
		this.mag = mag;
		network_ = network;
		graph_ = network.getGraph();
		isDirected_ = network.getIsDirected();
		numNodes_ = graph_.getVertexCount();
		numRefNodes_ = network.getNumRefNodes();
		
		nodeProperties_ = new LinkedHashMap<String, Number[]>();
	}
	
	
	// ============================================================================
	// PRIVATE METHODS


	
}
