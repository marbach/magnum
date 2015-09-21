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
package edu.mit.magnum.netops;

import java.util.ArrayList;

import edu.mit.magnum.MagnumUtils;
import edu.mit.magnum.Magnum;
import edu.mit.magnum.net.*;
import edu.uci.ics.jung.graph.AbstractTypedGraph;


/**
 * Abstract class, extended by classes that compute node / network properties:
 * - NodeProperties
 */
public class Union {

	/** The network directory */
	String networkDir_ = null;
	/** The input networks */
	ArrayList<String> networkFiles_ = null;
	/** The result network */
	Network network_ = null;

	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public Union(String networkDir) {
		
		networkDir_ = networkDir;
		networkFiles_ = MagnumUtils.listFiles(networkDir);
	}
	
	
	/** Constructor */
	public Union(String networkDir, ArrayList<String> networkFiles) {
		
		networkDir_ = networkDir;
		networkFiles_ = networkFiles;
	}

	
	// ----------------------------------------------------------------------------

	/** Compute union (max edges) across given network files */
	public Network run() {
		
		if (!Magnum.set.isDirected_)
			throw new RuntimeException("Not yet implemented for undirected networks");

		// Load the first network
		network_ = new Network(networkDir_ + "/" + networkFiles_.get(0), 
				Magnum.set.isDirected_, Magnum.set.removeSelfLoops_, Magnum.set.isWeighted_, Magnum.set.threshold_);
		
		// Add one network after the other
		for (int i=1; i<networkFiles_.size(); i++) {
			Network nextNet = new Network(networkDir_ + "/" + networkFiles_.get(i), 
					Magnum.set.isDirected_, Magnum.set.removeSelfLoops_, Magnum.set.isWeighted_, Magnum.set.threshold_);
			union(nextNet);
		}
		
		return network_;
	}
	
	
	// ============================================================================
	// PRIVATE METHODS

	/** Take union of network_ with the this network */
	private void union(Network nextNet) {

		AbstractTypedGraph<Node, Edge> graph = network_.getGraph();
		AbstractTypedGraph<Node, Edge> nextGraph = nextNet.getGraph(); 
				
		for (Edge nextEdge : nextGraph.getEdges()) {
			Node tf = nextGraph.getSource(nextEdge);
			Node gene = nextGraph.getDest(nextEdge);
			Edge edge = graph.findEdge(tf, gene);
			if (edge == null)
				graph.addEdge(nextEdge, tf, gene);
			else
				edge.w_ = Math.max(edge.w_, nextEdge.w_);
		}
	}

}

