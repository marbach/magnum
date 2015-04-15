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
package edu.mit.magnum.net;

import java.util.Comparator;


/**
 * A graph vertex / node
 */
public class Node {

	// The ID of the node
	protected String id_ = null;

	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public Node(String id) {

		id_ = id;
	}

	
	// ----------------------------------------------------------------------------

	/** Equality is defined by id_ only */
	public boolean equals(Object o) {
		return id_.equals(((Node) o).getId());
	}
	
	/** Equality is defined by id_ only */
	public int hashCode() {
		return id_.hashCode();
	}

	/** String representation, returns simply the id (JUNG uses this sometimes) */
	public String toString() {
		 return id_;
	} 
	
	


	
	// ============================================================================
	// PRIVATE METHODS
		
	// Used to sort nodes by ID
	static private class NodeComparator implements Comparator<Node> {
	    public int compare(Node n1, Node n2) {
	        return n1.getId().compareTo(n2.getId());
	    }
	}
	

	// ============================================================================
	// GETTERS AND SETTERS

	public String getId() { return id_; }
	
	static public NodeComparator getNodeComparator() { return new NodeComparator(); }
	
}
