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

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import org.apache.commons.collections15.BidiMap;
import org.apache.commons.collections15.bidimap.DualHashBidiMap;

import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import cern.colt.matrix.linalg.SeqBlas;
import edu.mit.magnum.FileExport;
import edu.mit.magnum.FileParser;
import edu.mit.magnum.Magnum;
import edu.mit.magnum.MagnumUtils;
import edu.uci.ics.jung.algorithms.util.Indexer;
import edu.uci.ics.jung.graph.AbstractTypedGraph;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import edu.uci.ics.jung.graph.util.Pair;


/**
 * Encapsulates a JUNG graph
 */
public class Network {

	/** The file from which this network was loaded */
	protected File file_ = null;
	
	/** The JUNG graph */
	protected AbstractTypedGraph<Node, Edge> graph_ = null;
	
	/** The nodes of the graph, indexed by their id label */
	protected HashMap<String, Integer> labelIndexMap_ = null;
	/** A bidirectional map assigning an index to each node (0, 1, ..., N-1) */
	protected BidiMap<Node, Integer> nodeIndexMap_ = null;
	/** A bidirectional map assigning an index to each prior node (0, 1, ..., M-1) */
	protected BidiMap<Node, Integer> refNodeIndexMap_ = null;
	
	/** Define an index for each edge */
	protected HashMap<Edge, Integer> edgeIndexMap_ = null;

	/** The number of nodes */
	protected int numNodes_ = -1;
	/** The number of reference nodes (default, all nodes of the network) */
	protected int numRefNodes_ = -1;
	/** True if reference nodes are given (numRefNodes_ < numNodes) */
	protected boolean useRefNodes_ = false;
	
	/** Defines whether the network is directed or undirected */
	protected boolean isDirected_ = true;
	/** Defines whether self-loops are allowed or removed when loading a network */	
	protected boolean removeSelfLoops_ = false;
	/** Defines whether edges are weighted */
	protected boolean isWeighted_ = false;
	/** Threshold for including weighted edges */
	protected double threshold_ = 0;
	
	/** The number of multi-edges that were removed when loading the network */
	protected int numRemovedMultiEdges_ = 0;
	/** The number of self-loops that were removed when loading the network */
	protected int numRemovedSelfEdges_ = 0;
	/** The number of edges that were discarded because their weight is below the threshold */
	protected int numBelowThreshold_ = 0;
	/** The number of super-hubs that were removed when loading the network */
	protected int numRemovedSuperHubs_ = 0;
	/** The number of isolated nodes (degree 0) that were removed */
	protected int numRemovedIsolatedNodes_ = 0;

	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Default constructor */
	public Network() { }

	/** Constructor loading network from file */
	public Network(File file, 
			boolean isDirected, boolean removeSelfLoops, boolean isWeighted, double threshold) {
		
		isDirected_ = isDirected;
		removeSelfLoops_ = removeSelfLoops;
		isWeighted_ = isWeighted;
		threshold_ = threshold;
		
		// Load the network
		loadNetwork(file);
	}

	/** Constructor loading network from file and specifying a set of reference nodes */
	public Network(File file, File refNodesFile, 
			boolean isDirected, boolean removeSelfLoops, boolean isWeighted, double threshold) {
		
		this(file, isDirected, removeSelfLoops, isWeighted, threshold);
		
		// Load the reference nodes
		if (refNodesFile != null)
			loadRefNodes(refNodesFile);
	}


	/** Constructor for unweighted networks with default threshold 0 */
	public Network(File file, 
			boolean isDirected, boolean removeSelfLoops) {
		
		this(file, isDirected, removeSelfLoops, false, 0);
	}

	
	// ----------------------------------------------------------------------------

	/** Load the Network */
	public void loadNetwork(File file) {
				
		// Create graph_
		this.file_ = file;
		loadGraph();
		// Remove super-hubs
		removeSuperHubs();
		// Remove nodes that became isolated after removing super-hubs
		removeIsolatedNodes();
		
		// Initialize nodeLabelMap_, nodeIndexMap_, priorNodeIndexMap_
		initializeNodeMaps();
		// Initialize edgeIndexMap_
		initializeEdgeIndexMap();
		
		// Initialize number of nodes and reference nodes
		numNodes_ = graph_.getVertexCount();
		numRefNodes_ = refNodeIndexMap_.size();
		useRefNodes_ = false;
		
		// Print info
		Magnum.log.println("- Treat as: " + (isDirected_ ? "DIRECTED" : "UNDIRECTED") + ", " + (isWeighted_ ? "WEIGHTED" : "UNWEIGHTED"));
		Magnum.log.println("- Remove self-loops: " + (removeSelfLoops_ ? "YES" : "NO"));
		Magnum.log.println("- Discard edges below threshold: " + threshold_);
		Magnum.log.println("");
		
		Magnum.log.println("Loaded network with:");
		Magnum.log.println("- " + graph_.getVertexCount() + " nodes");
		Magnum.log.println("- " + graph_.getEdgeCount() + " edges");
		Magnum.log.println("Removed:");
		Magnum.log.println("- " + numBelowThreshold_ + " edges below threshold");			
		Magnum.log.println("- " + numRemovedMultiEdges_ + " multi-edges" + (isWeighted_? ", taking MAX edge weight" : "") );
		if (removeSelfLoops_)
			Magnum.log.println("- " + numRemovedSelfEdges_ + " self-loops");
		if (Magnum.set.superHubThreshold_ > 0 && Magnum.set.superHubThreshold_ < 1)
			Magnum.log.println("- " + numRemovedSuperHubs_ + " super-hubs connecting > " + 100*Magnum.set.superHubThreshold_ + "% of all genes");
		if (numRemovedIsolatedNodes_ > 0)
			Magnum.log.println("- " + numRemovedIsolatedNodes_ + " isolated nodes (degree 0)");
		Magnum.log.println("");
	}

	
	// ----------------------------------------------------------------------------

	/** Save the Network */
	public void write(String filename) {

		FileExport writer = new FileExport(filename, true);
		for (Edge edge : graph_.getEdges()) {
			String nextLine = graph_.getSource(edge).id_ + "\t" + graph_.getDest(edge).id_;
			if (isWeighted_)
				nextLine += "\t" + MagnumUtils.toStringScientific10(edge.w_);
			writer.println(nextLine);
		}
		writer.close();
	}
	
		
	// ----------------------------------------------------------------------------

	/** Load the set of reference nodes */
	public void loadRefNodes(File file) {
		
		// Open the file
		FileParser parser = new FileParser(file);
		String[] nextLine = parser.readLine();
		refNodeIndexMap_ = new DualHashBidiMap<Node,Integer>();
		
		// For each node
		while (nextLine != null) {
			if (nextLine.length != 1)
				throw new RuntimeException("Line " + parser.getLineCounter() + " has " + nextLine.length + " columns (file format is one node ID per line)");
			
			// Get the node
			String id = nextLine[0];
			Node node = getNode(id);
			if (node == null)
				throw new RuntimeException("Line " + parser.getLineCounter() + ": node '" + id + "' is not part of the network");
			
			if (refNodeIndexMap_.containsKey(node))
				throw new RuntimeException("Line " + parser.getLineCounter() + ": node '" + id + "' is listed multiple times");
			
			// Add the ref node
			refNodeIndexMap_.put(node, parser.getLineCounter()-1);
			
			// Read the next line
			nextLine = parser.readLine();
		}
		assert parser.getLineCounter()-1 == refNodeIndexMap_.size();
		
		numRefNodes_ = refNodeIndexMap_.size();
		useRefNodes_ = true;
	}
		
		
	// ----------------------------------------------------------------------------

	/** Get the neighbors of the given vertex, do *not* include the node itself in case there is a self-loop */
	public Collection<Node> getNeighborsNoSelf(Node node) {
		
		// Get the neighbors, may include the given node itself if there is a self-loop
		Collection<Node> neighbors = graph_.getNeighbors(node);
		
		// Remove node if there is a self-loop
		if (graph_.isNeighbor(node, node)) {
			neighbors = new HashSet<Node>(neighbors);
			boolean found = neighbors.remove(node);
			assert(found);
		}
		return neighbors;
	}

	
	// ----------------------------------------------------------------------------

	/** Get the edge between the two nodes, return null if there is no such edge */
	public Edge getEdge(String reg, String tar) {

		return graph_.findEdge(new Node(reg), new Node(tar));
	}

	
	// ----------------------------------------------------------------------------

	/** Define a set of reference nodes */
	public void setRefNodes(Collection<Node> refNodes) {
					
		for (Node n : refNodes)
			if (!graph_.containsVertex(n))
				throw new IllegalArgumentException("The specified reference node (" + n.getId() + ") is not part of the network");
		
		numRefNodes_ = refNodes.size();
		refNodeIndexMap_ = Indexer.create(refNodes);
		if (numRefNodes_ < numNodes_)
			useRefNodes_ = true;
	}
	

	// ----------------------------------------------------------------------------

	/** Get the set of reference nodes */
	public Collection<Node> getRefNodes() {
		return refNodeIndexMap_.keySet();
	}

	/** Return true if this node is part of the reference node set */
	public boolean isRefNode(Node node) {
		return refNodeIndexMap_.containsKey(node);
	}

	
	// ----------------------------------------------------------------------------

	/** Get the sum of incident edge weights (i.e., the degree for unweighted networks) */
	public double weightedDegree(Node node) {
	
		double wdegree = 0;
		for (Edge e : graph_.getIncidentEdges(node)) {
			if (isWeighted_)
				wdegree += e.w_;
			else
				wdegree++;
		}
		return wdegree;
	}
	
	
	// ----------------------------------------------------------------------------

    /**
     * Returns the adjacency matrix A of the graph. A_ij corresponds to the edge
     * from vertex i to vertex j. Adapted from JUNG graphToSparseMatrix().
     * For weighted networks, the weighted adjacency matrix is given.
     */
    public SparseDoubleMatrix2D computeAdjacencyMatrix() {
    	
        SparseDoubleMatrix2D A = new SparseDoubleMatrix2D(numNodes_, numNodes_);
        
        for (int i=0; i<numNodes_; i++) {
        	Node node_i = getNode(i);
        	
            for (Edge e : graph_.getOutEdges(node_i)) {
                Node node_j = graph_.getOpposite(node_i, e);
                int j = getNodeIndex(node_j);
                if (isWeighted_)
                	A.set(i, j, e.w_); // JUNG uses getQuick()
                else
                	A.set(i, j, 1); // JUNG uses getQuick()
            }
        }
        return A;
    }

    
	// ----------------------------------------------------------------------------

    /**
     * Compute the Laplacian (only defined for simple graphs, i.e. no self-loops allowed):
     *         / degree(v_i)  if i==j
     * l_ij = |  -1           if i!=j and there is an edge (v_i, v_j)
     *         \ 0            otherwise
     */
	public SparseDoubleMatrix2D computeLaplacian() {
		
		// L = A
		SparseDoubleMatrix2D L = computeAdjacencyMatrix();
		// L: numNodes_ x numNodes_
		if (L.rows() != numNodes_ || L.columns() != numNodes_)
			throw new RuntimeException("Unexpected number of rows/cols in adjacency matrix");

		// L = -A
		SeqBlas.seqBlas.dscal(-1, L);

		// l_ii = degree(v_i)
		for (int i=0; i<numNodes_; i++) {
			Node n = getNode(i);
			if (graph_.isNeighbor(n, n))
				throw new RuntimeException("Laplacian is defined for simple graphs, no self-loops allowed (rerun with option removeSelfLoops)");
			
			L.set(i, i, weightedDegree(n));
		}
		return L;
	}

	
	// ----------------------------------------------------------------------------

    /**
     * Compute the normalized Laplacian:
     *         / 1                          if i==j and degree(v_i)!=0
     * l_ij = | -w_ij/sqrt(deg(v_i)*deg(v_j))  if i!=j and there is an edge (v_i, v_j)
     *         \ 0                          otherwise 
     */
	public SparseDoubleMatrix2D computeNormalizedLaplacian() {

		// The normalized Laplacian
		SparseDoubleMatrix2D Lnorm = new SparseDoubleMatrix2D(numNodes_, numNodes_);
		// The weighted degree of each node
		double[] weightedDegree = new double[numNodes_];
		
        // Case 1
		boolean isolatedNodes = false;
        for (int i=0; i<numNodes_; i++) {
        	Node n = getNode(i);
			if (graph_.isNeighbor(n, n))
				throw new RuntimeException("Laplacian is defined for simple graphs, no self-loops allowed (rerun with option removeSelfLoops)");

			if (graph_.degree(n) != 0)
        		Lnorm.set(i, i, 1);
        	else
        		isolatedNodes = true;
			
			// Pre-compute the weighted degree (needed below)
			weightedDegree[i] = weightedDegree(n);
        }
        if (isolatedNodes)
        	Magnum.log.warning("Isolated nodes in network (degree 0)");

		// Case 2
        for (Edge edge : graph_.getEdges()) {
        	Pair<Node> nodes = graph_.getEndpoints(edge);
        	int i = getNodeIndex(nodes.getFirst());
        	int j = getNodeIndex(nodes.getSecond());
        	
        	// Note, since the vertices are connected their degree cannot be 0
        	double l_ij = -edge.w_ / Math.sqrt(weightedDegree[i]*weightedDegree[j]);
        	Lnorm.set(i, j, l_ij); // JUNG uses setQuick()
        	Lnorm.set(j, i, l_ij); // JUNG uses setQuick()
        }		
		return Lnorm;
	}

	
	// ----------------------------------------------------------------------------

    /** Get the regulators (nodes with outdegree > 0) */
	public Collection<Node> getRegulatorNodes() {
		
		HashSet<Node> regulators = new HashSet<Node>();
		for (Node node : graph_.getVertices())
			if (graph_.outDegree(node) > 0)
				regulators.add(node);
		
		return regulators;
	}

	
    /** Get the targets (nodes with indegree > 0) */
	public Collection<Node> getTargetNodes() {
		
		HashSet<Node> targets = new HashSet<Node>();
		for (Node node : graph_.getVertices())
			if (graph_.inDegree(node) > 0)
				targets.add(node);
		
		return targets;
	}


    
	// ============================================================================
	// PRIVATE METHODS
		
	/** 
	 * Read the network from the given file, create JUNG graph instance.
	 * Returns the number of removed multi-edges and self-edges.
	 */
	private void loadGraph() {
		
		// New graph instance
		if (isDirected_)
			graph_ = new DirectedSparseGraph<Node, Edge>();
		else
			graph_ = new UndirectedSparseGraph<Node, Edge>();
		
		// Open the file
		FileParser parser = new FileParser(file_);
		if (Magnum.set.networkFileDelim_.equalsIgnoreCase("TAB"))
			parser.setSeparator("\t");
		else if (Magnum.set.networkFileDelim_.equalsIgnoreCase("SPACE"))
			parser.setSeparator(" ");
		else
			throw new IllegalArgumentException("Settings.networkFileDelim_ must be either 'tab' or 'space' (in words like this)");
			
		// Read first line
		String[] nextLine = parser.readLine();
		// Check number of columns
		int numCol = nextLine.length;
		if (isWeighted_ && numCol != 3)
			throw new RuntimeException("Weighted network must have 3 columns");
		else if (!isWeighted_ && (numCol < 2 || numCol > 3))
			throw new RuntimeException("Network file must have 2 or 3 columns");
		
		// The number of multi / self edges that were removed
		numRemovedMultiEdges_ = 0;
		numRemovedSelfEdges_ = 0;
		numBelowThreshold_ = 0;
		
		HashMap<String, Node> nodes = new HashMap<String, Node>();
		
		// For each line / edge
		while (nextLine != null) {
			if (nextLine.length != numCol)
				throw new RuntimeException("Line " + parser.getLineCounter() + " has " + nextLine.length + " columns");
			
			// Create the edge
			Edge edge;
			if (isWeighted_) {
				// Parse the weight
				double w = Double.parseDouble(nextLine[2]);
				// Skip if the weight is under the threshold
				if (w < threshold_) {
					numBelowThreshold_++;
					nextLine = parser.readLine();
					continue;
				}
				edge = new Edge(w);
			} else {
				edge = new Edge();
			}
			
			// The regulator and target node label
			String regId = nextLine[0];
			String tarId = nextLine[1];

			// Skip self loops if they are not allowed
			if (removeSelfLoops_ && regId.equals(tarId)) {
				numRemovedSelfEdges_++;
				nextLine = parser.readLine();
				continue;
			}
			
			// Get the nodes if they already exist (avoids creating new ones when calling
			// graph_.addEdge() below, it actually does although it doesn't seem to be a bug)
			Node reg = nodes.get(regId);
			Node tar = nodes.get(tarId);
			if (reg == null) {
				reg = new Node(regId);
				nodes.put(regId, reg);
			}
			if (tar == null) {
				tar = new Node(tarId);
				nodes.put(tarId, tar);
			}
			
			// Add the edge (also adds the nodes if they are new)
			// Returns false if the edge already exists
			boolean edgeAdded = graph_.addEdge(edge, reg, tar);
			if (!edgeAdded) {
				if (isWeighted_) {
					// For weighted networks, assign the max if there are multiple edges
					Edge existingEdge = getEdge(regId, tarId);
					if (existingEdge.w_ < edge.w_)
						existingEdge.w_ = edge.w_;
					//throw new RuntimeException("Multi-edges are not handled for weighted networks, remove multi-edges from file");
				}
				numRemovedMultiEdges_++;
			}
			
			// Return if thread was interrupted
			if (Thread.interrupted()) {
				parser.close();
				Magnum.setInterrupted();
				throw new RuntimeException();
			}
			
			// Read the next line
			nextLine = parser.readLine();
		}
		parser.close();
		assert nodes.size() == graph_.getVertexCount();
	}
	
	
	// ----------------------------------------------------------------------------

	/** Exclude "super-hubs" that connect to more than the given fraction of genes */
	private void removeSuperHubs() {
				
		if (Magnum.set.superHubThreshold_ <= 0 || Magnum.set.superHubThreshold_ >= 1)
			return;
		
		int maxDegree = (int) (Magnum.set.superHubThreshold_ * graph_.getVertexCount());
		ArrayList<Node> superHubs = new ArrayList<Node>();
		
		for (Node node : graph_.getVertices())
			if (graph_.degree(node) >= maxDegree)
				superHubs.add(node);

		numRemovedSuperHubs_ = superHubs.size();
		
		for (Node hub : superHubs)
			graph_.removeVertex(hub);
	}

	
	// ----------------------------------------------------------------------------

	/** Remove disconnected nodes */
	private void removeIsolatedNodes() {
				
		ArrayList<Node> isolatedNodes = new ArrayList<Node>();
		for (Node node : graph_.getVertices())
			if (graph_.degree(node) == 0)
				isolatedNodes.add(node);

		numRemovedIsolatedNodes_ = isolatedNodes.size();
		
		for (Node island : isolatedNodes)
			graph_.removeVertex(island);
	}

	
	// ----------------------------------------------------------------------------

	/** Put all the nodes from the graph into the nodes_ hash map */
	private void initializeNodeMaps() {
		
		// Create Node-Index bidi map
		nodeIndexMap_ = Indexer.create(graph_.getVertices());
		
		// Create a sorted list of nodes
		ArrayList<Node> nodes = new ArrayList<Node>(graph_.getVertices());
		Collections.sort(nodes, Node.getNodeComparator());
		
		// Create the index map
		nodeIndexMap_ = new DualHashBidiMap<Node,Integer>();
		for (int i=0; i<nodes.size(); i++)
			nodeIndexMap_.put(nodes.get(i), i);
		
		// All nodes are reference nodes by default
		refNodeIndexMap_ = nodeIndexMap_;

		// Create Label-Index map
		//nodeLabelMap_ = new HashMap<String, Node>();
		labelIndexMap_ = new HashMap<String, Integer>();
		
		// Loop over the nodes, put then into the map
		Iterator<Node> iter = graph_.getVertices().iterator();
		while (iter.hasNext()) {
			Node n = iter.next();
			//nodeLabelMap_.put(n.getId(), n);
			labelIndexMap_.put(n.getId(), nodeIndexMap_.get(n));
		}
		
		//assert(nodeLabelMap_.size() == graph_.getVertexCount());
		assert(nodeIndexMap_.size() == graph_.getVertexCount());
	}

	// ----------------------------------------------------------------------------

	/** Initialize edge indexes */
	private void initializeEdgeIndexMap() {
		
		edgeIndexMap_ = new HashMap<Edge, Integer>(graph_.getEdgeCount());
		int index = 0;
		for (Edge e : graph_.getEdges())
			edgeIndexMap_.put(e, index++);
	}

	
	// ============================================================================
	// SETTERS AND GETTERS
	
	public File getFile() { return file_; }
	
	public int getNumRegulators() {
		if (!isDirected_)
			return 0;
		else
			return getRegulatorNodes().size();
	}
	
	public int getNumNodes() { return numNodes_; }
	public int getNumEdges() { return graph_.getEdgeCount(); }
	public int getNumRefNodes() { return numRefNodes_; }
	public boolean getUseRefNodes() { return useRefNodes_; }
	public boolean getIsDirected() { return isDirected_; }
	public boolean getIsWeighted() { return isWeighted_; }
	
	public AbstractTypedGraph<Node, Edge> getGraph() { return graph_; }
	
	public Node getNode(String id) { return getNode(getNodeIndex(id)); }
	
	public Integer getNodeIndex(String id) { return labelIndexMap_.get(id); }
	public Integer getNodeIndex(Node n) { return nodeIndexMap_.get(n); }
	public Node getNode(Integer i) { return nodeIndexMap_.getKey(i); }
	
	public Integer getEdgeIndex(Edge e) { return edgeIndexMap_.get(e); }

	//public Integer getRefNodeIndex(String id) { return labelIndexMap_.get(id); }
	public Integer getRefNodeIndex(Node n) { return refNodeIndexMap_.get(n); }
	public Node getRefNode(Integer i) { return refNodeIndexMap_.getKey(i); }

}
