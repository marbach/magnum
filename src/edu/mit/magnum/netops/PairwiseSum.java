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

import java.io.File;
import java.util.ArrayList;

import edu.mit.magnum.*;
import edu.mit.magnum.net.*;
import edu.uci.ics.jung.graph.AbstractTypedGraph;


/**
 * Abstract class, extended by classes that compute node / network properties:
 * - NodeProperties
 */
public class PairwiseSum {

	/** The magnum instance */
	private Magnum mag;

	/** The first network directory */
	private File networkDir1_ = null;
	/** The second network directory */
	private File networkDir2_ = null;
	
	/** The input networks */
	private ArrayList<String> networkFiles1_ = null;
	/** The input networks */
	private ArrayList<String> networkFiles2_ = null;
	/** The sample IDs (files are given in this order) */
	private String[] sampleIds_ = null;
	
	/** Number of files */
	private int numFiles_ = -1;
	
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public PairwiseSum(Magnum mag, File networkDir1, File networkDir2) {
		
		this.mag = mag;
		networkDir1_ = networkDir1;
		networkDir2_ = networkDir2;
		
		initializeNetworkFiles();
	}
	
	
	// ----------------------------------------------------------------------------

	/**  */
	public Network[] run(boolean writeFiles) {
		
		if (!mag.set.isDirected_)
			throw new RuntimeException("Not yet implemented for undirected networks");
		if (!mag.set.isWeighted_)
			throw new RuntimeException("Not yet implemented for unweighted networks");
		
		// Output file prefix
		String outPrefix = mag.set.outputDirectory_ + System.getProperty("file.separator") + mag.set.outputFilename_ + ".";
		
		// Initialize array for results
		Network[] networks = null;
		if (!writeFiles)
			networks = new Network[numFiles_];
		
		for (int i=0; i<numFiles_; i++) {
			Network net1 = new Network(mag, new File(networkDir1_, networkFiles1_.get(i)), mag.set.isDirected_, mag.set.removeSelfLoops_, mag.set.isWeighted_, mag.set.threshold_);
			Network net2 = new Network(mag, new File(networkDir2_, networkFiles2_.get(i)), mag.set.isDirected_, mag.set.removeSelfLoops_, mag.set.isWeighted_, mag.set.threshold_);
			
			AbstractTypedGraph<Node, Edge> graph1 = net1.getGraph();
			AbstractTypedGraph<Node, Edge> graph2 = net2.getGraph();
			
			// Add edges of net2 to net1
			for (Edge edge2 : graph2.getEdges()) {
				Node tf = graph2.getSource(edge2);
				Node gene = graph2.getDest(edge2);
				Edge edge1 = graph1.findEdge(tf, gene);
				
				if (edge1 == null)
					graph1.addEdge(edge2, tf, gene);
				else
					edge1.w_ += edge2.w_;
			}
			
			// Write or save the file
			if (writeFiles)
				net1.write(outPrefix + sampleIds_[i] + ".txt");
			else
				networks[i] = net1;
		}
		
		return networks;
	}
	
	
	// ============================================================================
	// PRIVATE METHODS

	/** Initialize file names, check that pairwise networks are listed at corresponding positions */
	private void initializeNetworkFiles() {

		networkFiles1_ = mag.utils.listFiles(networkDir1_);
		networkFiles2_ = mag.utils.listFiles(networkDir2_);
		
		numFiles_ = networkFiles1_.size();		
		if (numFiles_ != networkFiles2_.size())
			throw new RuntimeException("Inconsistent number of files in directory");

		// Initialize sample IDs
		sampleIds_ = new String[numFiles_];
		for (int i=0; i<numFiles_; i++) {
			
			String name1 = mag.utils.extractBasicFilename(networkFiles1_.get(i), false);
			String name2 = mag.utils.extractBasicFilename(networkFiles2_.get(i), false);
			
			int index1 = name1.lastIndexOf(".");
			int index2 = name2.lastIndexOf(".");
			
			String sample1 = name1.substring(index1+1);
			String sample2 = name2.substring(index2+1);
			
			if (!sample1.equals(sample2))
				throw new RuntimeException("File pairs don't match");

			sampleIds_[i] = sample1;
		}
	}
}
