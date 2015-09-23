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
package edu.mit.magnum;

import java.io.File;
import java.util.ArrayList;

import edu.mit.magnum.net.*;
import edu.mit.magnum.netops.*;
import edu.mit.magnum.netprop.NetpropMain;
import edu.mit.magnum.enrich.*;

/**
 * Main class 
 * Gene Module Weaver (GMW) 
 * MAGNuM (Multi-scale Analysis of Gene
 * Network Modules) Multi-scale analysis of GWAS 
 * Adam (Adaptive DiseAse Modules)
 * Magneto: multi-scale module analysis of gene networks
 * (superhero theme)
 */
public class Magnum {
	
	
	// ============================================================================
	// STATIC METHODS
	
	/** The logger -- a different logger can be plugged in for custom logging */
	static public MagnumLogger log = new MagnumLogger();
	/** The settings (needs to be after log so that we can print stuff!) */
	public static MagnumOptionParser set = new MagnumOptionParser();

	/** Main function */
	static public void main(String[] args) {

		try {
			Magnum magnum = new Magnum(args);
			magnum.run();
			
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	
	// ----------------------------------------------------------------------------

	/** Flag indicates if thread was interrupted */
	static private boolean interrupted_ = false;

	/** Gets and resets the interrupted flag (same behaviour as Thread.interrupted() */
	static public boolean interrupted() {
		boolean returnValue = interrupted_;
		interrupted_ = false;
		return returnValue;
	}

	/** Set interrupted flag true */
	static public void setInterrupted() {
		interrupted_ = true;
	}

	/** Throws a runtime exception if the thread has been interrupted */
	static public void exitOnInterrupt() {
		// Throw exception on interrupt
		if (Thread.interrupted()) {
			interrupted_ = true;
			throw new RuntimeException();
		}
	}
	
	/** Checks if there has been an interrupt, useful if cleanup has to be done before exciting */
	static public boolean checkInterrupt() {
		// Return true on interrupt
		if (Thread.interrupted())
			interrupted_ = true;
		return interrupted_;
	}

	
	// ============================================================================
	// PUBLIC METHODS

	/** Constructor without args */
	public Magnum() {
		this(null);
	}

	
	/** Constructor, parse command-line arguments, initialize settings */
	public Magnum(String[] args) {

		// Print magnum version
		Magnum.log.println("Magnum v" + set.magnumVersion + "\n");
		
		// Parse command-line arguments and initialize settings
		if (args != null)
			set.parse(args);
		else  // TODO hmmm....
			set.resetToDefaults();
	}

	
	// ----------------------------------------------------------------------------

	/** Run */
	public void run() {

		// Create output directory
		set.outputDirectory_.mkdirs();

		if (set.mode_ == 1)
			runNetworkAnalysis();
		else if (set.mode_ == 2)
			runNetworkOperations();
		else if (set.mode_ == 3)
			runEnrichmentAnalysis();
		else if (set.mode_ == 4)
			runLinkModuleAnalysis();
		else {
			MagnumOptionParser.displayHelp();
			throw new IllegalArgumentException("--mode <int> must be between 1 and 3, found mode=" + set.mode_);
		}

		Magnum.log.println("Success!");
	}

	// ----------------------------------------------------------------------------

	/** Analyze network properties */
	public void runNetworkAnalysis() {

		// Load the network
		if (set.computePstepKernel_) {
			set.isDirected_ = false;
			set.removeSelfLoops_ = true;
		}
		
		if (set.networkFile_ != null) {
			runNetworkAnalysis(set.networkFile_);
			
		} else {
			// List all files in the given directory
			ArrayList<String> networkFiles = MagnumUtils.listFiles(set.networkDir_);
			// Run for each file
			ArrayList<ArrayList<Double>> networkMeans = new ArrayList<ArrayList<Double>>();
			for (String file_i : networkFiles)
				networkMeans.add(runNetworkAnalysis(new File(file_i)));
			
			// Print means to file
			FileExport writer = new FileExport("networkMeans.txt");
			for (int n=0; n<networkFiles.size(); n++) {
				ArrayList<Double> means = networkMeans.get(n);
				if (means == null || means.size() == 0)
					continue;
				
				writer.print(networkFiles.get(n));
				for (int i=0; i<means.size(); i++)
					writer.print("\t" + means.get(i));
				writer.print("\n");
			}
			writer.close();
		}

	}
	
	
	// ----------------------------------------------------------------------------

	/** Analyze network properties */
	public ArrayList<Double> runNetworkAnalysis(File networkFile) {

		// Load the specified network file
		Network network = loadInputNetwork(networkFile);

		// Analyze network properties
		Magnum.log.println("COMPUTING NETWORK PROPERTIES");
		Magnum.log.println("----------------------------\n");

		NetpropMain netprop = new NetpropMain(network);
		ArrayList<Double> networkMeans = netprop.runAll();
		Magnum.log.println();
		netprop.saveAll();
		
		return networkMeans;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Perform network operations */
	public void runNetworkOperations() {

		// Perform operations on networks
		Magnum.log.println("PERFORMING NETWORK OPERATIONS");
		Magnum.log.println("-----------------------------\n");

		if (set.computeUnion_) {
			GroupNetworks grouper = new GroupNetworks(set.networkDir_, set.networkGroupFile_, set.networkFilePrefix_);
			grouper.run();

		} else if (set.computePairwiseSum_) {
			PairwiseSum networkSum = new PairwiseSum(set.networkDir_, set.networkDir2_);
			networkSum.run(true);
		}
	}

	
	// ----------------------------------------------------------------------------

	/** Analyze enrichment */
	public void runEnrichmentAnalysis() {

		// Load input files
		Magnum.log.println("LOADING INPUT FILES");
		Magnum.log.println("-------------------\n");

		EnrichMain enrich = new EnrichMain();

		// Analyze enrichment
		Magnum.log.println("COMPUTING ENRICHMENT");
		Magnum.log.println("--------------------\n");

		enrich.run();
	}


	// ----------------------------------------------------------------------------

	/** Analyze link module association */
	public void runLinkModuleAnalysis() {

//		// Load the network
//		Network network = loadInputNetwork();
//
//		Ngsea.println("LINK COMMUNITY IDENTIFICATION");
//		Ngsea.println("-----------------------------\n");
//
//		// Detect communities
//		LinkComm linkComm = new LinkCommPartitionImpl(network);
//		linkComm.runLinkCommunityDetection();
//
//		// Visualize results
//		// LinkCommViz.visualize(linkComm);
	}

		
	// ============================================================================
	// PRIVATE METHODS

	/** Load the input network */
	private Network loadInputNetwork(File networkFile) {

		Magnum.log.println("LOADING INPUT NETWORK");
		Magnum.log.println("---------------------\n");

		if (set.networkDir_ != null)
			networkFile = new File(set.networkDir_, networkFile.getPath());
		
		return new Network(networkFile, set.refNodesFile_,
				set.isDirected_, set.removeSelfLoops_,
				set.isWeighted_, set.threshold_);
	}

}
