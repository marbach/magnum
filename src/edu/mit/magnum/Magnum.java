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
 * Main class Gene Module Weaver (GMW) MAGNuM (Multi-scale Analysis of Gene
 * Network Modules) Multi-scale analysis of GWAS 
 * Adam (Adaptive DiseAse Modules)
 * Magneto: multi-scale module analysis of gene networks
 * (superhero theme)
 */
public class Magnum {

	// ============================================================================
	// STATIC METHODS

	/** Main function */
	static public void main(String[] args) {

		try {
			Magnum magnum = new Magnum(args);
			magnum.run();
		} catch (Exception e) {
			error(e);
		}
	}

	// ----------------------------------------------------------------------------

	/** Print the stack trace of the exception and exit */
	static public void error(Exception e) {
		e.printStackTrace();
		System.exit(-1); // return -1 in case of error
	}

	/** Print the stack trace of the exception and exit with error message */
	static public void error(Exception e, String msg) {
		System.err.println("ERROR: " + msg);
		error(e);
	}

	/** Print error message and exit */
	static public void error(String msg) {
		System.err.println("ERROR: " + msg);
		System.exit(-1);
	}

	/** Print warning message to stderr */
	static public void warning(String msg) {
		System.out.println("WARNING: " + msg);
		System.out.flush();
	}

	/** Write line to stdout */
	static public void println(String msg) {
		System.out.println(msg);
	}

	/** Write empty line to stdout */
	static public void println() {
		System.out.println();
	}

	/** Write string to stdout */
	static public void print(String msg) {
		System.out.print(msg);
	}
	
	/** Write string to stdout only if in verbose mode */
	static public void printlnVerbose(String msg) {
		if (Settings.verbose_)
			println(msg);
	}


	// ============================================================================
	// PUBLIC METHODS

	/** Constructor, parse command-line arguments, initialize settings */
	public Magnum(String[] args) {

		// Print magnum version
		printlnVerbose("magnum " + Settings.magnumVersion_ + "\n");
		
		// Parse command-line arguments and initialize settings
		MagnumOptionParser optionParser = new MagnumOptionParser();
		optionParser.parse(args);

		// Create output directory
		File outputDir = new File(Settings.outputDirectory_);
		if (!outputDir.exists())
			outputDir.mkdirs();
	}

	// ----------------------------------------------------------------------------

	/** Run */
	public void run() {

		if (Settings.mode_ == 1)
			runNetworkAnalysis();
		else if (Settings.mode_ == 2)
			runNetworkOperations();
		else if (Settings.mode_ == 3)
			runEnrichmentAnalysis();
		else if (Settings.mode_ == 4)
			runLinkModuleAnalysis();
		else {
			MagnumOptionParser.displayHelp();
			throw new IllegalArgumentException("--mode <int> must be between 1 and 3, found mode=" + Settings.mode_);
		}

		System.out.println("Success!");
	}

	// ----------------------------------------------------------------------------

	/** Analyze network properties */
	public void runNetworkAnalysis() {

		// Load the network
		if (Settings.computePstepKernel_) {
			Settings.isDirected_ = false;
			Settings.removeSelfLoops_ = true;
		}
		
		if (Settings.networkFile_ != null && !Settings.networkFile_.isEmpty()) {
			runNetworkAnalysis(Settings.networkFile_);
			
		} else {
			// List all files in the given directory
			ArrayList<String> networkFiles = MagnumUtils.listFiles(Settings.networkDir_);
			// Run for each file
			ArrayList<ArrayList<Double>> networkMeans = new ArrayList<ArrayList<Double>>();
			for (String file_i : networkFiles)
				networkMeans.add(runNetworkAnalysis(file_i));
			
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
	public ArrayList<Double> runNetworkAnalysis(String networkFile) {

		// Load the specified network file
		Network network = loadInputNetwork(networkFile);

		// Analyze network properties
		Magnum.println("COMPUTING NETWORK PROPERTIES");
		Magnum.println("----------------------------\n");

		NetpropMain netprop = new NetpropMain(network);
		ArrayList<Double> networkMeans = netprop.runAll();
		Magnum.println();
		netprop.saveAll();
		
		return networkMeans;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Perform network operations */
	public void runNetworkOperations() {

		// Perform operations on networks
		Magnum.println("PERFORMING NETWORK OPERATIONS");
		Magnum.println("-----------------------------\n");

		if (Settings.computeUnion_) {
			GroupNetworks grouper = new GroupNetworks(Settings.networkDir_, Settings.networkGroupFile_, Settings.networkFilePrefix_);
			grouper.run();

		} else if (Settings.computePairwiseSum_) {
			PairwiseSum networkSum = new PairwiseSum(Settings.networkDir_, Settings.networkDir2_);
			networkSum.run(true);
		}
	}

	
	// ----------------------------------------------------------------------------

	/** Analyze enrichment */
	public void runEnrichmentAnalysis() {

		// Load input files
		Magnum.println("LOADING INPUT FILES");
		Magnum.println("-------------------\n");

		EnrichMain enrich = new EnrichMain();
		// Ngsea.println();

		// Analyze enrichment
		Magnum.println("COMPUTING ENRICHMENT");
		Magnum.println("--------------------\n");

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
	private Network loadInputNetwork(String networkFile) {

		Magnum.println("LOADING INPUT NETWORK");
		Magnum.println("---------------------\n");

		return new Network(Settings.networkDir_ + "/" + networkFile, Settings.refNodesFile_,
				Settings.isDirected_, Settings.removeSelfLoops_,
				Settings.isWeighted_, Settings.threshold_);
	}

}
