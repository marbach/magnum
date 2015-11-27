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

import ch.unil.gpsutils.FileExport;
import ch.unil.gpsutils.Logger;
import ch.unil.gpsutils.Utils;
import edu.mit.magnum.net.*;
import edu.mit.magnum.netops.*;
import edu.mit.magnum.netprop.NetpropMain;
import edu.mit.magnum.enrich.*;
import edu.mit.magnum.experiments.Experiments;

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
	
	/** Current version */
	final public static String version = "1.0";

	/** The logger -- a different logger can be plugged in for custom logging */
	public Logger log;
	/** The settings */
	public MagnumOptionParser set;
	/** The utilities */
	public Utils utils;
	
	/** Connectivity enrichment analysis */
	private EnrichMain enrichMain;
	

	/** Main function */
	static public void main(String[] args) {

		try {
			Magnum magnum = new Magnum(args, null);
			magnum.run();
			
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	
	// ============================================================================
	// PUBLIC METHODS

	/** Constructor with custom logger */
	public Magnum() {
		this(null, null);
	}

	/** Constructor with custom logger */
	public Magnum(Logger customLog) {
		this(null, customLog);
	}

	
	/** Constructor, parse command-line arguments, initialize settings */
	public Magnum(String[] args, Logger customLog) {

		// Initialize
		if (customLog != null)
			log = customLog;
		else
			log = new Logger(); // must be first
		utils = new Utils(log);
		set = new MagnumOptionParser(this); // sets defaults
		
		// Parse command-line arguments and initialize settings
		if (args != null)
			set.parse(args);
		
		// Set verbose flag in logger
		log.setVerbose(set.verbose_);
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
		else if (set.mode_ == 5)
			new Experiments(this).run();
		else {
			set.displayHelp();
			throw new IllegalArgumentException("--mode <int> must be between 1 and 3, found mode=" + set.mode_);
		}

		log.println("Success!\n"
				  + "--------\n");
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
			ArrayList<String> networkFiles = utils.listFiles(set.networkDir_);
			// Run for each file
			ArrayList<ArrayList<Double>> networkMeans = new ArrayList<ArrayList<Double>>();
			for (String file_i : networkFiles)
				networkMeans.add(runNetworkAnalysis(new File(file_i)));
			
			// Print means to file
			FileExport writer = new FileExport(log, "networkMeans.txt");
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
		log.println("COMPUTING NETWORK PROPERTIES");
		log.println("----------------------------\n");

		NetpropMain netprop = new NetpropMain(this, network);
		ArrayList<Double> networkMeans = netprop.runAll();
		log.println();
		netprop.saveAll();
		
		return networkMeans;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Perform network operations */
	public void runNetworkOperations() {

		// Perform operations on networks
		log.println("PERFORMING NETWORK OPERATIONS");
		log.println("-----------------------------\n");

		if (set.computeUnion_) {
			GroupNetworks grouper = new GroupNetworks(this, set.networkDir_, set.networkGroupFile_, set.networkFilePrefix_);
			grouper.run();

		} else if (set.computePairwiseSum_) {
			PairwiseSum networkSum = new PairwiseSum(this, set.networkDir_, set.networkDir2_);
			networkSum.run(true);
		}
	}

	
	// ----------------------------------------------------------------------------

	/** Analyze enrichment */
	public void runEnrichmentAnalysis() {

		// Load input files
		log.println("LOADING INPUT FILES");
		log.println("-------------------\n");

		enrichMain = new EnrichMain(this);

		// Analyze enrichment
		log.println("COMPUTING ENRICHMENT");
		log.println("--------------------\n");

		enrichMain.run();
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

	
	// ----------------------------------------------------------------------------

//	/** Convert dense Colt to Apache Commons matrix */
//	public RealMatrix colt2apache(DoubleMatrix2D colt) {
//		
//		RealMatrix apache = MatrixUtils.createRealMatrix(colt.rows(), colt.columns());
//		for (int i=0; i<colt.rows(); i++)
//			for (int j=0; j<colt.columns(); j++)
//				apache.setEntry(i, j, colt.get(i, j));
//		
//		return apache;
//	}
//
//	
//	// ----------------------------------------------------------------------------
//
//	/** Convert dense Apache Commons to Colt matrix */
//	public DoubleMatrix2D apache2colt(RealMatrix apache) {
//		
//		// Convert apache commons matrix to colt
//		DoubleMatrix2D colt = new DenseDoubleMatrix2D(apache.getRowDimension(), apache.getColumnDimension());
//		for (int i=0; i<apache.getRowDimension(); i++)
//			for (int j=0; j<apache.getColumnDimension(); j++)
//				colt.set(i, j, apache.getEntry(i, j));
//		
//		return colt;
//	}

		
	// ============================================================================
	// PRIVATE METHODS

	/** Load the input network */
	private Network loadInputNetwork(File networkFile) {

		log.println("LOADING INPUT NETWORK");
		log.println("---------------------\n");

		if (set.networkDir_ != null)
			networkFile = new File(set.networkDir_, networkFile.getPath());
		
		return new Network(this, networkFile, set.refNodesFile_,
				set.isDirected_, set.removeSelfLoops_,
				set.isWeighted_, set.threshold_);
	}

	
	// ============================================================================
	// GETTERS AND SETTERS
	
	public EnrichMain getEnrichMain() { return enrichMain; }

}
