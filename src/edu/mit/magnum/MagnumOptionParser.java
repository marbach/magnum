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

import java.util.ArrayList;

import joptsimple.OptionParser;
import joptsimple.OptionSet;

/**
 * Defines and parses the command-line arguments
 */
public class MagnumOptionParser extends MagnumSettings {

	/** The option parser */
	OptionParser parser_ = null;

	// ============================================================================
	// PUBLIC METHODS

	/** Parse the command-line arguments, initialize all settings */
	public MagnumOptionParser() {

		// Defines the arguments
		defineArgs();
	}

	// ----------------------------------------------------------------------------

	/** Load settings file and parse command-line arguments (defined by defineArgs()) */
	public void parse(String[] args) {

		// (1) Parse the options
		OptionSet options = null;
		try {
			options = parser_.parse(args);
		} catch (Exception e) {
			displayHelp();
			throw new RuntimeException(e);
		}

		// Display help
		if (options.has("help") || options.has("h") || options.has("?")) {
			displayHelp();
			System.exit(0);
		}
		// verbose
		if (options.has("verbose"))
			verbose_ = true;

		// (2) Set the settings file
		if (options.has("set"))
			settingsFile = (String) options.valueOf("set");

		// (3) Load the settings file
		Magnum.log.printlnVerbose("SETTINGS FILE");
		Magnum.log.printlnVerbose("-------------\n");
		loadSettings();

		// (4) Set command-line options (default is what has been loaded from settings file)
		
		// GENERAL OPTIONS
		if (options.has("mode"))
			mode_ = (Integer) options.valueOf("mode");
		if (options.has("seed"))
			randomSeed_ = (Integer) options.valueOf("seed");
		if (options.has("verbose"))
			verbose_ = true;
		if (options.has("outdir"))
			outputDirectory_ = (String) options.valueOf("outdir");
		if (options.has("netdir"))
			networkDir_ = (String) options.valueOf("netdir");				
		if (options.has("net"))
			networkFile_ = (String) options.valueOf("net");
		if (options.has("directed"))
			isDirected_ = true;
		if (options.has("weighted"))
			isWeighted_ = true;
		if (options.has("noself"))
			removeSelfLoops_ = true;
		if (options.has("cutoff"))
			threshold_ = (Double) options.valueOf("cutoff");
		
		// MODE 1
		if (options.has("pstep"))
			computePstepKernel_ = true;
		if (options.has("nsteps")) {
			int nsteps = (Integer) options.valueOf("nsteps");
			pstepKernelP_ = new ArrayList<Integer>();
			pstepKernelP_.add(nsteps);
		}
		if (options.has("degree"))
			computeDegree_ = true;
		if (options.has("betweenness"))
			computeBetweenness_ = true;
		if (options.has("clustcoeff"))
			computeClusteringCoefficient_ = true;
		if (options.has("shortestpath"))
			computeShortestPathLengths_ = true;

		// MODE 2
		if (options.has("union"))
			computeUnion_ = true;
		
		// MODE 3
		if (options.has("genes"))
			geneCoordFile_ = (String) options.valueOf("genes");				
		if (options.has("scores"))
			geneScoreFile_ = (String) options.valueOf("scores");				
		if (options.has("cmatrix"))
			functionalDataFile_ = (String) options.valueOf("cmatrix");				
		if (options.has("excl"))
			excludedGenesFile_ = (String) options.valueOf("excl");				
		if (options.has("neighbors"))
			excludedGenesDistance_ = (Double) options.valueOf("neighbors");
		if (options.has("bins"))
			numBins_ = (Integer) options.valueOf("bins");
		if (options.has("permut"))
			numPermutations_ = (Integer) options.valueOf("permut");
		if (options.has("curve"))
			curveCutoff_ = (Double) options.valueOf("curve");

		// PRIVATE OPTIONS (not included in documentation)		
		if (options.has("scalekernel"))
			scaleKernel_ = true;
		if (options.has("idtypefdata"))
			idTypeFunctionalData_ = (String) options.valueOf("idtypefdata");
		if (options.has("chr"))
			chromosome_ = "chr" + options.valueOf("chr");
		if (options.has("gsstart"))
			geneScoreIndexStart_ = (Integer) options.valueOf("gsstart");
		if (options.has("gsend"))
			geneScoreIndexEnd_ = (Integer) options.valueOf("gsend");
		
		// TBD, write a method that checks consistency / if everything has been
		// defined that we need
		checkOptions();
	}

	// ============================================================================
	// PRIVATE METHODS

	/** Display help on console */
	static public void displayHelp() {

		Magnum.log.println("Running magnum " + MagnumSettings.magnumVersion);
		Magnum.log.println();

		Magnum.log.println("1. USAGE");
		Magnum.log.println("--------");
		Magnum.log.println("   java [JAVA OPTIONS] -jar magnum.jar --mode <int> [OPTIONS]\n");
		
		Magnum.log.println("2. JAVA OPTIONS (see java documentation for details)");
		Magnum.log.println("---------------");
		Magnum.log.println("   -Xmx<memory>    Increase memory, e.g., -Xmx8g for 8GB memory");
		Magnum.log.println("                   (more may be necessary depending on network size)");
		Magnum.log.println("   -ea             Enable assertions (not necessary, activates");
		Magnum.log.println("                   'debugging' tests built into the code)");
		Magnum.log.println();

		Magnum.log.println("3. GENERAL OPTIONS");
		Magnum.log.println("------------------");
		Magnum.log.println("   --help | -h     Display help");
		Magnum.log.println("   --verbose       Use verbose output (useful for debugging)");
		Magnum.log.println("   --mode <int>    Select the mode (REQUIRED):");
		Magnum.log.println("                      1 = Compute network properties (diffusion kernels,");
		Magnum.log.println("                          shortest paths, clustering coefficients)");
		Magnum.log.println("                      2 = Perform network operations (union)");
		Magnum.log.println("                      3 = Connectivity enrichment analysis");
		Magnum.log.println("   --seed <int>    Random number generator seed (default: 42; current time: -1)");
		Magnum.log.println("   --outdir <dir>  Output directory (default: working directory)");
		Magnum.log.println("   --netdir <dir>  Directory of input networks (default: working directory)");
		Magnum.log.println("   --net <file>    Input network filename");
		Magnum.log.println("   --directed      Input network is directed (default: undirected)");
		Magnum.log.println("   --weighted      Input network is weighted (default: unweighted)");
		Magnum.log.println("   --noself        Remove self-loops from input network");
		Magnum.log.println("   --cutoff <value> Remove edges with weight < cutoff from input network");
		Magnum.log.println();

		Magnum.log.println("4. NETWORK PROPERTIES");
		Magnum.log.println("---------------------");
		Magnum.log.println("   --pstep         P-step random walk kernel (Smola & Kondor, 2003; allows for");
		Magnum.log.println("                   weighted networks)");
		Magnum.log.println("   --nsteps <int>  Number of steps for p-step random walk kernel (default: 4)");
		Magnum.log.println("   --degree        Node degree (directed networks, also indegree and outdegree)");
		Magnum.log.println("   --betweenness   Node betweenness centrality (allows for directed networks)");
		Magnum.log.println("   --clustcoeff    Node clustering coefficient (allows for directed networks)");
		Magnum.log.println("   --shortestpath  Shortest path lengths and closeness centrality");
		Magnum.log.println();
		
		Magnum.log.println("5. NETWORK OPERATIONS");
		Magnum.log.println("---------------------");
		Magnum.log.println("   --union         Union (max edge weight) of all networks in network directory");
		Magnum.log.println("                   (see option: --netdir <dir>)");
		Magnum.log.println();

		Magnum.log.println("6. NETWORK CONNECTIVITY ENRICHMENT");
		Magnum.log.println("----------------------------------");
		Magnum.log.println("   --genes <file>  The gene coordinates (REQUIRED)");
		Magnum.log.println("   --scores <file> The GWAS gene scores (REQUIRED)");
		Magnum.log.println("   --cmatrix <file> The connectivity matrix (e.g., diffusion kernel; REQUIRED)");
		Magnum.log.println("   --excl <file>   Genes to be excluded (e.g., HLA region)");
		Magnum.log.println("   --neighbors <X> Ignore connectivity between genes with distance < X mega-bases");
		Magnum.log.println("                   (default: 1 [mega-base])");
		Magnum.log.println("   --bins <int>    The number of bins for within-degree permutation (default: 100)");
		Magnum.log.println("   --permut <int>  No. permutations to compute empirical p-values (default: 10000)");
		Magnum.log.println("   --curve <X>     Compute curves only for the top part of the ranked gene list");
		Magnum.log.println("                   (default: 0.2 [top 20%])");
		Magnum.log.println();

		Magnum.log.println("7. EXAMPLES");
		Magnum.log.println("-----------");
		Magnum.log.println("See the step-by-step tutorial in the user guide for details.");
		Magnum.log.println();		
		Magnum.log.println("Compute random walk kernel for the network included in the tutorial directory:");
		Magnum.log.println("   >> java -Xmx6g -ea -jar magnum_v1.0.jar --mode 1 --pstep --netdir tutorial_data --net smooth_muscle_cells_-_umbilical_vein.txt.gz --weighted");
		Magnum.log.println();		
		Magnum.log.println("Perform network connectivity enrichment analysis:");
		Magnum.log.println("   >> java -Xmx6g -ea -jar magnum_v1.0.jar --mode 3 --genes tutorial_data/gene_coord.bed --excl tutorial_data/excluded_genes.txt --scores tutorial_data/macular_degeneration_neovascular.txt --cmatrix smooth_muscle_cells_-_umbilical_vein_4stepKernel_alpha2.0_weighted.txt.gz --permut 10000");
		Magnum.log.println();		
	}

	
	// ----------------------------------------------------------------------------

	/** TODO Check validity of arguments 
	 * TODO check default values */
	private void checkOptions() {
		
		if (excludedGenesDistance_ >= 1000)
			throw new IllegalArgumentException("excludedGenesDistance is given in mega bases and cannot exceed a value of 1000mb");
	}

		
		
	// ----------------------------------------------------------------------------

	/** Defines the command-line arguments */
	private void defineArgs() {

		// The command-line parser
		parser_ = new OptionParser();

		// Display help
		parser_.accepts("help");
		parser_.accepts("h");
		parser_.accepts("?");
		// settingsFile_
		parser_.accepts("set").withRequiredArg();
		
		// verbose_
		parser_.accepts("verbose");
		// mode_
		parser_.accepts("mode").withRequiredArg().ofType(Integer.class);
		// randomSeed_
		parser_.accepts("seed").withRequiredArg().ofType(Integer.class);
		// outputDirectory_
		parser_.accepts("outdir").withRequiredArg();
		
		// networkDir_
		parser_.accepts("netdir").withRequiredArg();
		// networkFile_
		parser_.accepts("net").withRequiredArg();
		// isDirected_
		parser_.accepts("directed");
		// isWeighted_
		parser_.accepts("weighted");
		// removeSelfLoops_
		parser_.accepts("noself");
		// threshold_
		parser_.accepts("cutoff").withRequiredArg().ofType(Double.class);
		
		// kernels
		parser_.accepts("pstep");
		parser_.accepts("nsteps").withRequiredArg().ofType(Integer.class);
		// network properties
		parser_.accepts("degree");
		parser_.accepts("betweenness");
		parser_.accepts("clustcoeff");
		parser_.accepts("shortestpath");
		
		// network operations
		parser_.accepts("union");
		
		// Enrichment / gene scores
		parser_.accepts("genes").withRequiredArg();
		parser_.accepts("scores").withRequiredArg();
		parser_.accepts("cmatrix").withRequiredArg();
		parser_.accepts("excl").withRequiredArg();
		parser_.accepts("neighbors").withRequiredArg().ofType(Double.class);
		parser_.accepts("bins").withRequiredArg().ofType(Integer.class);
		parser_.accepts("permut").withRequiredArg().ofType(Integer.class);
		parser_.accepts("curve").withRequiredArg().ofType(Double.class);
		
		// idTypeFunctionalData_
		parser_.accepts("idtypefdata").withRequiredArg();
		// chromosome_
		parser_.accepts("chr").withRequiredArg().ofType(Integer.class);
		// scaleKernel_
		parser_.accepts("scalekernel");
		// geneScoreIndexStart_
		parser_.accepts("gsstart").withRequiredArg().ofType(Integer.class);
		// geneScoreIndexEnd_
		parser_.accepts("gsend").withRequiredArg().ofType(Integer.class);

		// Example
		// parser_.accepts("cut").withRequiredArg().ofType(Integer.class).defaultsTo(-1);
	}

}
