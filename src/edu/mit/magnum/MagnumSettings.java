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
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Properties;
import java.util.Random;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.Well19937c;



/** 
 * Offers global parameters (settings) and functions used by all classes of the
 * ngsea package.
 */
public class MagnumSettings extends Settings {	
	
	/** The configuration file with the settings (set null for default settings) */
	public String settingsFile = null;
	/** Flag indicates whether the parser should expect the settings file to be complete */
	public boolean requireAllOptions = false;

	/** Colt Mersenne Twister random engine (should be used by all other random number generators) */
	public MersenneTwister mersenneTwisterRng_ = null;
	/** Apache Commons random engine */
	public Well19937c wellRng_ = null;
	/** Java random engine */
	public Random jdkRng_ = null;

	// ----------------------------------------------------------------------------
	// VARIOUS
	
	/** Mode: 1 => Network analysis; 2 => Enrichment analysis */
	public int mode_ = 0;
	/** PRIVATE, NEEDS TO BE SET WITH setRandomSeed(), which initializes the random number generators. Set to -1 to use current time */
	private int randomSeed_ = 42;
	/** Set true to use verbose mode (print more information) */
	public boolean verbose_ = false;
	/** Output directory to save stuff */
	public String outputDirectory_ = ".";
	/** Output filename */
	public String outputFilename_ = "";
	/** Compress output files (gzip) */
	public boolean compressFiles_ = true;

	// ----------------------------------------------------------------------------
	// NETWORK PROPERTIES

	// INPUT NETWORK
	/** Directory containing the networks */
	public String networkDir_ = null;
	/** The input network file */
	public String networkFile_ = null;
	/** Delimiter used to separate columns (default 'tab' */
	public String networkFileDelim_ = "TAB";  
	/** Defines if the network should be interpreted as directed or undirected */
	public boolean isDirected_ = true;
	/** Defines if self loops should be removed from the network */
	public boolean removeSelfLoops_ = true;
	/** Set true to treat the network as weighted */
	public boolean isWeighted_ = true;
	/** Threshold for including weighted edges */
	public double threshold_ = 0;
	/** Exclude "super-hubs" that connect to more than the given fraction of genes (set 1 to include all) */
	public double superHubThreshold_ = 0;
	/** Optional file specifying a set of reference nodes */
	public String refNodesFile_ = null;

	// NETWORKOPS
	/** Take union (max edge) over all networks in networkDir or the sets specified in the file below */
	public boolean computeUnion_ = false;
	/** Define the network sets that should be combined (leave empty to combine all networks) */
	public String networkGroupFile_ = null;
	/** Prefix of the files in the network dir */
	public String networkFilePrefix_ = null;
	/** Add networks of the same cell type */
	public boolean computePairwiseSum_ = false;
	/** The network directory of the second networks */
	public String networkDir2_ = null;
	
	// BASIC NETWORK PROPERTIES
	/** Node degree (directed networks, also indegree and outdegree) */
	public boolean computeDegree_ = false;
	/** Node betweenness centrality (edge directionality observed for directed networks) */
	public boolean computeBetweenness_ = false;
	/** Node clustering coefficient (edge directionality observed for directed networks) */
	public boolean computeClusteringCoefficient_ = false;
	/** For each node, distance to all other nodes (or all reference nodes) and closeness centrality */
	public boolean computeShortestPathLengths_ = false;

	// KERNELS
	/** P-step random walk kernel (Smola & Kondor, 2003) */
	public boolean computePstepKernel_ = false;
	/** alpha parameter of p-step random walk kernel (alpha >= 2) */
	public double pstepKernelAlpha_ = 2;
	/** Number of steps p of random walk kernel (p >= 1) */
	public ArrayList<Integer> pstepKernelP_ = new ArrayList<Integer>(Arrays.asList(4));
	/** Normalize the kernel matrix (divide by the max) */
	public boolean pstepKernelNormalize_ = true;
	
	// TANIMOTO COEFFICIENT
	/** Tanimoto coefficient between target genes */
	public boolean computeTargetTanimoto_ = false;
	/** Tanimoto coefficient between regulators */
	public boolean computeTfTanimoto_ = false;

	// OUTPUT FILES
	/** A suffix/ending that is appended to all output files for this run (use to distinguish output files from multiple runs) */
	public String outputSuffix_ = "";
	/** Export all computed pairwise node properties (e.g., similarity, distance matrices) */
	public boolean exportPairwiseNodeProperties_ = true;
	/** Export all computed node properties (e.g., avg. similarity, distance for each node) */
	public boolean exportNodeProperties_ = true;

	// ----------------------------------------------------------------------------
	// GENOME ANNOTATION

	/** Set of genes to be considered (leave empty to use all genes from the annotation file) */
	public String genesToBeLoadedFile_ = null;

	/** The chromosome to be considered (chr1, ..., chr22, chrX, chrY), leave empty for all chromosomes */ 
	public String chromosome_ = null;
	/** Ignore sex chromosomes */
	public boolean ignoreAllosomes_ = true;

	/** The file with the gencode annotation */
	public String gencodeAnnotationFile_ = null;
	/** UCSC genome browser annotation (use for Entrez IDs) */ 
	public String ucscAnnotationFile_ = null;
	/** Set true to load only protein-coding genes */
	public boolean loadOnlyProteinCodingGenes_ = true;

	/** Mapping file to convert Entrez IDs, ENSEMBL IDs and gene symbols */
	public String geneIdMappingFile_ = null;
		
	// ----------------------------------------------------------------------------
	// ENRICHMENT CURVES
	
	// INPUT
	/** The gene coordinates (custom annotation) */
	public String geneCoordFile_ = null;
	/** The gene scores */
	public String geneScoreFile_ = null;
	/** Cutoff for genome-wide significance of gene scores */
	public double genomeWideSignificanceThreshold_ = 1e-6;
	/** Exclude genome-wide significant genes (below threshold) */
	public boolean excludeGenomeWideSignificantGenes_ = false;

	/** The file with the functional data, e.g. network kernels (cols: gene id, property 1, property 2, ...) */ 
	public String functionalDataFile_ = null; 
	/** Specify which columns should be loaded (-1: all columns; 1: first gene property column) */
	public ArrayList<Integer> functionalDataCols_ = null;

	/** Genes to be excluded from enrichment analysis (e.g., MHC region) */
	public String excludedGenesFile_ = null;
	/** Gene pairs to be excluded from enrichment analysis (e.g., genes in LD) */
	public String excludedGenePairsFile_ = null;
	/** Exclude gene pairs with windows smaller than the given distance apart (given in megabases; -1: no exclusion; 1000000: all genes on same chromosome) */
	public double excludedGenesDistance_ = 1; 

	/** Gene IDs used in geneScoreFile, excludedGenesFile, excludedGenePairsFile ('ensembl', 'entrez', 'hugo') */
	public String idTypeGeneScores_ = "custom";
	/** Gene IDs used in functionalDataFile */
	public String idTypeFunctionalData_ = "custom";

	// PARAMETERS
	/** Number of random permutations to compute confidence intervals */
	public int numPermutations_ = 10000;
	/** The number of bins for within-degree permutation */
	public int numBins_ = 100;
	/** Scale kernels: K'(i,j) = K(i,j)/sqrt(rowSums(K)[i] * colSums(K)[j]) */
	public boolean scaleKernel_ = false;

	/** Equidistant curve resolution, e.g., set 10 to compute every 10th point on the curves */
	public int constCurveResolution_ = 10;
	/** Varying curve resolution, e.g., set 2 to compute points: 2, 6, 12, 20, 30, 42, ... (takes precedence over constCurveResolution, set -1 to disable) */
	public int varCurveResolution_ = -1;
	/** Compute curves only for the top part of the list (e.g., 0.1 for top 10%) */
	public double curveCutoff_ = 0.2;
	/** Sliding window size */
	public int slidingWindowSize_ = -1;

	/** Draw boundaries for given p-values (e.g., set 0.05 to draw the upper/lower boundary where only 5% of random curves above/below) */ 
	public ArrayList<Double> pval_ = new ArrayList<Double>(Arrays.asList(0.01, 0.05));
	/** Indicates whether the boundaries are one-sided or two-sided (equivalent to dividing pval by two) */
	public boolean twoSidedTest_ = false;
	/** Control FDR over all points of the curve together (i.e., correct for multiple hypothesis testing across curve) */
	public boolean controlFDR_ = false;
	/** The top X snps will be ignored for FDR control (too noisy at start of the list) */
	public int FDRStart_ = 100;
	/** Start to integrate the AUC only at k=10 (reduce noise at the start of the list) */
	public int AUCStart_ = 10;
	/** Index of gene scores for which enrichment is to be computed, e.g., (0,9) for the first ten gene scores */ 
	public int geneScoreIndexStart_ = 0;
	public int geneScoreIndexEnd_ = 0;

	
	// OUTPUT FILES
	/** Number of random permutations for which enrichment curves are exported (smaller or equal numPermutations) */
	public int numPermutationsExport_ = 0;

	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public MagnumSettings() {
		resetToDefaults();
	}
	
	
	// ----------------------------------------------------------------------------
	
	/** Set default values for all settings */
	public void resetToDefaults() {
		
		settingsFile = null;
		requireAllOptions = false;

		// Initializes the RNGs
		setRandomSeed(42);

		mode_ = 0;
		verbose_ = false;
		outputDirectory_ = ".";
		outputFilename_ = "";
		compressFiles_ = true;

		networkDir_ = null;
		networkFile_ = null;
		networkFileDelim_ = "TAB";  
		isDirected_ = true;
		removeSelfLoops_ = true;
		isWeighted_ = true;
		threshold_ = 0;
		superHubThreshold_ = 0;
		refNodesFile_ = null;

		computeUnion_ = false;
		networkGroupFile_ = null;
		networkFilePrefix_ = "";
		computePairwiseSum_ = false;
		networkDir2_ = null;
		
		computeDegree_ = false;
		computeBetweenness_ = false;
		computeClusteringCoefficient_ = false;
		computeShortestPathLengths_ = false;

		computePstepKernel_ = false;
		pstepKernelAlpha_ = 2;
		pstepKernelP_ = new ArrayList<Integer>();
		pstepKernelP_.add(4);
		pstepKernelNormalize_ = true;
		
		computeTargetTanimoto_ = false;
		computeTfTanimoto_ = false;

		outputSuffix_ = "";
		exportPairwiseNodeProperties_ = true;
		exportNodeProperties_ = true;

		genesToBeLoadedFile_ = null;

		chromosome_ = null;
		ignoreAllosomes_ = true;

		gencodeAnnotationFile_ = null;
		ucscAnnotationFile_ = null;
		loadOnlyProteinCodingGenes_ = true;

		geneIdMappingFile_ = null;
			
		geneCoordFile_ = null;
		geneScoreFile_ = null;
		genomeWideSignificanceThreshold_ = 1e-6;
		excludeGenomeWideSignificantGenes_ = false;

		functionalDataFile_ = null; 
		functionalDataCols_ = null;

		excludedGenesFile_ = null;
		excludedGenePairsFile_ = null;
		excludedGenesDistance_ = 1; 

		idTypeGeneScores_ = "custom";
		idTypeFunctionalData_ = "custom";

		numPermutations_ = 10000;
		numBins_ = 100;
		scaleKernel_ = false;

		constCurveResolution_ = 10;
		varCurveResolution_ = -1;
		curveCutoff_ = 0.2;
		slidingWindowSize_ = -1;

		pval_ = new ArrayList<Double>();
		pval_.add(0.01);
		pval_.add(0.05);
		twoSidedTest_ = false;
		controlFDR_ = false;
		FDRStart_ = 100;
		AUCStart_ = 10;
		geneScoreIndexStart_ = 0;
		geneScoreIndexEnd_ = 0;
		
		numPermutationsExport_ = 0;		
	}

	
	// ----------------------------------------------------------------------------
	
	/** Load settings from given file */
	public void loadSettings(String settingsFile) {
		
		Magnum.log.printlnVerbose("SETTINGS FILE");
		Magnum.log.printlnVerbose("-------------\n");
		
		try {
			// Open the input stream
			//InputStream in = MagnumSettings.class.getClassLoader().getResourceAsStream("edu/mit/magnum/settings.txt");
				
			// Check that the specified settings file exists
			if (settingsFile == null || settingsFile.isEmpty())
				throw new RuntimeException("No settings file specified");
			else if (!new File(settingsFile).exists())
				throw new RuntimeException("Settings file not found: " + settingsFile);

			// Open file input stream
			Magnum.log.println("- Loading settings file: " + settingsFile + "\n");
			InputStream in = new FileInputStream(settingsFile);

			// Load the settings
			prop = new Properties();
			prop.load(new InputStreamReader(in));
			
			// Get the param values
			setParameterValues();
			
		} catch (Exception e) {
			Magnum.log.warning(e.getMessage());
			Magnum.log.error("Failed to load settings file (a parameter may be missing or malformed): " + settingsFile);
		}		
	}
	
	
	// ----------------------------------------------------------------------------
	
	/**
	 * Set the user path. Could be with or without "/" terminal.
	 * @param absPath Absolute path
	 */
	public void setOutputDirectory(String absPath) {
		
		outputDirectory_ = absPath;
		String sep = System.getProperty("file.separator");
		if (outputDirectory_.charAt(outputDirectory_.length()-1) != sep.charAt(0))
			outputDirectory_ += sep;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Create new instances for the random number generators, initialize with randomSeed_ */
	public void setRandomSeed(int seed) {
		
		randomSeed_ = seed;
		if (randomSeed_ == -1) {
			mersenneTwisterRng_ = new MersenneTwister();
			wellRng_ = new Well19937c();
			jdkRng_ = new Random();
		} else {
			mersenneTwisterRng_ = new MersenneTwister(randomSeed_);
			wellRng_ = new Well19937c(randomSeed_);
			jdkRng_ = new Random(randomSeed_);
		}
		
		//uniformDistribution_ = new Uniform(mersenneTwister_);
		//normalDistribution_ = new Normal(0, 1, mersenneTwister_); // mean 0, stdev 1
	}
	

	// ============================================================================
	// PRIVATE METHODS

	/** Set ngsea parameters based on the loaded properties */
	private void setParameterValues() throws Exception {

		// VARIOUS
		mode_ = getSettingInt("mode");
		setRandomSeed(getSettingInt("randomSeed"));
		verbose_ = getSettingBoolean("verbose");
		outputDirectory_ = getSetting("outputDirectory");
		if (outputDirectory_.equals("")) 
			outputDirectory_ = System.getProperty("user.dir");
		outputFilename_ = getSetting("outputFilename");
		
		// INPUT NETWORK
		networkDir_ = getSetting("networkDir");
		networkFile_ = getSetting("networkFile");
		networkFileDelim_ = getSetting("networkFileDelim");
		isDirected_ = getSettingBoolean("isDirected");
		removeSelfLoops_ = getSettingBoolean("removeSelfLoops");
		isWeighted_ = getSettingBoolean("isWeighted");
		threshold_ = getSettingDouble("threshold");
		superHubThreshold_ = getSettingDouble("superHubThreshold");
		refNodesFile_ = getSetting("refNodesFile");

		// OUTPUT FILES
		outputSuffix_ = getSetting("outputSuffix");
		exportPairwiseNodeProperties_ = getSettingBoolean("exportPairwiseNodeProperties");
		exportNodeProperties_ = getSettingBoolean("exportNodeProperties");
		compressFiles_ = getSettingBoolean("compressFiles");
		
		// NETWORKOPS
		computeUnion_ = getSettingBoolean("computeUnion");
		networkGroupFile_ = getSetting("networkGroupFile");
		networkFilePrefix_ = getSetting("networkFilePrefix");
		computePairwiseSum_ = getSettingBoolean("computePairwiseSum");
		networkDir2_ = getSetting("networkDir2");
		
		// BASIC NETWORK PROPERTIES
		computeDegree_ = getSettingBoolean("computeDegree");
		computeBetweenness_ = getSettingBoolean("computeBetweenness");
		computeClusteringCoefficient_ = getSettingBoolean("computeClusteringCoefficient");
		
		// SHORTEST PATHS
		computeShortestPathLengths_ = getSettingBoolean("computeShortestPathLengths");

		// KERNELS
		computePstepKernel_ = getSettingBoolean("computePstepKernel");
		pstepKernelAlpha_ = getSettingDouble("pstepKernelAlpha");
		pstepKernelP_ = getSettingIntArray("pstepKernelP", true);
		pstepKernelNormalize_ = getSettingBoolean("pstepKernelNormalize");
		
		// TANIMOTO
		computeTargetTanimoto_ = getSettingBoolean("computeTargetTanimoto");
		computeTfTanimoto_ = getSettingBoolean("computeTfTanimoto");
		
		// ----------------------------------------------------------------------------
		// GENOME ANNOTATION
		
		genesToBeLoadedFile_ = getSetting("genesToBeLoadedFile");
		chromosome_ = getSetting("chromosome");
		ignoreAllosomes_ = getSettingBoolean("ignoreAllosomes");
		gencodeAnnotationFile_ = getSetting("genecodeAnnotationFile");
		ucscAnnotationFile_ = getSetting("ucscAnnotationFile");
		loadOnlyProteinCodingGenes_ = getSettingBoolean("loadOnlyProteinCodingGenes");
		geneIdMappingFile_ = getSetting("geneIdMappingFile");
		
		// ----------------------------------------------------------------------------
		// ENRICHMENT ANALYSIS

		geneCoordFile_ = getSetting("geneCoordFile");
		if (!geneCoordFile_.equals("")) {
			idTypeFunctionalData_ = "custom";
			idTypeGeneScores_ = "custom";
		}

		geneScoreFile_ = getSetting("geneScoreFile");
		genomeWideSignificanceThreshold_ = getSettingDouble("genomeWideSignificanceThreshold");
		excludeGenomeWideSignificantGenes_ = getSettingBoolean("excludeGenomeWideSignificantGenes");

		functionalDataFile_ = getSetting("functionalDataFile"); 
		functionalDataCols_ = getSettingIntArray("functionalDataCols", true);		
		
		excludedGenesFile_ = getSetting("excludedGenesFile");
		excludedGenePairsFile_ = getSetting("excludedGenePairsFile");
		excludedGenesDistance_ = getSettingInt("excludedGenesDistance");
		
		idTypeGeneScores_ = getSetting("idTypeGeneScores");
		idTypeFunctionalData_ = getSetting("idTypeFunctionalData");
		
		// ENRICHMENT
		numPermutations_ = getSettingInt("numPermutations");
		numBins_ = getSettingInt("numBins");
		scaleKernel_ = getSettingBoolean("scaleKernel");

		constCurveResolution_ = getSettingInt("constCurveResolution");
		varCurveResolution_ = getSettingInt("varCurveResolution");
		curveCutoff_ = getSettingDouble("curveCutoff");
		slidingWindowSize_ = getSettingInt("slidingWindowSize");
		
		numPermutationsExport_ = getSettingInt("numPermutationsExport");
		if (numPermutationsExport_ > numPermutations_)
			throw new IllegalArgumentException("Invalid settings: 'numPermutationsExport' must be smaller or equal 'numPermutations'");
		pval_ = getSettingDoubleArray("pval", true);
		twoSidedTest_ = getSettingBoolean("twoSidedTest");
		controlFDR_ = getSettingBoolean("controlFDR");
		FDRStart_ = getSettingInt("FDRStart");
		AUCStart_ = getSettingInt("AUCStart");
		geneScoreIndexStart_ = getSettingInt("geneScoreIndexStart");
		geneScoreIndexEnd_ = getSettingInt("geneScoreIndexEnd");
		
	}
	
	
	// ============================================================================
	// GETTERS AND SETTERS
	
	public int getRandomSeed() { return randomSeed_; }

	
}
