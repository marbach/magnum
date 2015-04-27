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
import java.util.Properties;
import java.util.Random;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.Well19937c;



/** 
 * Offers global parameters (settings) and functions used by all classes of the
 * ngsea package.
 * 
 * Code adapted from GnwSettings.java by Thomas Schaffter and Daniel Marbach (gnw.sf.net)
 */
public class Settings {	
	
	/** The configuration file with the settings (leave empty for default settings) */
	static public String settingsFile_ = null;
	/** The properties (settings file) */
	static private Properties set_ = null;
	
	/** Colt Mersenne Twister random engine (should be used by all other random number generators) */
	static public MersenneTwister mersenneTwisterRng_ = null;
	/** Apache Commons random engine */
	static public Well19937c wellRng_ = null;
	/** Java random engine */
	static public Random jdkRng_ = null;

	// ----------------------------------------------------------------------------
	// VARIOUS
	
	/** Current version of N-GSEA */
	static public String magnumVersion_ = "1.0 alpha";
	/** Mode: 1 => Network analysis; 2 => Enrichment analysis */
	static public int mode_ = -1;
	/** Seed for the random number generator. Set to -1 to use current time */
	static public int randomSeed_ = -1;
	/** Set true to use verbose mode (print more information) */
	static public boolean verbose_ = false;
	/** Output directory to save stuff */
	static public String outputDirectory_ = "";
	/** Output filename */
	static public String outputFilename_ = "";
	/** Compress output files (gzip) */
	static public boolean compressFiles_ = true;

	// ----------------------------------------------------------------------------
	// NETWORK PROPERTIES

	// INPUT NETWORK
	/** Directory containing the networks */
	static public String networkDir_ = null;
	/** The input network file */
	static public String networkFile_ = null;
	/** Delimiter used to separate columns (default 'tab' */
	static public String networkFileDelim_ = null;  
	/** Defines if the network should be interpreted as directed or undirected */
	static public boolean isDirected_ = true;
	/** Defines if self loops should be removed from the network */
	static public boolean removeSelfLoops_ = false;
	/** Set true to treat the network as weighted */
	static public boolean isWeighted_ = false;
	/** Threshold for including weighted edges */
	static public double threshold_ = 0;
	/** Exclude "super-hubs" that connect to more than the given fraction of genes (set 1 to include all) */
	static public double superHubThreshold_ = 1;
	/** Optional file specifying a set of reference nodes */
	static public String refNodesFile_ = null;

	// NETWORKOPS
	/** Take union (max edge) over all networks in networkDir or the sets specified in the file below */
	static public boolean computeUnion_ = false;
	/** Define the network sets that should be combined (leave empty to combine all networks) */
	static public String networkGroupFile_ = null;
	/** Prefix of the files in the network dir */
	static public String networkFilePrefix_ = null;
	/** Add networks of the same cell type */
	static public boolean computePairwiseSum_ = false;
	/** The network directory of the second networks */
	static public String networkDir2_ = null;
	
	// BASIC NETWORK PROPERTIES
	/** Node degree (directed networks, also indegree and outdegree) */
	static public boolean computeDegree_ = true;
	/** Node betweenness centrality (edge directionality observed for directed networks) */
	static public boolean computeBetweenness_ = false;
	/** Node clustering coefficient (edge directionality observed for directed networks) */
	static public boolean computeClusteringCoefficient_ = false;
	/** For each node, distance to all other nodes (or all reference nodes) and closeness centrality */
	static public boolean computeShortestPathLengths_ = false;

	// KERNELS
	/** P-step random walk kernel (Smola & Kondor, 2003) */
	static public boolean computePstepKernel_ = false;
	/** alpha parameter of p-step random walk kernel (alpha >= 2) */
	static public ArrayList<Double> pstepKernelAlpha_ = null;
	/** Number of steps p of random walk kernel (p >= 1) */
	static public ArrayList<Integer> pstepKernelP_ = null;
	/** Normalize the kernel matrix (divide by the max) */
	static public boolean pstepKernelNormalize_ = true;
	
	// TANIMOTO COEFFICIENT
	/** Tanimoto coefficient between target genes */
	static public boolean computeTargetTanimoto_ = false;
	/** Tanimoto coefficient between regulators */
	static public boolean computeTfTanimoto_ = false;

	// OUTPUT FILES
	/** A suffix/ending that is appended to all output files for this run (use to distinguish output files from multiple runs) */
	static public String outputSuffix_ = "";
	/** Export all computed pairwise node properties (e.g., similarity, distance matrices) */
	static public boolean exportPairwiseNodeProperties_ = true;
	/** Export all computed node properties (e.g., avg. similarity, distance for each node) */
	static public boolean exportNodeProperties_ = true;

	// ----------------------------------------------------------------------------
	// GENOME ANNOTATION

	/** Set of genes to be considered (leave empty to use all genes from the annotation file) */
	static public String genesToBeLoadedFile_ = null;

	/** The chromosome to be considered (chr1, ..., chr22, chrX, chrY), leave empty for all chromosomes */ 
	static public String chromosome_ = null;
	/** Ignore sex chromosomes */
	static public boolean ignoreAllosomes_ = true;

	/** The file with the gencode annotation */
	static public String gencodeAnnotationFile_ = null;
	/** UCSC genome browser annotation (use for Entrez IDs) */ 
	static public String ucscAnnotationFile_ = null;
	/** Set true to load only protein-coding genes */
	static public boolean loadOnlyProteinCodingGenes_ = false;

	/** Mapping file to convert Entrez IDs, ENSEMBL IDs and gene symbols */
	static public String geneIdMappingFile_ = null;
		
	// ----------------------------------------------------------------------------
	// ENRICHMENT CURVES
	
	// INPUT
	/** The gene coordinates (custom annotation) */
	static public String geneCoordFile_ = null;
	/** The gene scores */
	static public String geneScoreFile_ = null;
	/** Cutoff for genome-wide significance of gene scores */
	static public double genomeWideSignificanceThreshold_ = -1;
	/** Exclude genome-wide significant genes (below threshold) */
	static public boolean excludeGenomeWideSignificantGenes_ = false;

	/** The file with the functional data, e.g. network kernels (cols: gene id, property 1, property 2, ...) */ 
	static public String functionalDataFile_ = null; 
	/** Specify which columns should be loaded (-1: all columns; 1: first gene property column) */
	static public ArrayList<Integer> functionalDataCols_ = null;

	/** Genes to be excluded from enrichment analysis (e.g., MHC region) */
	static public String excludedGenesFile_ = null;
	/** Gene pairs to be excluded from enrichment analysis (e.g., genes in LD) */
	static public String excludedGenePairsFile_ = null;
	/** Exclude gene pairs with windows smaller than the given distance apart (given in megabases; -1: no exclusion; 1000000: all genes on same chromosome) */
	static public double excludedGenesDistance_ = -1; 

	/** Gene IDs used in geneScoreFile, excludedGenesFile, excludedGenePairsFile ('ensembl', 'entrez', 'hugo') */
	static public String idTypeGeneScores_ = null;
	/** Gene IDs used in functionalDataFile */
	static public String idTypeFunctionalData_ = null;

	// PARAMETERS
	/** Number of random permutations to compute confidence intervals */
	static public int numPermutations_ = -1;
	/** The number of bins for within-degree permutation */
	static public int numBins_ = 1;
	/** Scale kernels: K'(i,j) = K(i,j)/sqrt(rowSums(K)[i] * colSums(K)[j]) */
	static public boolean scaleKernel_ = false;

	/** Equidistant curve resolution, e.g., set 10 to compute every 10th point on the curves */
	static public int constCurveResolution_ = 1;
	/** Varying curve resolution, e.g., set 2 to compute points: 2, 6, 12, 20, 30, 42, ... (takes precedence over constCurveResolution, set -1 to disable) */
	static public int varCurveResolution_ = -1;
	/** Compute curves only for the top part of the list (e.g., 0.1 for top 10%) */
	static public double curveCutoff_ = 1;
	/** Sliding window size */
	static public int slidingWindowSize_ = -1;

	/** Draw boundaries for given p-values (e.g., set 0.05 to draw the upper/lower boundary where only 5% of random curves above/below) */ 
	static public ArrayList<Double> pval_ = null;
	/** Indicates whether the boundaries are one-sided or two-sided (equivalent to dividing pval by two) */
	static public boolean twoSidedTest_ = false;
	/** Control FDR over all points of the curve together (i.e., correct for multiple hypothesis testing across curve) */
	static public boolean controlFDR_ = false;
	/** The top X snps will be ignored for FDR control (too noisy at start of the list) */
	static public int FDRStart_ = 0;
	/** Start to integrate the AUC only at k=10 (reduce noise at the start of the list) */
	static public int AUCStart_ = 0;
	/** Index of gene scores for which enrichment is to be computed, e.g., (0,9) for the first ten gene scores */ 
	static public int geneScoreIndexStart_ = -1;
	static public int geneScoreIndexEnd_ = -1;

	
	// OUTPUT FILES
	/** Number of random permutations for which enrichment curves are exported (smaller or equal numPermutations) */
	static public int numPermutationsExport_ = -1;
	/** Prefix for output files (usually the network name) */
	static public String outputPrefix_ = null;

	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Initialize settings */
	static public void initialize() {
		
		initializeRandomNumberGenerators();
	}
	
	// ----------------------------------------------------------------------------
	
	/** Load and initialize settings */
	static public void loadSettings() {
		
		try {
			InputStream in;
			if (settingsFile_ == null || settingsFile_.compareTo("") == 0) {
				Magnum.printlnVerbose("- No settings file specified");
				Magnum.printlnVerbose("- Loading default settings\n");
				in = Settings.class.getClassLoader().getResourceAsStream("edu/mit/magnum/settings.txt");
			} else {
				Magnum.println("- Loading settings file: " + settingsFile_ + "\n");
				in = new FileInputStream(settingsFile_);
			}
			set_ = new Properties();
			set_.load(new InputStreamReader(in));
			
			setParameterValues();
			
		} catch (Exception e) {
			Magnum.warning(e.getMessage());
			Magnum.error("Failed to load settings file (a parameter may be missing or malformed): " + settingsFile_);
		}
		
		// Reinitialize the random number generators
		initializeRandomNumberGenerators();		
	}
	
	
	// ----------------------------------------------------------------------------
	
	/** Get folder of the default ngsea directory (in home) */
	public File getGseaDirectory() {
		
		File folder = new File(ngseaDirectoryPath());
		return folder;
	}
	
	
	// ----------------------------------------------------------------------------
	
	/** Get file associated to default settings file */
	public File getCustomNgseaSettings() {
		
		File file = new File(personalNgseaSettingsPath());
		return file;
	}
	
	
	// ----------------------------------------------------------------------------
	
	/** Get path to default ngsea directory (in home) */
	public String ngseaDirectoryPath() {
		
		return System.getProperty("user.home")
				+ System.getProperty("file.separator")
				+ "ngsea";
	}
	
	
	// ----------------------------------------------------------------------------
	
	/** Get path to default settings file */
	public String personalNgseaSettingsPath() {
		
		return ngseaDirectoryPath()
				+ System.getProperty("file.separator")
				+ "settings.txt";
	}
	
	
	// ----------------------------------------------------------------------------
	
	/** Return true if a custom default ngsea settings file exists */
	public boolean personalNgseaSettingsExist() {
		
		return (new File(personalNgseaSettingsPath())).exists();
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
	static public void initializeRandomNumberGenerators() {
		
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
	static private void setParameterValues() throws Exception {

		// VARIOUS
		mode_ = getSettingInt("mode");
		randomSeed_ = getSettingInt("randomSeed");
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
//tmp
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
		pstepKernelAlpha_ = getSettingDoubleArray("pstepKernelAlpha", false);
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
			Settings.idTypeFunctionalData_ = "custom";
			Settings.idTypeGeneScores_ = "custom";
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
		
		// OUTPUT FILES
		outputPrefix_ = getSetting("outputPrefix");

	}
	
	
	// ----------------------------------------------------------------------------

	/** Get the string value of a parameter from the setting file */
	static private String getSetting(String param) {
		
		String value = set_.getProperty(param);
		if (value == null)
			Magnum.error("Parameter not found in setting file: " + param);
		
		return value; 
	}

	
	// ----------------------------------------------------------------------------

	/** Get the integer value of a parameter from the setting file */
	static private int getSettingInt(String param) {
		return Integer.valueOf(getSetting(param)); 
	}

	/** Get the double value of a parameter from the setting file */
	static private double getSettingDouble(String param) {
		return Double.valueOf(getSetting(param)); 
	}

	// ----------------------------------------------------------------------------

	/** Parse a boolean property */
	static private boolean getSettingBoolean(String name) {
		
		String value = getSetting(name);
		if (value.equals("1") || value.equalsIgnoreCase("true") || value.equalsIgnoreCase("t"))
			return true;
		else if (value.equals("0") || value.equalsIgnoreCase("false") || value.equalsIgnoreCase("f"))
			return false;
		else
			throw new IllegalArgumentException("Parse error for boolean parameter '" + name + "': expected '1' or '0', found '" + value + "'");
	}

	
	// ----------------------------------------------------------------------------

	/** Parse an int array property */
	static private ArrayList<Integer> getSettingIntArray(String name, boolean positiveSorted) {
		
		String[] propStr = getSetting(name).split(",");
		ArrayList<Integer> prop = new ArrayList<Integer>();
		
		if (propStr.length == 1 && propStr[0].compareTo("") == 0)
			return prop;

		for (int i=0; i<propStr.length; i++)
			prop.add(Integer.valueOf(propStr[i]));
			
		if (positiveSorted && !MagnumUtils.posIntIncreasing(prop))
			Magnum.error("Error parsing settings file, " + name + " has to be an ordered list of positive integers, given in increasing order");
		
		return prop;
	}


	// ----------------------------------------------------------------------------

	/** Parse an int array property */
	static private ArrayList<Double> getSettingDoubleArray(String name, boolean positiveSorted) {
		
		String[] propStr = getSetting(name).split(",");
		ArrayList<Double> prop = new ArrayList<Double>();
		
		if (propStr.length == 1 && propStr[0].compareTo("") == 0)
			return prop;

		for (int i=0; i<propStr.length; i++)
			prop.add(Double.valueOf(propStr[i]));
			
		if (positiveSorted && !MagnumUtils.posDoubleIncreasing(prop))
			Magnum.error("Error parsing settings file, " + name + " has to be an ordered list of positive numbers, given in increasing order");
		
		return prop;
	}

	
	// ----------------------------------------------------------------------------

	/** Parse a string array property */
	@SuppressWarnings("unused")
	static private String[] getStringArraySetting(Properties set, String name) {
		
		return set.getProperty(name).split(",");
	}

}
