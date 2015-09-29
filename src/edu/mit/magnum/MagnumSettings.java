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
 */
public class MagnumSettings extends Settings {	
	
	/** The annotation file (default file included in jar) */
	public final String annotationRsc = "edu/mit/magnum/gene/rsc/gene_coord.bed";
	/** The HLA+TF genes file (default file included in jar) */
	public final String hlaTfsRsc = "edu/mit/magnum/gene/rsc/tfs_mhc.txt";
	/** The TF genes file (default file included in jar) */
	public final String tfsRsc = "edu/mit/magnum/gene/rsc/tfs.txt";

	/** Colt Mersenne Twister random engine (should be used by all other random number generators) */
	public MersenneTwister mersenneTwisterRng_;
	/** Apache Commons random engine */
	public Well19937c wellRng_;
	/** Java random engine */
	public Random jdkRng_;

	// ----------------------------------------------------------------------------
	// VARIOUS
	
	/** Mode: 1 => Network analysis; 2 => Enrichment analysis */
	public int mode_;
	/** PRIVATE, NEEDS TO BE SET WITH setRandomSeed(), which initializes the random number generators. Set to -1 to use current time */
	private int randomSeed_;
	/** Output directory to save stuff */
	public File outputDirectory_;
	/** Output filename */
	public String outputFilename_;
	/** Compress output files (gzip) */
	public boolean compressFiles_;
	/** Verbose console output */
	public boolean verbose_;

	// ----------------------------------------------------------------------------
	// NETWORK PROPERTIES

	// INPUT NETWORK
	/** Directory containing the networks */
	public File networkDir_;
	/** The input network file */
	public File networkFile_;
	/** Delimiter used to separate columns (default 'tab' */
	public String networkFileDelim_;  
	/** Defines if the network should be interpreted as directed or undirected */
	public boolean isDirected_;
	/** Defines if self loops should be removed from the network */
	public boolean removeSelfLoops_;
	/** Set true to treat the network as weighted */
	public boolean isWeighted_;
	/** Threshold for including weighted edges */
	public double threshold_;
	/** Exclude "super-hubs" that connect to more than the given fraction of genes (set 1 to include all) */
	public double superHubThreshold_;
	/** Optional file specifying a set of reference nodes */
	public File refNodesFile_;

	// NETWORKOPS
	/** Take union (max edge) over all networks in networkDir or the sets specified in the file below */
	public boolean computeUnion_;
	/** Define the network sets that should be combined (leave empty to combine all networks) */
	public File networkGroupFile_;
	/** Prefix of the files in the network dir */
	public String networkFilePrefix_;
	/** Add networks of the same cell type */
	public boolean computePairwiseSum_;
	/** The network directory of the second networks */
	public File networkDir2_;
	
	// BASIC NETWORK PROPERTIES
	/** Node degree (directed networks, also indegree and outdegree) */
	public boolean computeDegree_;
	/** Node betweenness centrality (edge directionality observed for directed networks) */
	public boolean computeBetweenness_;
	/** Node clustering coefficient (edge directionality observed for directed networks) */
	public boolean computeClusteringCoefficient_;
	/** For each node, distance to all other nodes (or all reference nodes) and closeness centrality */
	public boolean computeShortestPathLengths_;

	// KERNELS
	/** P-step random walk kernel (Smola & Kondor, 2003) */
	public boolean computePstepKernel_;
	/** alpha parameter of p-step random walk kernel (alpha >= 2) */
	public double pstepKernelAlpha_;
	/** Number of steps p of random walk kernel (p >= 1) */
	public ArrayList<Integer> pstepKernelP_;
	/** Normalize the kernel matrix (divide by the max) */
	public boolean pstepKernelNormalize_;
	
	// TANIMOTO COEFFICIENT
	/** Tanimoto coefficient between target genes */
	public boolean computeTargetTanimoto_;
	/** Tanimoto coefficient between regulators */
	public boolean computeTfTanimoto_;

	// OUTPUT FILES
	/** A suffix/ending that is appended to all output files for this run (use to distinguish output files from multiple runs) */
	public String outputSuffix_ = "";
	/** Export all computed pairwise node properties (e.g., similarity, distance matrices) */
	public boolean exportPairwiseNodeProperties_;
	/** Export all computed node properties (e.g., avg. similarity, distance for each node) */
	public boolean exportNodeProperties_;

	// ----------------------------------------------------------------------------
	// GENOME ANNOTATION

	/** Set of genes to be considered (leave empty to use all genes from the annotation file) */
	public String genesToBeLoadedFile_;

	/** The chromosome to be considered (chr1, ..., chr22, chrX, chrY), leave empty for all chromosomes */ 
	public String chromosome_;
	/** Ignore sex chromosomes */
	public boolean excludeXYChromosomes_;
	/** Ignore HLA genes (only works for default annotation, custom annotations need to provide their own excludedGenesFile) */
	public boolean excludeHlaGenes_;
	
	/** The file with the gencode annotation */
	public File gencodeAnnotationFile_;
	/** UCSC genome browser annotation (use for Entrez IDs) */ 
	public File ucscAnnotationFile_;
	/** Set true to load only protein-coding genes */
	public boolean loadOnlyProteinCodingGenes_;

	/** Mapping file to convert Entrez IDs, ENSEMBL IDs and gene symbols */
	public String geneIdMappingFile_;
		
	// ----------------------------------------------------------------------------
	// ENRICHMENT CURVES
	
	// INPUT
	/** The gene coordinates (custom annotation) */
	public File geneCoordFile_;
	/** The gene scores */
	public File geneScoreFile_;
	/** Cutoff for genome-wide significance of gene scores */
	public double genomeWideSignificanceThreshold_;
	/** Exclude genome-wide significant genes (below threshold) */
	public boolean excludeGenomeWideSignificantGenes_;

	/** The file with the functional data, e.g. network kernels (cols: gene id, property 1, property 2, ...) */ 
	public File functionalDataFile_; 
	/** Specify which columns should be loaded (-1: all columns; 1: first gene property column) */
	public ArrayList<Integer> functionalDataCols_;

	/** Genes to be excluded from enrichment analysis (e.g., MHC region) */
	public File excludedGenesFile_;
	/** Gene pairs to be excluded from enrichment analysis (e.g., genes in LD) */
	public File excludedGenePairsFile_;
	/** Exclude gene pairs with windows smaller than the given distance apart (given in megabases; -1: no exclusion; 1000000: all genes on same chromosome) */
	public double excludedGenesDistance_; 

	/** Gene IDs used in geneScoreFile, excludedGenesFile, excludedGenePairsFile ('ensembl', 'entrez', 'hugo') */
	public String idTypeGeneScores_;
	/** Gene IDs used in functionalDataFile */
	public String idTypeFunctionalData_;

	/** Use precomputed network kernels if available in networkKernelDir */
	public boolean usePrecomputedKernels = true;
	/** Directory for network kernels (default: <outputDir>/network_kernels/) */
	public File networkKernelDir;
	/** Save network kernels for use in subsequent runs (takes a lot of disk space!) */
	public boolean exportKernels;

	// PARAMETERS
	/** Number of random permutations to compute confidence intervals */
	public int numPermutations_;
	/** The number of bins for within-degree permutation */
	public int numBins_;
	/** Scale kernels: K'(i,j) = K(i,j)/sqrt(rowSums(K)[i] * colSums(K)[j]) */
	public boolean scaleKernel_;

	/** Equidistant curve resolution, e.g., set 10 to compute every 10th point on the curves */
	public int constCurveResolution_;
	/** Varying curve resolution, e.g., set 2 to compute points: 2, 6, 12, 20, 30, 42, ... (takes precedence over constCurveResolution, set -1 to disable) */
	public int varCurveResolution_;
	/** Compute curves only for the top part of the list (e.g., 0.1 for top 10%) */
	public double curveCutoff_;
	/** Sliding window size */
	public int slidingWindowSize_;

	/** Draw boundaries for given p-values (e.g., set 0.05 to draw the upper/lower boundary where only 5% of random curves above/below) */ 
	public ArrayList<Double> pval_;
	/** Indicates whether the boundaries are one-sided or two-sided (equivalent to dividing pval by two) */
	public boolean twoSidedTest_;
	/** Control FDR over all points of the curve together (i.e., correct for multiple hypothesis testing across curve) */
	public boolean controlFDR_;
	/** The top X snps will be ignored for FDR control (too noisy at start of the list) */
	public int FDRStart_;
	/** Start to integrate the AUC only at k=10 (reduce noise at the start of the list) */
	public int AUCStart_;
	/** Index of gene scores for which enrichment is to be computed, e.g., (0,9) for the first ten gene scores */ 
	public int geneScoreIndexStart_;
	public int geneScoreIndexEnd_;

	
	// OUTPUT FILES
	/** Number of random permutations for which enrichment curves are exported (smaller or equal numPermutations) */
	public int numPermutationsExport_;

	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public MagnumSettings(Magnum mag) {
		super(mag);
		resetToDefaults();
	}
	
	
	// ----------------------------------------------------------------------------
	
	/** Set default values for all settings */
	public void resetToDefaults() {
		
		// Initializes the RNGs
		setRandomSeed(42);

		mode_ = 0;
		outputDirectory_ = new File(System.getProperty("user.dir"));
		outputFilename_ = "";
		compressFiles_ = true;
		verbose_ = true;

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
		excludeXYChromosomes_ = true;
		excludeHlaGenes_ = true;
		
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

		usePrecomputedKernels = true;
		networkKernelDir = null;
		exportKernels = false;

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
	public void loadSettings(String settingsFile, boolean requireAll) {
		
		mag.log.println("SETTINGS FILE");
		mag.log.println("-------------\n");
				
		try {
			// Open the input stream
			//InputStream in = MagnumSettings.class.getClassLoader().getResourceAsStream("edu/mit/magnum/settings.txt");
				
			// Check that the specified settings file exists
			if (settingsFile == null || settingsFile.isEmpty())
				throw new RuntimeException("No settings file specified");
			else if (!new File(settingsFile).exists())
				throw new RuntimeException("Settings file not found: " + settingsFile);

			// Open file input stream
			mag.log.println("- Loading settings file: " + settingsFile + "\n");
			InputStream in = new FileInputStream(settingsFile);

			// Load the settings
			prop = new Properties();
			prop.load(new InputStreamReader(in));
			
			// Get the param values
			setParameterValues(requireAll);
			
		} catch (Exception e) {
			mag.log.warning(e.getMessage());
			throw new RuntimeException("Failed to load settings file (a parameter may be missing or malformed): " + settingsFile);
		}		
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
	private void setParameterValues(boolean requireAll) throws Exception {

		// VARIOUS
		if (requireAll || prop.containsKey("mode"))
			mode_ = getSettingInt("mode");
		if (requireAll || prop.containsKey("randomSeed"))
			setRandomSeed(getSettingInt("randomSeed"));
		if (requireAll || prop.containsKey("outputDirectory")) {
			outputDirectory_ = getFileSetting("outputDirectory");
			if (outputDirectory_.equals("")) 
				outputDirectory_ = new File(System.getProperty("user.dir"));
		}
		if (requireAll || prop.containsKey("outputFilename"))
			outputFilename_ = getSetting("outputFilename");
		if (requireAll || prop.containsKey("verbose"))
			verbose_ = getSettingBoolean("verbose");
		mag.log.setVerbose(verbose_);

		// INPUT NETWORK
		if (requireAll || prop.containsKey("networkDir"))
			networkDir_ = getFileSetting("networkDir");
		if (requireAll || prop.containsKey("networkFile"))
			networkFile_ = getFileSetting("networkFile");
		if (requireAll || prop.containsKey("networkFileDelim"))
			networkFileDelim_ = getSetting("networkFileDelim");
		if (requireAll || prop.containsKey("isDirected"))
			isDirected_ = getSettingBoolean("isDirected");
		if (requireAll || prop.containsKey("removeSelfLoops"))
			removeSelfLoops_ = getSettingBoolean("removeSelfLoops");
		if (requireAll || prop.containsKey("isWeighted"))
			isWeighted_ = getSettingBoolean("isWeighted");
		if (requireAll || prop.containsKey("threshold"))
			threshold_ = getSettingDouble("threshold");
		if (requireAll || prop.containsKey("superHubThreshold"))
			superHubThreshold_ = getSettingDouble("superHubThreshold");
		if (requireAll || prop.containsKey("refNodesFile"))
			refNodesFile_ = getFileSetting("refNodesFile");

		// OUTPUT FILES
		if (requireAll || prop.containsKey("outputSuffix"))
			outputSuffix_ = getSetting("outputSuffix");
		if (requireAll || prop.containsKey("exportPairwiseNodeProperties"))
			exportPairwiseNodeProperties_ = getSettingBoolean("exportPairwiseNodeProperties");
		if (requireAll || prop.containsKey("exportNodeProperties"))
			exportNodeProperties_ = getSettingBoolean("exportNodeProperties");
		if (requireAll || prop.containsKey("compressFiles"))
			compressFiles_ = getSettingBoolean("compressFiles");

		// NETWORKOPS
		if (requireAll || prop.containsKey("computeUnion"))
			computeUnion_ = getSettingBoolean("computeUnion");
		if (requireAll || prop.containsKey("networkGroupFile"))
			networkGroupFile_ = getFileSetting("networkGroupFile");
		if (requireAll || prop.containsKey("networkFilePrefix"))
			networkFilePrefix_ = getSetting("networkFilePrefix");
		if (requireAll || prop.containsKey("computePairwiseSum"))
			computePairwiseSum_ = getSettingBoolean("computePairwiseSum");
		if (requireAll || prop.containsKey("networkDir2"))
			networkDir2_ = getFileSetting("networkDir2");

		// BASIC NETWORK PROPERTIES
		if (requireAll || prop.containsKey("computeDegree"))
			computeDegree_ = getSettingBoolean("computeDegree");
		if (requireAll || prop.containsKey("computeBetweenness"))
			computeBetweenness_ = getSettingBoolean("computeBetweenness");
		if (requireAll || prop.containsKey("computeClusteringCoefficient"))
			computeClusteringCoefficient_ = getSettingBoolean("computeClusteringCoefficient");

		// SHORTEST PATHS
		if (requireAll || prop.containsKey("computeShortestPathLengths"))
			computeShortestPathLengths_ = getSettingBoolean("computeShortestPathLengths");

		// KERNELS
		if (requireAll || prop.containsKey("computePstepKernel"))
			computePstepKernel_ = getSettingBoolean("computePstepKernel");
		if (requireAll || prop.containsKey("pstepKernelAlpha"))
			pstepKernelAlpha_ = getSettingDouble("pstepKernelAlpha");
		if (requireAll || prop.containsKey("pstepKernelP"))
			pstepKernelP_ = getSettingIntArray("pstepKernelP", true);
		if (requireAll || prop.containsKey("pstepKernelNormalize"))
			pstepKernelNormalize_ = getSettingBoolean("pstepKernelNormalize");

		// TANIMOTO
		if (requireAll || prop.containsKey("computeTargetTanimoto"))
			computeTargetTanimoto_ = getSettingBoolean("computeTargetTanimoto");
		if (requireAll || prop.containsKey("computeTfTanimoto"))
			computeTfTanimoto_ = getSettingBoolean("computeTfTanimoto");

		// ----------------------------------------------------------------------------
		// GENOME ANNOTATION

		if (requireAll || prop.containsKey("genesToBeLoadedFile"))
			genesToBeLoadedFile_ = getSetting("genesToBeLoadedFile");
		if (requireAll || prop.containsKey("chromosome"))
			chromosome_ = getSetting("chromosome");
		if (requireAll || prop.containsKey("excludeXYChromosomes"))
			excludeXYChromosomes_ = getSettingBoolean("excludeXYChromosomes");
		if (requireAll || prop.containsKey("excludeHlaGenes"))
			excludeHlaGenes_ = getSettingBoolean("excludeHlaGenes");
		if (requireAll || prop.containsKey("genecodeAnnotationFile"))
			gencodeAnnotationFile_ = getFileSetting("genecodeAnnotationFile");
		if (requireAll || prop.containsKey("ucscAnnotationFile"))
			ucscAnnotationFile_ = getFileSetting("ucscAnnotationFile");
		if (requireAll || prop.containsKey("loadOnlyProteinCodingGenes"))
			loadOnlyProteinCodingGenes_ = getSettingBoolean("loadOnlyProteinCodingGenes");
		if (requireAll || prop.containsKey("geneIdMappingFile"))
			geneIdMappingFile_ = getSetting("geneIdMappingFile");

		// ----------------------------------------------------------------------------
		// ENRICHMENT ANALYSIS

		if (requireAll || prop.containsKey("geneCoordFile")) {
			geneCoordFile_ = getFileSetting("geneCoordFile");
			idTypeFunctionalData_ = "custom";
			idTypeGeneScores_ = "custom";
		}

		if (requireAll || prop.containsKey("geneScoreFile"))
			geneScoreFile_ = getFileSetting("geneScoreFile");
		if (requireAll || prop.containsKey("genomeWideSignificanceThreshold"))
			genomeWideSignificanceThreshold_ = getSettingDouble("genomeWideSignificanceThreshold");
		if (requireAll || prop.containsKey("excludeGenomeWideSignificantGenes"))
			excludeGenomeWideSignificantGenes_ = getSettingBoolean("excludeGenomeWideSignificantGenes");

		if (requireAll || prop.containsKey("functionalDataFile"))
			functionalDataFile_ = getFileSetting("functionalDataFile"); 
		if (requireAll || prop.containsKey("functionalDataCols"))
			functionalDataCols_ = getSettingIntArray("functionalDataCols", true);		

		if (requireAll || prop.containsKey("excludedGenesFile"))
			excludedGenesFile_ = getFileSetting("excludedGenesFile");
		if (requireAll || prop.containsKey("excludedGenePairsFile"))
			excludedGenePairsFile_ = getFileSetting("excludedGenePairsFile");
		if (requireAll || prop.containsKey("excludedGenesDistance"))
			excludedGenesDistance_ = getSettingInt("excludedGenesDistance");

		if (requireAll || prop.containsKey("idTypeGeneScores"))
			idTypeGeneScores_ = getSetting("idTypeGeneScores");
		if (requireAll || prop.containsKey("idTypeFunctionalData"))
			idTypeFunctionalData_ = getSetting("idTypeFunctionalData");

		if (requireAll || prop.containsKey("usePrecomputedKernels"))
			usePrecomputedKernels = getSettingBoolean("usePrecomputedKernels");
		if (requireAll || prop.containsKey("networkKernelDir"))
			networkKernelDir = getFileSetting("networkKernelDir");
		if (requireAll || prop.containsKey("exportKernels"))
			exportKernels = getSettingBoolean("exportKernels");
		
		// ENRICHMENT
		if (requireAll || prop.containsKey("numPermutations"))
			numPermutations_ = getSettingInt("numPermutations");
		if (requireAll || prop.containsKey("numBins"))
			numBins_ = getSettingInt("numBins");
		if (requireAll || prop.containsKey("scaleKernel"))
			scaleKernel_ = getSettingBoolean("scaleKernel");

		if (requireAll || prop.containsKey("constCurveResolution"))
			constCurveResolution_ = getSettingInt("constCurveResolution");
		if (requireAll || prop.containsKey("varCurveResolution"))
			varCurveResolution_ = getSettingInt("varCurveResolution");
		if (requireAll || prop.containsKey("curveCutoff"))
			curveCutoff_ = getSettingDouble("curveCutoff");
		if (requireAll || prop.containsKey("slidingWindowSize"))
			slidingWindowSize_ = getSettingInt("slidingWindowSize");

		if (requireAll || prop.containsKey("numPermutationsExport")) {
			numPermutationsExport_ = getSettingInt("numPermutationsExport");
			if (numPermutationsExport_ > numPermutations_)
				throw new IllegalArgumentException("Invalid settings: 'numPermutationsExport' must be smaller or equal 'numPermutations'");
		}
		if (requireAll || prop.containsKey("pval"))
			pval_ = getSettingDoubleArray("pval", true);
		if (requireAll || prop.containsKey("twoSidedTest"))
			twoSidedTest_ = getSettingBoolean("twoSidedTest");
		if (requireAll || prop.containsKey("controlFDR"))
			controlFDR_ = getSettingBoolean("controlFDR");
		if (requireAll || prop.containsKey("FDRStart"))
			FDRStart_ = getSettingInt("FDRStart");
		if (requireAll || prop.containsKey("AUCStart"))
			AUCStart_ = getSettingInt("AUCStart");
		if (requireAll || prop.containsKey("geneScoreIndexStart"))
			geneScoreIndexStart_ = getSettingInt("geneScoreIndexStart");
		if (requireAll || prop.containsKey("geneScoreIndexEnd"))
			geneScoreIndexEnd_ = getSettingInt("geneScoreIndexEnd");
	}
	
	
	// ============================================================================
	// GETTERS AND SETTERS
	
	public int getRandomSeed() { return randomSeed_; }

	
}
