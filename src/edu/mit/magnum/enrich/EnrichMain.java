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
package edu.mit.magnum.enrich;

import java.io.File;
import java.util.ArrayList;

import cern.colt.matrix.DoubleMatrix2D;
import edu.mit.magnum.*;
import edu.mit.magnum.gene.GeneIdMapping;
import edu.mit.magnum.net.Network;
import edu.mit.magnum.netprop.*;

/**
 *
 */
public class EnrichMain {

	/** The magnum instance */
	private Magnum mag;

	/** The name, used for output files (default, extracted from gwas and functional data filenames) */
	protected String name_ = null; 
	/** Genes with their scores */
	protected GeneScoreList geneScores_ = null;
	/** The functional / network data */
	protected FunctionalData functData_ = null;
	/** Provides functionality for (within-degree) label permutation */
	protected LabelPermuter permuter_ = null;	

	/** The enrichment analyzer */
	private Enrichment enrichment_ = null;
	
	/** The network */
	private Network network;
	
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public EnrichMain(Magnum mag) {
				
		this.mag = mag;
		// Initialize gene mapping
		if (!mag.set.idTypeFunctionalData_.equalsIgnoreCase(mag.set.idTypeGeneScores_))
			GeneIdMapping.getInstance(mag).load(mag.set.geneIdMappingFile_);
		
		// Load the gene scores, excluding genes from the excludedGenesFile
		geneScores_ = new GeneScoreList(mag, mag.set.geneScoreFile_, mag.set.excludedGenesFile_);
		
		// Initialize functional data (kernel) -- compute it or load from file
		File functionalDataFile = mag.set.functionalDataFile_;
		if (functionalDataFile == null) {
			// Change output dir to kernel dir (used to check if kernels are present or to export them)
			File outDirBkp = mag.set.outputDirectory_;
			if (mag.set.networkKernelDir != null)
				mag.set.outputDirectory_ = mag.set.networkKernelDir;
			else
				mag.set.outputDirectory_ = new File(outDirBkp, "network_kernels");
			
			// The default file
			File defaultFile = PstepKernel.getKFile(mag);
			// If it exists, load it
			if (defaultFile.exists() && mag.set.usePrecomputedKernels) {
				functionalDataFile = defaultFile;
			
			// Else, compute kernel
			} else {
				DoubleMatrix2D kernel = computeSimilarityNetwork();
				functData_ = new FunctionalData(mag, network, kernel, mag.set.excludedGenePairsFile_, geneScores_.getGenes());
				network = null; // Not needed anymore
				name_ = extractName(mag.set.geneScoreFile_, mag.set.networkFile_);
			}
			// Change the output dir back
			mag.set.outputDirectory_ = outDirBkp;
		}
		if (functionalDataFile != null) {
			if (!functionalDataFile.exists())
				throw new RuntimeException("File not found: " + functionalDataFile.getPath());
			functData_ = new FunctionalData(mag, functionalDataFile, mag.set.excludedGenePairsFile_, mag.set.functionalDataCols_, geneScores_.getGenes());
			name_ = extractName(mag.set.geneScoreFile_, functionalDataFile);
		}

		// The genes that were not loaded because they have no scores
		ArrayList<String> genesMissingScores = functData_.getGenesMissingScores();
		// Remove gene scores that are not in the funct data
		ArrayList<String> genesMissingFunctData = geneScores_.intersect(functData_.getGenes().keySet());
		
		// Exclude neighboring gene pairs (this is here because we first want to remove scores that are not in funct data (see previous line)
		if (mag.set.excludedGenesDistance_ >= 0)
			functData_.excludeNeighbors(geneScores_.getGenes());
		
		// This is only true if we don't map gene scores to entrez, otherwise the same entrez gene can appear multiple
		// times in the ranked score list because of many-to-many mappings.
		//assert geneScores_.getNumGenes() == functData_.getNumGenes();
		int numDuplicateGenes = geneScores_.getNumGenes() - functData_.getNumGenes();
		if (numDuplicateGenes < 0)
			throw new RuntimeException("More genes in functional data than gene score list");
		else if (numDuplicateGenes > 0)
			mag.log.warning(numDuplicateGenes + " duplicate genes in gene score list because of many-to-many mappings from Ensembl to Entrez gene IDs");
		
		// Print info on gene overlap
		mag.log.println("\nLoaded:");
		mag.log.println("- " + functData_.getNumGenes() + " genes with both GWAS and network data");
		mag.log.printlnVerbose("Removed:");
		mag.log.printlnVerbose("- " + genesMissingScores.size() + " genes have no GWAS data or were specifically excluded (e.g., HLA genes)");
		mag.log.printlnVerbose("- " + genesMissingFunctData.size() + " genes with no network data");
		mag.log.println();

//		// Export genes without scores
//		if (genesMissingScores.size() > 0) {
//			FileExport writer = new FileExport(Settings.outputDirectory_ + "/genesWithoutScores.txt");
//			for (String gene : genesMissingScores)
//				writer.println(gene);
//			writer.close();
//		}
//
//		// Export genes without funct data
//		if (genesMissingFunctData.size() > 0) {
//			FileExport writer = new FileExport(Settings.outputDirectory_ + "/genesWithoutFunctionalData.txt");
//			for (String gene : genesMissingFunctData)
//				writer.println(gene);
//			writer.close();
//		}
	}

	
	// ----------------------------------------------------------------------------

	/** Run enrichment analysis */
	public void run() {
		
		int numScoresPerGene = geneScores_.getNumScoresPerGene();
		
		// Compute enrichment for each gene score
		if (mag.set.geneScoreIndexEnd_ >= numScoresPerGene)
			throw new IllegalArgumentException("geneScoreIndexEnd_=" + mag.set.geneScoreIndexEnd_ + " >= numScoresPerGene=" + numScoresPerGene);
		if (mag.set.geneScoreIndexEnd_ < mag.set.geneScoreIndexStart_)
			throw new IllegalArgumentException("Settings.geneScoreIndexEnd_ < Settings.geneScoreIndexStart_");
		
		for (int i=mag.set.geneScoreIndexStart_; i<=mag.set.geneScoreIndexEnd_; i++) {
			// Rank genes according to i'th score
			geneScores_.sortGeneList(i);
			
			// Run enrichment
			if (functData_.getIsPairwiseData())
				runPairwiseData(i);
			else
				runPerGeneData(i);
		}
	}

	
	// ============================================================================
	// PRIVATE METHODS
		
	/** Compute kernel / tanimoto similarity network (when no precomputed similarity matrix is given) */
	private DoubleMatrix2D computeSimilarityNetwork() {
		
		// Load the input network
		File networkFile; 
		if (mag.set.networkDir_ != null)
			networkFile = new File(mag.set.networkDir_, mag.set.networkFile_.getPath());
		else
			networkFile = mag.set.networkFile_;
		
		mag.set.exportPairwiseNodeProperties_ = mag.set.exportKernels;
		mag.set.exportNodeProperties_ = false;
		mag.set.compressFiles_ = true;
		
		// Create kernel dir
		if (mag.set.exportKernels)
			mag.set.outputDirectory_.mkdirs();

		// Doesn't make sense to save more than one step
		if (mag.set.pstepKernelP_.size() != 1)
			throw new RuntimeException("Specify only one pstepKernelP when computing kernels on the fly within enrichment analysis");
			
		// Load network: p-step kernel only defined for undirected networks without self-loops
		mag.log.println();
		network = new Network(mag, networkFile, false, true, mag.set.isWeighted_, mag.set.threshold_);

		mag.log.println("COMPUTING RANDOM-WALK KERNEL");
		mag.log.println("----------------------------\n");
		
		// Compute and return K
		PairwiseProperties netprop = new PstepKernel(mag, network, mag.set.pstepKernelAlpha_, mag.set.pstepKernelP_, mag.set.pstepKernelNormalize_, false);
		netprop.run();
		DoubleMatrix2D kernel = netprop.getK();
		return kernel;		
	}


	// ----------------------------------------------------------------------------

	/** Run enrichment analysis for per gene functional data */
	private void runPerGeneData(int geneScoreIndex) {

		// For each property, run enrichment analysis and save results
		ArrayList<String> colNames = functData_.getColNames();
		for (int i=0; i<colNames.size(); i++) {
			// Reinitialize correct funct data indexes (they have been shuffled at prev iteration)
			permuter_ = new LabelPermuter(mag, functData_, geneScores_.getGenes(), mag.set.numBins_, i);

			// Compute mean enrichment
			enrichment_ = new EnrichmentIndividual(mag, functData_, geneScores_, permuter_, i);
			enrichment_.run();

			// Save results
			String col = colNames.get(i).replace("_weighted", "");
			String filename = new File(mag.set.outputDirectory_, name_ + "_" + col + "_bin" + mag.set.numBins_).getPath();
			if (geneScoreIndex > 0)
				filename += "." + geneScoreIndex;
			enrichment_.save(filename);
			mag.log.println();
		}
	}
	
	
	// ----------------------------------------------------------------------------

	/** Run enrichment analysis for pairwise functional data, e.g. network kernel */
	private void runPairwiseData(int geneScoreIndex) {
		
		// Map gwas genes to functional data
		permuter_ = new LabelPermuter(mag, functData_, geneScores_.getGenes(), mag.set.numBins_);

		enrichment_ = new EnrichmentPairwise(mag, functData_, geneScores_, permuter_);
		enrichment_.run();
		
		// Results
		mag.log.println("RESULTS\n" + 
				        "-------\n");
		// Save results
		mag.log.println("Writing enrichment curves and AUCs ...");
		String filename = new File(mag.set.outputDirectory_, name_).getPath();
		if (geneScoreIndex > 0)
			filename += "." + geneScoreIndex;
		enrichment_.save(filename);
		mag.log.println();

		// Print pvals
		enrichment_.printPvals();
	}

	
	// ----------------------------------------------------------------------------
	
	/** Extract the name for this run from the snp and functional data files */
	private String extractName(File geneScoreFile, File functDataFile) {
		
		String gwasName = mag.utils.extractBasicFilename(geneScoreFile.getName(), false);
		String functName = mag.utils.extractBasicFilename(functDataFile.getName(), false);
		
		functName = functName.replace("_nodeProperties", "");
		//functName = functName.replace("_undir", "");
		
		return gwasName + "--" + functName;
	}
	
	
	// ============================================================================
	// SETTERS AND GETTERS
	
	public Enrichment getEnrichment() { return enrichment_; }
	public GeneScoreList getGeneScores() { return geneScores_; }
	
}
