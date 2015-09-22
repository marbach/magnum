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
	/** The corresponding kernel (as is before any filtering for genes with scores) */
	private DoubleMatrix2D kernel;
	
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public EnrichMain() {
				
		// Initialize gene mapping
		if (!Magnum.set.idTypeFunctionalData_.equalsIgnoreCase(Magnum.set.idTypeGeneScores_))
			GeneIdMapping.getInstance().load(Magnum.set.geneIdMappingFile_);
		
		// Load the gene scores, excluding genes from the excludedGenesFile
		geneScores_ = new GeneScoreList(Magnum.set.geneScoreFile_, Magnum.set.excludedGenesFile_);
		//LinkageDisequilibrium.exportNeighboringGeneWindows(geneScores_.getGenes());

		// Initialize functional data (kernel) -- compute it or load from file
		File functionalDataFile = Magnum.set.functionalDataFile_;
		if (Magnum.set.functionalDataFile_ == null || !Magnum.set.functionalDataFile_.exists()) {
			// Compute kernel
			computeSimilarityNetwork();
			functData_ = new FunctionalData(network, kernel, Magnum.set.excludedGenePairsFile_, geneScores_.getGenes());
			name_ = extractName(Magnum.set.geneScoreFile_, Magnum.set.networkFile_);

		} else {
			functData_ = new FunctionalData(functionalDataFile, Magnum.set.excludedGenePairsFile_, Magnum.set.functionalDataCols_, geneScores_.getGenes());
			name_ = extractName(Magnum.set.geneScoreFile_, functionalDataFile);
		}
		
		// The genes that were not loaded because they have no scores
		ArrayList<String> genesMissingScores = functData_.getGenesMissingScores();
		// Remove gene scores that are not in the funct data
		ArrayList<String> genesMissingFunctData = geneScores_.intersect(functData_.getGenes().keySet());
		
		// Exclude neighboring gene pairs (this is here because we first want to remove scores that are not in funct data (see previous line)
		if (Magnum.set.excludedGenesDistance_ >= 0)
			functData_.excludeNeighbors(geneScores_.getGenes());
		
		// This is only true if we don't map gene scores to entrez, otherwise the same entrez gene can appear multiple
		// times in the ranked score list because of many-to-many mappings.
		//assert geneScores_.getNumGenes() == functData_.getNumGenes();
		int numDuplicateGenes = geneScores_.getNumGenes() - functData_.getNumGenes();
		if (numDuplicateGenes < 0)
			Magnum.log.error("More genes in functional data than gene score list");
		else if (numDuplicateGenes > 0)
			Magnum.log.warning(numDuplicateGenes + " duplicate genes in gene score list because of many-to-many mappings from Ensembl to Entrez gene IDs");
		
		// Print info on gene overlap
		Magnum.log.println();
		Magnum.log.println("Loaded:");
		Magnum.log.println("- " + functData_.getNumGenes() + " genes with scores and functional data");
		Magnum.log.println("Removed:");
		Magnum.log.println("- " + genesMissingScores.size() + " genes with functional data but no scores");
		Magnum.log.println("- " + genesMissingFunctData.size() + " genes with scores but no functional data");
		
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
		Magnum.log.println();
	}

	
	// ----------------------------------------------------------------------------

	/** Run enrichment analysis */
	public void run() {
		
		int numScoresPerGene = geneScores_.getNumScoresPerGene();
		
		// Compute enrichment for each gene score
		if (Magnum.set.geneScoreIndexEnd_ >= numScoresPerGene)
			throw new IllegalArgumentException("geneScoreIndexEnd_=" + Magnum.set.geneScoreIndexEnd_ + " >= numScoresPerGene=" + numScoresPerGene);
		if (Magnum.set.geneScoreIndexEnd_ < Magnum.set.geneScoreIndexStart_)
			throw new IllegalArgumentException("Settings.geneScoreIndexEnd_ < Settings.geneScoreIndexStart_");
		
		for (int i=Magnum.set.geneScoreIndexStart_; i<=Magnum.set.geneScoreIndexEnd_; i++) {
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
	private void computeSimilarityNetwork() {
		
		// Load the input network
		File networkFile; 
		if (Magnum.set.networkDir_ != null)
			networkFile = new File(Magnum.set.networkDir_, Magnum.set.networkFile_.getPath());
		else
			networkFile = Magnum.set.networkFile_;
		
		// TODO make this a setting in the app
		//Magnum.set.exportPairwiseNodeProperties_ = true;
		Magnum.set.exportPairwiseNodeProperties_ = false;
		Magnum.set.exportNodeProperties_ = false;
		Magnum.set.compressFiles_ = true;

		// Doesn't make sense to save more than one step
		if (Magnum.set.pstepKernelP_.size() != 1)
			throw new RuntimeException("Specify only one pstepKernelP when computing kernels on the fly within enrichment analysis");
			
		// Load network: p-step kernel only defined for undirected networks without self-loops
		Magnum.log.println();
		network = new Network(networkFile, false, true, Magnum.set.isWeighted_, Magnum.set.threshold_);

		Magnum.log.println("COMPUTING RANDOM-WALK KERNEL");
		Magnum.log.println("----------------------------\n");
		
		// Compute and save K
		PairwiseProperties netprop = new PstepKernel(network, Magnum.set.pstepKernelAlpha_, Magnum.set.pstepKernelP_, Magnum.set.pstepKernelNormalize_, false);
		netprop.run();
		kernel = netprop.getK();
		
		// filename where K was written
		//String filename = MagnumUtils.extractBasicFilename(network.getFile().getName(), false);
		//File kfile = new File(Magnum.set.outputDirectory_, filename + "_" + netprop.getName() + ".txt.gz");
		//return kfile;
	}


	// ----------------------------------------------------------------------------

	/** Run enrichment analysis for per gene functional data */
	private void runPerGeneData(int geneScoreIndex) {

		// For each property, run enrichment analysis and save results
		ArrayList<String> colNames = functData_.getColNames();
		for (int i=0; i<colNames.size(); i++) {
			// Reinitialize correct funct data indexes (they have been shuffled at prev iteration)
			permuter_ = new LabelPermuter(functData_, geneScores_.getGenes(), Magnum.set.numBins_, i);

			// Compute mean enrichment
			enrichment_ = new EnrichmentIndividual(functData_, geneScores_, permuter_, i);
			enrichment_.run();

			// Save results
			String col = colNames.get(i).replace("_weighted", "");
			String filename = new File(Magnum.set.outputDirectory_, name_ + "_" + col + "_bin" + Magnum.set.numBins_).getPath();
			if (geneScoreIndex > 0)
				filename += "." + geneScoreIndex;
			enrichment_.save(filename);
			Magnum.log.println();
		}
	}
	
	
	// ----------------------------------------------------------------------------

	/** Run enrichment analysis for pairwise functional data, e.g. network kernel */
	private void runPairwiseData(int geneScoreIndex) {
		
		// Map gwas genes to functional data
		permuter_ = new LabelPermuter(functData_, geneScores_.getGenes(), Magnum.set.numBins_);

		enrichment_ = new EnrichmentPairwise(functData_, geneScores_, permuter_);
		enrichment_.run();
		
		// Print pvals
		enrichment_.printPvals();
		
		// Save results
		String filename = new File(Magnum.set.outputDirectory_, name_).getPath();
		if (geneScoreIndex > 0)
			filename += "." + geneScoreIndex;
		enrichment_.save(filename);
		Magnum.log.println();
	}

	
	// ----------------------------------------------------------------------------
	
	/** Extract the name for this run from the snp and functional data files */
	private String extractName(File geneScoreFile, File functDataFile) {
		
		String gwasName = MagnumUtils.extractBasicFilename(geneScoreFile.getName(), false);
		String functName = MagnumUtils.extractBasicFilename(functDataFile.getName(), false);
		
		functName = functName.replace("_nodeProperties", "");
		//functName = functName.replace("_undir", "");
		
		return gwasName + "--" + functName;
	}
	
	
	// ============================================================================
	// SETTERS AND GETTERS
	
	public Enrichment getEnrichment() { return enrichment_; }
	public GeneScoreList getGeneScores() { return geneScores_; }
	
}
