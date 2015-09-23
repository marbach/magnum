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
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import edu.mit.magnum.FileParser;
import edu.mit.magnum.Magnum;
import edu.mit.magnum.gene.Gene;
import edu.mit.magnum.gene.GeneAnnotation;
import edu.mit.magnum.gene.GeneAnnotationCustom;
import edu.mit.magnum.gene.GeneAnnotationGencode;
import edu.mit.magnum.gene.GeneIdMapping;

/**
 * 
 */
public class GeneScoreList {

	/** Genes to be excluded */
	private HashSet<String> excludedGenes_ = null;
	/** Genes with their scores */
	private ArrayList<Gene> rankedGenes_ = null;
	/** The number of genes with score below genome-wide significance threshold */
	private int numGenomeWideSignificant_ = -1;
	/** Number of scores per gene */
	private int numScoresPerGene_ = -1;
	/** Flag indicates if HLA genes have been excluded */
	private boolean excludingHlaGenes_ = false;

	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public GeneScoreList(File geneScoreFile, File excludedGenesFile) {
			
		// Load genes to be excluded
		loadExcludedGenes(excludedGenesFile);
		// Load genes and scores
		loadGeneScores(geneScoreFile);
		// Expand genes using the given window size and neighbor distance (used to exclude genes within a given distance) 
		if (Magnum.set.excludedGenesDistance_ > 0)
			expandGeneWindows();
	}

	
	// ----------------------------------------------------------------------------

	/** Remove genes that are not in the given gene set, return list of removed genes */
	public ArrayList<String> intersect(Set<String> geneSet) {
		
		// The genes not in the given set
		ArrayList<String> removedGenes = new ArrayList<String>();
		
		Iterator<Gene> iter = rankedGenes_.iterator();
		while (iter.hasNext()) {
			Gene gene = iter.next();
			
			if (!geneSet.contains(gene.id_)) {
				removedGenes.add(gene.id_);
				iter.remove();
			}
		}
		return removedGenes;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Sort the ranked gene list according to the i'th gene score */
	public void sortGeneList(int i) {
		
		// Comparator to sort genes by score
		final class GeneScoreComparator implements Comparator<Gene> {
			// The index of the gene score used for the comparision
			private int scoreIndex_ = -1;
			// Constructor
			GeneScoreComparator(int scoreIndex) {
				scoreIndex_ = scoreIndex;
			}
			// Compare
			@Override
			public int compare(Gene g1, Gene g2) {
				return Double.compare(g1.getScore(scoreIndex_), g2.getScore(scoreIndex_));
			}
		}
		Collections.sort(rankedGenes_, new GeneScoreComparator(i));
		
		// Count number of genes above genome-wide significance threshold
		numGenomeWideSignificant_ = 0;
		for (int k=0; k<rankedGenes_.size(); k++) {
			double score_k = rankedGenes_.get(k).getScore(i);
		
			if (score_k <= Magnum.set.genomeWideSignificanceThreshold_)
				numGenomeWideSignificant_++;
			else
				break;
		}
	}

		
	// ============================================================================
	// PRIVATE METHODS
		
	/** Load the gene scores, initialize excludedGenes_ */
	private void loadExcludedGenes(File excludedGenesFile) {
		
		excludedGenes_ = new HashSet<String>();
		excludingHlaGenes_ = false;
		GeneIdMapping idMapping = GeneIdMapping.getInstance();
		
		FileParser parser;
		if (excludedGenesFile == null) {
			// Use default
			if (Magnum.set.geneCoordFile_ != null)
				return;
			
			InputStream in;
			// TODO maybe allow also for no genes excluded
			if (Magnum.set.excludeHlaGenes_) {
				excludingHlaGenes_ = true;
				in = Magnum.class.getClassLoader().getResourceAsStream(Magnum.set.hlaTfsRsc);
			} else {
				in = Magnum.class.getClassLoader().getResourceAsStream(Magnum.set.tfsRsc);
			}
			parser = new FileParser(in);
		
		} else {
			parser = new FileParser(excludedGenesFile);
		}

		// Parse header
		String[] header = parser.readLine();
		if (!header[0].equals("gene_id"))
			Magnum.log.error("Expected header line with first field 'gene_id' (tab-separated)");
		
		// First line
		while (true) {
			String[] nextLine = parser.readLine();
			if (nextLine == null)
				break;
			
			// Gene id
			String id = nextLine[0];
			id = idMapping.removeEnsemblVersion(id);
			excludedGenes_.add(id);
		}
		parser.close();
	}

	
	// ----------------------------------------------------------------------------

	/** 
	 * Load the gene scores, initialize genes_. 
	 * Removes genes in excludedGenes_. 
	 * Translates the gene score ids to the idTypeFunctionalData if necessary. 
	 */
	private void loadGeneScores(File geneScoreFile) {
		
		GeneIdMapping idMapping = GeneIdMapping.getInstance();
		boolean translateToEntrez = Magnum.set.idTypeFunctionalData_.equalsIgnoreCase("entrez");

		// Load gene positions so that gene pairs on same chromosome can be excluded
		GeneAnnotation geneAnnot;
		if (Magnum.set.idTypeGeneScores_.equals("ensembl"))
			geneAnnot = new GeneAnnotationGencode(null, Magnum.set.loadOnlyProteinCodingGenes_);
		else if (Magnum.set.idTypeGeneScores_.equals("custom"))
			geneAnnot = new GeneAnnotationCustom(null);
		else
			throw new RuntimeException("Invalid idTypeGeneScores: " + Magnum.set.idTypeGeneScores_);
		geneAnnot.loadAnnotation();
		
		rankedGenes_ = new ArrayList<Gene>();
		FileParser parser = new FileParser(geneScoreFile);
		
		// Parse header
		String[] header = parser.readLine();
		
		int geneIdCol = -1;
		int geneSymbolCol = -1;
		int firstPvalCol = -1;
		int lastPvalCol = -1;
		for (int i=0; i<header.length; i++) {
			if (header[i].contains("pvalue")) {
				if (firstPvalCol == -1)
					firstPvalCol = i;
				lastPvalCol = i;
			} else if (header[i].equals("gene_symbol")) {
				geneSymbolCol = i;
			} else if (header[i].equals("gene_id")) {
				geneIdCol = i;
			}
		}
		if (geneIdCol == -1)
			throw new RuntimeException("Didn't find column 'gene_id' in header");
		if (firstPvalCol == -1)
			throw new RuntimeException("Didn't find column 'pvalue' in header");
		numScoresPerGene_ = lastPvalCol - firstPvalCol + 1;
		
		int numExcluded = 0;
		int numNoScore = 0;
		int numNoAnnot = 0;
		int numSignificantExcluded = 0;
		
		while (true) {
			// Read next line
			String[] nextLine = parser.readLine();
			if (nextLine == null)
				break;
			
			// Gene id
			String geneId = nextLine[geneIdCol];
			if (Magnum.set.idTypeGeneScores_.equals("ensembl"))
				geneId = idMapping.removeEnsemblVersion(geneId);
			
			// Check if it should be excluded
			if (excludedGenes_.contains(geneId)) {
				numExcluded++;
				continue;
			}

			// Parse scores
			double[] scores = new double[numScoresPerGene_];
			try {
				for (int i=0; i<numScoresPerGene_; i++)
					scores[i] = Double.parseDouble(nextLine[firstPvalCol + i]);
			} catch (Exception e) {
				numNoScore++;
				continue;
			}

			// Get the genome coordinates from the annotation
			Gene gencodeGene = geneAnnot.getGene(geneId);
			if (gencodeGene == null) {
				numNoAnnot++;
				continue;
				//gencodeGene = new Gene(null);
			}
			if (Magnum.set.excludeXYChromosomes_ && geneAnnot.isAllosome(gencodeGene.chr_))
				continue;
			
			// Exclude genome-wide significant genes
			if (Magnum.set.excludeGenomeWideSignificantGenes_ && scores[0] <= Magnum.set.genomeWideSignificanceThreshold_) {
				numSignificantExcluded++;
				continue;
			}
			
			HashSet<String> idSet;
			if (translateToEntrez) {
				idSet = idMapping.ensembl2entrez(geneId);
			} else {
				idSet = new HashSet<String>(1);
				idSet.add(geneId);
			}
			
			for (String id : idSet) {
				Gene gene = new Gene(id);
				gene.setPosition(gencodeGene.chr_, gencodeGene.start_, gencodeGene.end_, gencodeGene.posStrand_);
				rankedGenes_.add(gene);
			
				// Parse gene symbol
				if (geneSymbolCol != -1)
					gene.symbol_ = nextLine[geneSymbolCol];
			
				gene.setScore(scores);
			}
		}
		sortGeneList(0);
		
		if (Magnum.set.excludeXYChromosomes_)
			Magnum.log.println("- Excluding genes on X, Y chromosomes");
		if (excludingHlaGenes_)
			Magnum.log.println("- Excluding HLA genes");
		//if (numExcluded > 0)
		//	Magnum.log.println("- " + numExcluded + " genes excluded");
		if (numNoScore > 0)
			Magnum.log.println("- " + numNoScore + " genes with NA score");
		if (numNoAnnot > 0)
			Magnum.log.println("- " + numNoAnnot + " genes not found in the annotation");
		if (numSignificantExcluded > 0)
			Magnum.log.println("- " + numSignificantExcluded + " genes with p-val <" + Magnum.set.genomeWideSignificanceThreshold_ + " excluded");
	}
	
	
	// ----------------------------------------------------------------------------

	/** Expand genes using the given neighbor distance (used to exclude genes within a given distance) */
	private void expandGeneWindows() {
		
		int d = (int) (1000000 * Magnum.set.excludedGenesDistance_);

		for (Gene g : rankedGenes_)
			g.expand(d, d);
	}

	
	// ============================================================================
	// SETTERS AND GETTERS
	
	public ArrayList<Gene> getGenes() { return rankedGenes_; }
	public Gene getGene(int i) { return rankedGenes_.get(i); }
	public int getNumGenes() { return rankedGenes_.size(); }
	public int getNumScoresPerGene() { return numScoresPerGene_; }
	public int getNumGenomeWideSignificant() { return numGenomeWideSignificant_; }
}
