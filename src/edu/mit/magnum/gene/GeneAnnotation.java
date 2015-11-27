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
package edu.mit.magnum.gene;

import java.io.File;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import ch.unil.gpsutils.FileExport;
import ch.unil.gpsutils.FileParser;
import edu.mit.magnum.Magnum;


/**
 * Gene annotation, provides functionality to load IDs and coordinates of genes
 */
abstract public class GeneAnnotation {

	/** The magnum instance */
	protected Magnum mag;

	/** The file with the genome annotation */
	protected File annotationFile_ = null;
	
	/** User specified set of genes to be loaded (leave empty for all genes, boolean indicates if gene was found in the annotation) */
	protected HashMap<String, Boolean> genesToBeLoaded_ = null;
	/** Load only genes from this chromosome (leave empty for all chromosomes) */
	protected String chromosomeToBeLoaded_ = null;
	/** Flag indicates whether only protein coding genes should be loaded */
	protected boolean loadOnlyProteinCoding_ = false;
	
	/** The genes that were loaded from the annotation file */
	protected LinkedHashMap<String, Gene> genes_ = null;
	
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public GeneAnnotation(Magnum mag, File annotationFile, String chromosomeToBeLoaded, boolean loadOnlyProteinCoding) {
		
		this.mag = mag;
		annotationFile_ = annotationFile; 
		chromosomeToBeLoaded_ = chromosomeToBeLoaded;
		loadOnlyProteinCoding_ = loadOnlyProteinCoding;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Load the specified genes (genesToBeLoaded_) from the annotation file */
	public LinkedHashMap<String, Gene> loadAnnotation(String genesToBeLoadedFile) {
		
		// Load the set of genes to be considered
		loadGenesToBeLoaded(genesToBeLoadedFile);
		// Load the specified genes from the annotation
		loadAnnotation();
		
		// If not all specified genes were found, print warning
		if (genesToBeLoaded_ != null && genesToBeLoaded_.size() != genes_.size()) {
			int numNotFound = genesToBeLoaded_.size() - genes_.size();
			
			String genesNotFound = "";
			for (Entry<String,Boolean> entry : genesToBeLoaded_.entrySet())
				if (!entry.getValue())
					genesNotFound += entry.getKey() + " ";

			mag.log.println("   - " + numNotFound + " genes were not found in the annotation: " + genesNotFound);
		}
		return genes_;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Initialize the coordinates (chromosome and position) for the given genes */
	public void loadCoordinates(Collection<Gene> genes) {
		
		// Load the gene annotation
		loadAnnotation();
		
		// Annotate the given genes
		int notFound = 0;
		for (Gene g : genes) {
			Gene g_annot = genes_.get(g.id_);
			if (g_annot == null)
				notFound++;
			else
				g.setPosition(g_annot.chr_, g_annot.start_, g_annot.end_, g_annot.posStrand_);
		}
		if (notFound > 0)
			mag.log.warning(notFound + " genes not found in annotation");
	}


	// ----------------------------------------------------------------------------

	/** Load the specified genes (genesToBeLoaded_) from the annotation file */
	abstract public HashMap<String, Gene> loadAnnotation();

	
	// ----------------------------------------------------------------------------

	/** Load the specified set of genes (genesToBeLoaded_) */
	public void loadGenesToBeLoaded(String filename) {
				
		// Return if no gene file was specified (all genes from the annotation will be loaded)
		if (filename == null ||
				filename.equals(" ") ||
				filename.equals(""))
			return;

		genesToBeLoaded_ = new HashMap<String, Boolean>();
		FileParser parser = new FileParser(mag.log, filename);
		
		while (true) {
			String[] nextLine = parser.readLine();
			if (nextLine == null)
				break;
			if (nextLine.length != 1)
				parser.error("Expected one column (gene ID)");
			
			genesToBeLoaded_.put(nextLine[0].toUpperCase(), false);
		}
		parser.close();
	}


	// ----------------------------------------------------------------------------

	/** Write the genes to a file with gene id, symbol and position */
	public void writeGeneList(String filename) {
		
		FileExport writer = new FileExport(mag.log, filename);
		String prevChr = null;
		int prevStart = -1;
		
		String header = "geneId\tgeneSymbol\tchromosome\tstart\tend\tstrand";
		writer.println(header);
		
		for (Gene gene : genes_.values()) {
			if (prevStart == -1 || prevChr != gene.chr_) {
				prevStart = gene.start_;
				prevChr = gene.chr_;
			}
			
			if (gene.start_ < prevStart)
				throw new RuntimeException("Genes are not ordered by genomic position");
			
			String nextLine = gene.id_ + "\t" + gene.symbol_ + "\t" + 
					gene.chr_ + "\t" + gene.start_ + "\t" + gene.end_ + "\t" + 
					(gene.posStrand_ ? "+" : "-");
			
			writer.println(nextLine);
		}
		writer.close();
	}
	

	// ----------------------------------------------------------------------------

	/** Find the gene with the given id, returns null if not found */
	public Gene getGene(String id) {
		
		return genes_.get(id);
	}

		
	// ----------------------------------------------------------------------------

	/** Return true for '+', false for '-' and error otherwise */
	public boolean isAllosome(String ch) {
	
		return ch.equals("chrX") || ch.equals("chrY");
	}


		
	// ============================================================================
	// PROTECTED METHODS

	/** Return true for '+', false for '-' and error otherwise */
	protected boolean isPosStrand(String ch) {
	
		boolean posStrand = true;
		if (ch.equals("-"))
			posStrand = false;
		else if (!ch.equals("+"))
			throw new RuntimeException("Strand has to be '+' or '-'");
		
		return posStrand;
	}
	
	

}
