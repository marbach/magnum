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

import java.util.LinkedHashMap;

import edu.mit.magnum.FileParser;
import edu.mit.magnum.Magnum;
import edu.mit.magnum.MagnumSettings;


/**
 * UCSC genome browser annotation with Entrez IDs
 * => There are inconsistent entries in the file (same entrez ids, different strand / location) 
 */
public class GeneAnnotationUcsc extends GeneAnnotation {

	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public GeneAnnotationUcsc(String chromosomeToBeLoaded, boolean loadProteinCodingOnly) {
		
		super(MagnumSettings.gencodeAnnotationFile_, chromosomeToBeLoaded, loadProteinCodingOnly);
	}

	
	// ----------------------------------------------------------------------------

	/** 
	 * Load gene coordinates for the given gene set (load all genes if the set is empty) 
	 * => There are inconsistent entries in the file (same entrez ids, different strand / location) 
	 */
	public LinkedHashMap<String, Gene> loadAnnotation() {
		
		Magnum.log.error("TBD: Check, debug, adapt so that genesToBeLoaded_ is used");
		
		genes_ = new LinkedHashMap<String, Gene>();
		FileParser parser = new FileParser(MagnumSettings.ucscAnnotationFile_);
		
		// Skip the header lines (start with #)
		String[] nextLine = parser.readLine();
		while (nextLine[0].startsWith("#"))
			nextLine = parser.readLine();

		while (nextLine != null) {
			// Check number of columns
			if (nextLine.length != 8)
				parser.error("Expected 8 columns");

			// Chromosome
			String chr = nextLine[1];
			// Continue if not the specified chromosome
			if (!chr.equals(chromosomeToBeLoaded_)) {
				nextLine = parser.readLine();
				continue;
			}
			
			// Gene id and name
			String gene_id = nextLine[7];
			String gene_name = nextLine[5];

			// Skip entries without entrez id
			if (gene_id.equals("n/a")) {
				nextLine = parser.readLine();
				continue;				
			}
			
			// Strand
			boolean posStrand = isPosStrand(nextLine[6]);
			// Start and end
			int start = Integer.parseInt(nextLine[3]);
			int end = Integer.parseInt(nextLine[4]);

			Gene nextGene = genes_.get(gene_id);
			// Create new gene
			if (nextGene == null) {
				nextGene = new Gene(gene_id, gene_name);
				nextGene.setPosition(chr, start, end, posStrand);
				genes_.put(gene_id, nextGene);
			
			// Update existing gene
			} else {
				// Check that strand and symbol are consistent
				if (nextGene.posStrand_ != posStrand || !nextGene.symbol_.equals(gene_name))
					parser.error("Inconsistent entry (strand or gene symbol doesn't agree)");
				// Update start and end
				if (start < nextGene.start_)
					nextGene.start_ = start;
				if (end > nextGene.end_)
					nextGene.end_ = end;
			}
				
			// Read next line
			nextLine = parser.readLine();
		}
		parser.close();		

		return genes_;
	}
	
}
