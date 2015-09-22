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

import java.util.HashMap;
import java.util.LinkedHashMap;

import edu.mit.magnum.*;


/**
 * Custom gene annotation (bed file)
 */
public class GeneAnnotationCustom extends GeneAnnotation {

	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public GeneAnnotationCustom(String chromosomeToBeLoaded) {
		
		super(Magnum.set.geneCoordFile_, chromosomeToBeLoaded, false);
	}
	

	// ----------------------------------------------------------------------------

	/** Load gene coordinates for the given gene set (load all genes if the set is empty) */
	public HashMap<String, Gene> loadAnnotation() {
				
		genes_ = new LinkedHashMap<String, Gene>();
		
		// Open the file
		FileParser parser;
		if (annotationFile_ == null)
			parser = new FileParser(MagnumSettings.class.getClassLoader().getResourceAsStream("edu/mit/magnum/gene/rsc/gene_coord.bed"));
		else
			parser = new FileParser(annotationFile_);
						
		while (true) {
			// Read next line
			String[] nextLine = parser.readLine();
			if (nextLine == null)
				break;

			// Check number of columns
			if (nextLine.length != 6)
				parser.error("Expected 6 columns");
			
			// Chromosome
			String chr = nextLine[0];
			// Continue if not the specified chromosome
			if (chromosomeToBeLoaded_ != null && chromosomeToBeLoaded_.length() > 0 && !chr.equals(chromosomeToBeLoaded_))
				continue;

			// Gene id
			String gene_id = nextLine[3];
			
			// If a gene set to be loaded was specified and this gene is NOT in this set, continue
			if (genesToBeLoaded_ != null) {
				if (genesToBeLoaded_.containsKey(gene_id))
					// Flag this gene as found
					genesToBeLoaded_.put(gene_id, true);
				else
					continue;
			}
					
			// Check that the id is unique
			if (genes_.containsKey(gene_id))
				parser.error("Duplicate gene id: " + gene_id);
			
			// Create the gene
			Gene nextGene = new Gene(gene_id);
			genes_.put(gene_id, nextGene);
				
			// Position
			int start = Integer.parseInt(nextLine[1]) + 1;
			int end = Integer.parseInt(nextLine[2]);
			boolean posStrand = isPosStrand(nextLine[5]);	
			nextGene.setPosition(chr, start, end, posStrand);

		}		
		parser.close();		

		return genes_;
	}
	
	
	// ============================================================================
	// PRIVATE METHODS

	
}
