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

import ch.unil.gpsutils.FileParser;
import edu.mit.magnum.*;


/**
 * Gencode gene annotation.
 * Gene types are (attribute 'gene_type'): 
 *    antisense
 *    snRNA
 *    miRNA
 *    snoRNA
 *    polymorphic_pseudogene
 *    lincRNA
 *    protein_coding
 *    3prime_overlapping_ncrna
 *    misc_RNA
 *    rRNA
 *    sense_intronic
 *    processed_transcript
 *    pseudogene
 *    sense_overlapping
 */
public class GeneAnnotationGencode extends GeneAnnotation {

	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public GeneAnnotationGencode(Magnum mag, String chromosomeToBeLoaded, boolean loadProteinCodingOnly) {
		
		super(mag, mag.set.gencodeAnnotationFile_, chromosomeToBeLoaded, loadProteinCodingOnly);
	}
	

	// ----------------------------------------------------------------------------

	/** Load gene coordinates for the given gene set (load all genes if the set is empty) */
	public HashMap<String, Gene> loadAnnotation() {
				
		genes_ = new LinkedHashMap<String, Gene>();
		// To print gene types
		//HashSet<String> geneType = new HashSet<String>();
		
		// Open the file
		FileParser parser = new FileParser(mag.log, annotationFile_);
		GeneIdMapping mapping = GeneIdMapping.getInstance(mag.log);
		
		// Skip the first 5 lines (start with #)
		String[] nextLine = parser.readLine();
		while (nextLine[0].startsWith("#"))
			nextLine = parser.readLine();
				
		while (nextLine != null) {
			// Check number of columns
			if (nextLine.length != 9)
				parser.error("Expected 9 columns");
			
			// Check that this is a gene
			if (!nextLine[2].equals("gene"))
				parser.error("Third column expected to be 'gene'");

			// Chromosome
			String chr = nextLine[0];
			// Continue if not the specified chromosome
			if (chromosomeToBeLoaded_ != null && chromosomeToBeLoaded_.length() > 0 && !chr.equals(chromosomeToBeLoaded_)) {
				nextLine = parser.readLine();
				continue;
			}

			if (loadOnlyProteinCoding_) {
				// Check that it's a protein coding gene
				String gene_type = getGencodeKeyValue(nextLine[8], "gene_type");
				//geneType.add(gene_type);
				
				if (!gene_type.equalsIgnoreCase("protein_coding")) {
					nextLine = parser.readLine();
					continue;
				}
				
			}
			// Gene id
			String gene_id = mapping.removeEnsemblVersion(getGencodeKeyValue(nextLine[8], "gene_id").toUpperCase());
			String gene_name = getGencodeKeyValue(nextLine[8], "gene_name").toUpperCase();
			
			// If a gene set to be loaded was specified and this gene is NOT in this set, continue
			if (genesToBeLoaded_ != null) {
				// Check if this gene is part of the specified gene set
				String specifiedGene = null;
				if (genesToBeLoaded_.containsKey(gene_id))
					specifiedGene = gene_id;
				else if (genesToBeLoaded_.containsKey(gene_name))
					specifiedGene = gene_name;

				// Flag this gene as found
				if (specifiedGene != null) {
					genesToBeLoaded_.put(specifiedGene, true);
				} else {
					nextLine = parser.readLine();
					continue;
				}
			}
					
			// Check that the id is unique
			if (genes_.containsKey(gene_id))
				parser.error("Duplicate gene id: " + gene_id);
			
			// Create the gene
			Gene nextGene = new Gene(gene_id, gene_name);
			genes_.put(gene_id, nextGene);
				
			// Position
			int start = Integer.parseInt(nextLine[3]);
			int end = Integer.parseInt(nextLine[4]);
			boolean posStrand = isPosStrand(nextLine[6]);	
			nextGene.setPosition(chr, start, end, posStrand);

			// Read next line
			nextLine = parser.readLine();
		}		
		parser.close();		

		// Print gene types
		//for (String type : geneType)
		//Ngsea.println(type);
	
		return genes_;
	}
	
	
	// ============================================================================
	// PRIVATE METHODS

	/** Get the value of the given key, throw exception if not found */
	private String getGencodeKeyValue(String keyValueList, String key) {
		
		int start = keyValueList.indexOf(key + " \"");
		if (start == -1)
			throw new RuntimeException("Key not found: '" + key + "\"");
		
		start = start + key.length() + 2;
		int end = keyValueList.indexOf("\"", start);
		
		String value = keyValueList.substring(start, end);
		return value;
	}
	
}
