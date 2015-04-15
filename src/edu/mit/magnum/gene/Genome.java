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

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;

import edu.mit.magnum.FileExport;
import edu.mit.magnum.Settings;


/**
 * A set of genomic elements organized by chromosome
 */
public class Genome {

	/** The chromosomes (chr1, ..., chr22, chrX, chrY, chrM) */
	private LinkedHashMap<String, Chromosome> chromosomes_ = null;

	/** The total number of elements */
	private int numElements_ = 0;

	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public Genome() {

		// initialize
		chromosomes_ = new LinkedHashMap<String, Chromosome>();
		numElements_ = 0;
		
		for (int i=1; i<=22; i++)
			chromosomes_.put("chr" + i, new Chromosome());
		if (!Settings.ignoreAllosomes_) {
			chromosomes_.put("chrX", new Chromosome());
			chromosomes_.put("chrY", new Chromosome());
		}
		//chromosomes_.put("chrM", new Chromosome());
	}

	
	/** Constructor */
	@SuppressWarnings("rawtypes")
	public Genome(Collection elements) {
		
		this();
		addElements(elements);
	}

	
	// ----------------------------------------------------------------------------

	/** Add elements to the genome */
	@SuppressWarnings({ "rawtypes" })
	public void addElements(Collection elements) {
		
		addElements(elements, null);
	}

	// ----------------------------------------------------------------------------

	/** Add the elements of the specified chromosome to the genome (add all elements if chromosome is null or empty) */
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public void addElements(Collection elements, String chromosome) {
		
		boolean addAll = (chromosome == null) || (chromosome.length() == 0);
		Iterator<GenomicElement> iter = elements.iterator();
		
		while (iter.hasNext()) {
			GenomicElement nextElement = iter.next();

			// Skip if not specified chromosome
			if (nextElement.chr_ == null || 
					(!addAll && !nextElement.chr_.equals(chromosome)))
				continue;
			
			// Add element to the corresponding chromosome
			chromosomes_.get(nextElement.chr_).addElement(nextElement);
		}
		
		// Count number of elements
		if (addAll)
			numElements_ = elements.size();
		else
			numElements_ = chromosomes_.get(chromosome).getNumElements();
	}
	
	
	// ----------------------------------------------------------------------------

	/** Get the elements in the given window */
	public ArrayList<GenomicElement> getElementsIn(String chr, int start, int end) {
		
		return chromosomes_.get(chr).getElementsIn(start, end);
	}

	
	// ----------------------------------------------------------------------------

	/** Get all elements of the given chromosome */
	public Collection<GenomicElement> getElements(String chr) {
		
		return chromosomes_.get(chr).getElements();
	}

	
	// ----------------------------------------------------------------------------

	/** Get the nearest element (in either direction) */
	public GenomicElement getNearestElement(String chr, int pos) {
		
		return chromosomes_.get(chr).getNearestElement(pos);
	}

	
	// ----------------------------------------------------------------------------

	/** Write a BED file with all the elements */
	public void writeBedFile(String filename) {
		
		FileExport writer = new FileExport(filename);
		
		for (Chromosome chr : chromosomes_.values()) {
			for (GenomicElement el : chr.getElements()) {
				String nextLine = el.bedFormatString();
				writer.println(nextLine);
			}
		}
		writer.close();
	}

	
	// ============================================================================
	// PRIVATE METHODS
		

	// ============================================================================
	// GETTERS AND SETTERS

	public int getNumElements() { return numElements_; }
	public int getNumChromosomes() { return chromosomes_.size(); }
	
	public HashMap<String, Chromosome> getChromosomes() { return chromosomes_; }
	
}
