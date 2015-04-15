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
import java.util.Map.Entry;
import java.util.TreeMap;


/**
 * A position on the genome (e.g., SNP)
 */
public class Chromosome {

	/** The genomic elements of this chromosome, sorted by starting position */
	private TreeMap<Integer, GenomicElement> elements_ = null;
	
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public Chromosome() {

		elements_ = new TreeMap<Integer, GenomicElement>();
	}

	
	// ----------------------------------------------------------------------------

	/** Add the given element to this chromosome */
	public void addElement(GenomicElement element) {
		
		elements_.put(element.start_, element);
	}

	
	// ----------------------------------------------------------------------------

	/** Get the elements in the given window */
	public ArrayList<GenomicElement> getElementsIn(int start, int end) {
				
		return new ArrayList<GenomicElement>(elements_.subMap(start, true, end, true).values());
	}

	
	// ----------------------------------------------------------------------------

	/** Get all elements */
	public Collection<GenomicElement> getElements() {
				
		return elements_.values();
	}

	
	// ----------------------------------------------------------------------------

	/** Get the nearest neighboring element (based on the start_ positions) */
	public GenomicElement getNearestElement(int pos) {

		Entry<Integer, GenomicElement> prev = elements_.floorEntry(pos);
		Entry<Integer, GenomicElement> next = elements_.ceilingEntry(pos);
		
		// The case where either one or both are null
		if (prev == null && next == null)
			return null;
		if (prev == null)
			return next.getValue();
		if (next == null)
			return prev.getValue();

		// Distance to the two neighbors
		int deltaPrev = pos - prev.getKey();
		int deltaNext = next.getKey() - pos;
		assert deltaPrev >= 0 && deltaNext >= 0;
		
		// Return the closer
		if (deltaPrev <= deltaNext)
			return prev.getValue();
		else
			return next.getValue();
	}


	// ============================================================================
	// PRIVATE METHODS
		

	// ============================================================================
	// GETTERS AND SETTERS

	public TreeMap<Integer, GenomicElement> getElementTreeMap() { return elements_; }
	
	public int getNumElements() { return elements_.size(); }
}
