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


/**
 * A position on the genome (e.g., SNP)
 */
public class GenomicElement {

	/** The ID of the element */
	public String id_ = null;	
	/** The chromosome (chr1, ..., chr22, chrX, chrY, chrM) */
	public String chr_ = null;
	/** The start position */
	public int start_ = -1;
	/** The end position */
	public int end_ = -1;
	/** The strand (true for +, false for -) */
	public boolean posStrand_ = true;

	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public GenomicElement(String id) {

		id_ = id;
	}

	
	// ----------------------------------------------------------------------------

	/** Clone this element*/
	public GenomicElement clone() {
		
		GenomicElement clone = new GenomicElement(id_);
		clone.setPosition(chr_, start_, end_, posStrand_);
		
		return clone;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Get a string representation in BED format */
	public String bedFormatString() {
		
		return chr_ + "\t" + start_ + "\t" + end_ + "\t" + id_;
	}

	
	// ----------------------------------------------------------------------------

	/** Get a window around the given element (up/down depends on posStrand_) */
	public int[] getWindow(int downstreamOffset, int upstreamOffset) {
		
		int[] window = new int[2];
		window[0] = start_;
		window[1] = end_;
		
		if (posStrand_) {
			window[0] = start_ - upstreamOffset;
			window[1] = end_ + downstreamOffset;
		} else {
			window[0] = start_ - downstreamOffset;
			window[1] = end_ + upstreamOffset;
		}

		return window;
	}


	// ----------------------------------------------------------------------------

	/** Same as getWindow(), but sets the new start and end position in this genomic element */
	public void expand(int downstreamOffset, int upstreamOffset) {
		
		int[] window = getWindow(downstreamOffset, upstreamOffset);
		start_ = window[0];
		end_ = window[1];
	}

	
	// ----------------------------------------------------------------------------

	/** Returns true if the given element overlaps this element */
	public boolean overlaps(GenomicElement element) {
		
		return chr_.equals(element.chr_) && start_ <= element.end_ && end_ >= element.start_;
	}


	// ============================================================================
	// GETTERS AND SETTERS

	public String getId() { return id_; }
	
	public void setPosition(String chr, int start) {
		chr_ = chr;
		start_ = start;
		end_ = start;
	}

	public void setPosition(String chr, int start, int end, boolean posStrand) {
		chr_ = chr;
		start_ = start;
		end_ = end;
		posStrand_ = posStrand;
	}

}
