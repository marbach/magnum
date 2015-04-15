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

import java.util.ArrayList;
import java.util.Iterator;
import java.util.TreeSet;



/**
 * 
 */
public class Percentile {
	
	/** The points of all curves at the given position (k) */
	private TreeSet<Point> points_ = null;
		
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	public Percentile() {

		points_ = new TreeSet<Point>();
	}

	
	// ----------------------------------------------------------------------------
	
	/** Add a point */
	public void addPoint(Point point) {
		
		points_.add(point);
	}
	
	
	// ----------------------------------------------------------------------------

	/** Get the values corresponding to the elements at the given offsets from both ends of the list */
	public ArrayList<Double> getValues(ArrayList<Integer> offsets) {

		// First, do always the median
		ArrayList<Integer> midPoint = new ArrayList<Integer>();
		midPoint.add((int) Math.round(0.5*points_.size()));
		ArrayList<Double> median = getValues(midPoint, true);
		assert median.size() == 1;
		
		ArrayList<Double> values = new ArrayList<Double>();
		values.add(median.get(0));
		
		ArrayList<Double> valuesForward = getValues(offsets, true);
		for (int i=0; i<valuesForward.size(); i++)
			values.add(valuesForward.get(i));
		
		ArrayList<Double> valuesBackward = getValues(offsets, false);
		for (int i=valuesBackward.size()-1; i>=0; i--)
			values.add(valuesBackward.get(i));
		
		return values;
	}

	
	// ----------------------------------------------------------------------------

	/** Get the empirical p-value for the given point */
	public double pValue(double x) {
		
		int numSmaller = 0;
		int numTies = 0;
		double epsilon = 1e-12;
		
		for (Point point : points_) {
			// the current point is smaller
			if (point.value_ + epsilon < x)
				numSmaller++;
			// the current point is equal
			else if (Math.abs(x - point.value_) < epsilon)
				numTies++;
			// the current point is greater
			else
				break;
		}
		return (numSmaller + (numTies/2.0)) / (double)points_.size();
	}

		
	// ============================================================================
	// PRIVATE METHODS
		
	/** 
	 * Get the values corresponding to the elements at the given indexes
	 * (equivalent to ArrayList.get()). Forward indicates whether the index
	 * is counted from the start or end of the list.
	 */
	private ArrayList<Double> getValues(ArrayList<Integer> indexes, boolean forward) {
		
		ArrayList<Double> values = new ArrayList<Double>();
		
		Iterator<Point> iter = (forward ? points_.iterator() : points_.descendingIterator());
		int k = 0;
		
		// For each index
		for (int i=0; i<indexes.size(); i++) {
			int index_i = indexes.get(i);
			assert index_i > k;
			
			// Get the point at position index_i
			Point p = null;
			do {
				p = iter.next();
				k++;
			} while (k < index_i);
			
			// Add the value
			values.add(p.value_);
		}
		return values;
	}

	
	
}
