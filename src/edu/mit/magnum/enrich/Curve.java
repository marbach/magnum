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


/**
 * 
 */
public class Curve {

	/** The enrichment values of the curve / the y-axis */
	private ArrayList<Point> points_ = null;
	
	
	// ============================================================================
	// PUBLIC METHODS
	
	/** Constructor */
	Curve(int initialCapacity) {
		
		points_ = new ArrayList<Point>(initialCapacity);
	}

	// ----------------------------------------------------------------------------

	/** Add a point with the given value */
	public void addPoint(double value) {
		
		points_.add(new Point(this, value));
	}
	
	
	// ============================================================================
	// PRIVATE METHODS

	
	// ============================================================================
	// GETTERS AND SETTERS
	
	public int getNumPoints() { return points_.size(); }
	public ArrayList<Point> getPoints() { return points_; }
	
	public Point getPoint(int k) { return points_.get(k); }
	public double getValue(int k) { return points_.get(k).value_; }
	
}
