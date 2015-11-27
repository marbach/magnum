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
package ch.unil.gpsutils;

import java.io.PrintWriter;
import java.io.StringWriter;

/**
 * Logger -- can be extended to customize/redirect outputs
 */
public class Logger {

	/** Verbose output */
	private boolean verbose;
	
	// ============================================================================
	// PUBLIC FUNCTIONS
	
	/** Write empty line to stdout */
	public void println() {
		print("\n");
	}

	/** Write line to stdout */
	public void println(String msg) {
		print(msg + "\n");
	}

	/** Write line to stdout */
	public void printlnVerbose(String msg) {
		if (verbose)
			println(msg);
	}

	/** Write line to stdout */
	public void printlnVerbose(String verboseMsg, String nonVerboseMsg) {
		if (verbose)
			println(verboseMsg);
		else
			println(nonVerboseMsg);
	}

	/** Print warning message */
	public void warning(String msg) {
		print("WARNING: " + msg + "\n");
	}

	/** Write string to stdout */
	public void print(String msg) {
		System.out.print(msg);
	}
	
	/** Get stack trace as string */
	public void printStackTrace(Throwable e) {
		
		StringWriter sw = new StringWriter();
		PrintWriter pw = new PrintWriter(sw);
		e.printStackTrace(pw);
		println(sw.toString());
	}

	
	// ============================================================================
	// SETTERS AND GETTERS

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}
	public boolean getVerbose() {
		return verbose;
	}
}
