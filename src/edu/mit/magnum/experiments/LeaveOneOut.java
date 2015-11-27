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
package edu.mit.magnum.experiments;

import java.util.ArrayList;

import ch.unil.gpsutils.FileExport;
import edu.mit.magnum.Magnum;
import edu.mit.magnum.net.Network;
import edu.mit.magnum.netops.GroupNetworks;
import edu.mit.magnum.netops.Union;

/**
 * For each of the networks, leave one of the samples out and see how many edges disappear.
 */
public class LeaveOneOut extends GroupNetworks {


	
	// ============================================================================
	// PUBLIC METHODS

	/** Constructor */
	public LeaveOneOut(Magnum mag) {
		super(mag, mag.set.networkDir_, mag.set.networkGroupFile_, mag.set.networkFilePrefix_);
	}
	
	
	// ----------------------------------------------------------------------------

	/** Perform the analysis */
	public void run() {
		
		// Percentage of network edges unique to the given sample
		FileExport writer = new FileExport(mag.log, "leaveOneOut-cutoff0.1.txt");
		
		// For each network
		for (String name : networkSets_.keySet()) {
			
			// The samples
			ArrayList<String> filenames = networkSets_.get(name);
			mag.log.println("- " + name + " (" + filenames.size() + " networks)\n");
			
			// Union of all
			Union union = new Union(mag, networkDir_, filenames);
			Network netAll = union.run();
			int numEdgesAll = netAll.getNumEdges();

			if (filenames.size() == 1) {
				writer.println(filenames.get(0) + "\t" + 1);
				writer.flush();
				continue;
			}

			// Leave one out
			for (String leaveOutFile : filenames) {
				// The remaining samples
				ArrayList<String> remainingFiles = new ArrayList<String>(filenames);
				remainingFiles.remove(leaveOutFile);
				
				// Union of the remaining samples
				union = new Union(mag, networkDir_, remainingFiles);
				Network leaveOneOut = union.run();
				
				// By construction, the difference are the edges contributed by the left out net
				int numEdges_i = numEdgesAll - leaveOneOut.getNumEdges();
				
				// Write result
				writer.println(leaveOutFile + "\t" + numEdges_i/(double)numEdgesAll);
				writer.flush();
			}
			
		}
		writer.close();
	}
	
	

	
	// ============================================================================
	// STATIC METHODS

	
	// ============================================================================
	// SETTERS AND GETTERS

		
}
