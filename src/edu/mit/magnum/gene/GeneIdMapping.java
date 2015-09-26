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
import java.util.HashSet;

import edu.mit.magnum.FileParser;
import edu.mit.magnum.Magnum;


/**
 * 
 */
public class GeneIdMapping {

	/** The magnum instance */
	private Magnum mag;

	/** The unique instance of the mapping (Singleton design pattern) */
	static private GeneIdMapping instance_ = null;
	
	/** Mapping ensembl to entrez */
	private HashMap<String, HashSet<String>> ensembl2entrez_ = null;
	
	
	// ============================================================================
	// PUBLIC METHODS
	
	private GeneIdMapping(Magnum mag) {
		this.mag = mag;
	}
	
	/** Get the unique instance */
	public static GeneIdMapping getInstance(Magnum mag) {
	
		if (instance_ == null)
			instance_ = new GeneIdMapping(mag);
		
		return instance_;
	}
	
	
	// ----------------------------------------------------------------------------

	/** Map ensembl to entrez ids */
	public HashSet<String> ensembl2entrez(String ensemblId) {
		
		return ensembl2entrez_.get(ensemblId);
	}

	
	// ----------------------------------------------------------------------------

	/** Load the mapping */
	public void load(String filename) {
		
		if (ensembl2entrez_ != null) {
			mag.log.warning("Gene mapping already loaded");
			return;
		}
		
		ensembl2entrez_ = new HashMap<String, HashSet<String>>();
		FileParser parser = new FileParser(mag.log, filename);
		
		while(true) {
			// Read next line
			String[] nextLine = parser.readLine();
			if (nextLine == null)
				break;
			
			// Check number of columns
			if (nextLine.length != 3)
				throw new RuntimeException("Expected three columns (ensembl id, entrez id, gene symbol)");
			
			// Parse ensembl id
			String ensg = nextLine[0];
			if (!(ensg.length() > 4 && ensg.substring(0, 4).equals("ENSG")))
				throw new RuntimeException("Invalid ENSEMBL gene ID (expected 'ENSG...'): " + ensg);
			ensg = removeEnsemblVersion(ensg);

			// Parse entrez id
			String entrez = nextLine[1];
			if (entrez.length() > 0) {
				try {
					Integer.valueOf(entrez);
				} catch (NumberFormatException e) {
					throw new RuntimeException("Invalid Entrez gene ID (expected an integer number): " + entrez);
				}
			}
			
			// Parse gene symbol
			//String symbol = nextLine[2];
			
			HashSet<String> entrezSet = ensembl2entrez_.get(ensg);
			if (entrezSet == null) {
				entrezSet = new HashSet<String>(1);
				ensembl2entrez_.put(ensg, entrezSet);
			}
			entrezSet.add(entrez);
		}
	}

	
	// ----------------------------------------------------------------------------

	/** If this is an ensembl id, remove the version number */
	public String removeEnsemblVersion(String id) {
		
		if (id.length() > 4 && id.substring(0, 4).equals("ENSG")) {
			int dot = id.lastIndexOf(".");
			if (dot != -1)
				id = id.substring(0, dot);
		}
		return id;
	}

	
	
	// ============================================================================
	// PRIVATE METHODS
	


	
	// ============================================================================
	// GETTERS AND SETTERS
		
	
}
