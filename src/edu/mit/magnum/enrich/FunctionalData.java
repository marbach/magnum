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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import edu.mit.magnum.*;
import edu.mit.magnum.gene.*;
import edu.mit.magnum.net.Network;


/**
 * The functional data available for the genes (e.g., network kernels)
 */
public class FunctionalData {

	/** The original data matrix TODO delete and clean memory */
	private DoubleMatrix2D unfilteredData;
	/** The rows of the original data matrix */
	private ArrayList<String> unfilteredDataRows;
	/** The cols of the original data matrix */
	private ArrayList<String> unfilteredDataCols;
	
	/** The data matrix (genes in rows) */
	private DoubleMatrix2D data_ = null;
	/** The number of genes (rows) */
	private int numGenes_ = -1;
	/** Genes (rows of the data matrix) */
	private LinkedHashMap<String, Integer> genes_ = null;
	/** Column names of the data matrix */
	private ArrayList<String> colNames_ = null;
	/** Indicates whether the data matrix is a kernel / similarity matrix */
	private boolean isPairwiseData_ = false;

	/** The number of genes in the functional data file */
	//private int unfilteredDataNumGenes_ = -1;
	/** The indexes of the columns that were loaded in data_ */
	private ArrayList<Integer> functDataColIndexes_ = null;

	/** Genes that were not loaded from the functional data because they are not in geneScores */
	private ArrayList<String> genesMissingScores_ = null;

	
	// ============================================================================
	// PUBLIC METHODS

	/** Constructor */
	public FunctionalData(Network network, DoubleMatrix2D kernel,
			File excludedGenePairsFile,	ArrayList<Gene> geneScores) {
		
		this.unfilteredData = kernel;
		
		// Get the row/col names of the kernel matrix from the network
		unfilteredDataRows = new ArrayList<String>(network.getNumNodes());
		for (int i=0; i<network.getNumNodes(); i++)
			unfilteredDataRows.add(network.getNode(i).getId());
		unfilteredDataCols = unfilteredDataRows;
		
		initialize(excludedGenePairsFile, geneScores);
	}
	
	/** Constructor */
	public FunctionalData(File functionalDataFile, File excludedGenePairsFile,
			ArrayList<Integer> functionalDataCols, ArrayList<Gene> geneScores) {

		// Load the kernel from the file
		functDataColIndexes_ = functionalDataCols;
		loadUnfilteredData(functionalDataFile);			
		initialize(excludedGenePairsFile, geneScores);
	}

	// ============================================================================
	// PRIVATE METHODS

	/** Load genes and their properties */
	private void initialize(File excludedGenePairsFile, ArrayList<Gene> geneScores) {

		// A hashmap with all genes that have scores
		HashMap<String, Gene> genesWithScores = new HashMap<String, Gene>();
		for (Gene g : geneScores)
			genesWithScores.put(g.id_, g);

		// Initialize genes_ with the set of overlapping genes between the
		// functional data and the gene scores
		initializeGenes(genesWithScores);
		// Load the data for the overlapping genes
		initializeData();
		// Normalize by row/col sums to adjust for hubs
		if (Magnum.set.scaleKernel_)
			scaleKernel();
		// Load the gene pairs that should be excluded from enrichment analysis, set corresponding data entries to NaN
		loadExcludedGenePairs(excludedGenePairsFile);
	}

//	/** Load genes and their properties */
//	private void load(File functionalDataFile, File excludedGenePairsFile, ArrayList<Gene> geneScores) {
//
//		// A hashmap with all genes that have scores
//		HashMap<String, Gene> genesWithScores = new HashMap<String, Gene>();
//		for (Gene g : geneScores)
//			genesWithScores.put(g.id_, g);
//
//		// Initialize genes_ with the set of overlapping genes between the
//		// functional data and the gene scores
//		loadGenes(functionalDataFile, genesWithScores);
//		// Load the data for the overlapping genes
//		loadData(functionalDataFile);
//		// Normalize by row/col sums to adjust for hubs
//		if (Magnum.set.scaleKernel_)
//			scaleKernel();
//		// Load the gene pairs that should be excluded from enrichment analysis, set corresponding data entries to NaN
//		loadExcludedGenePairs(excludedGenePairsFile);
//	}

	// ----------------------------------------------------------------------------

	/**
	 * Initialize genes_ with the set of overlapping genes between the
	 * functional data and the gene scores
	 */
	private void initializeGenes(HashMap<String, Gene> genesWithScores) {

		// A hashmap with the overlapping genes and their index, in the order in
		// which they occur in the funct data file
		genes_ = new LinkedHashMap<String, Integer>();
		ArrayList<Integer> colIndexes = new ArrayList<Integer>();
		colNames_ = new ArrayList<String>();
		genesMissingScores_ = new ArrayList<String>();
		int count = 0;
		
		for (int i=0; i<unfilteredDataRows.size(); i++) {
			String id = unfilteredDataRows.get(i);

			// Add this gene if it has a score			
			if (genesWithScores.containsKey(id)) {
				// Check that the gene doesn't exist already
				if (genes_.containsKey(id))
					throw new RuntimeException("Gene '" + id + "' is listed twice");

				genes_.put(id, count++);
				colIndexes.add(i);
				colNames_.add(id);
				
			} else {
				genesMissingScores_.add(id);
			}
		}

		numGenes_ = genes_.size();
		assert colIndexes.size() == numGenes_;

		// Check if it's a gene x gene matrix (kernel)
		if (unfilteredData.columns() == unfilteredData.rows()) {
			isPairwiseData_ = true;
			functDataColIndexes_ = colIndexes;
		} else {
			isPairwiseData_ = false;
		}
	}

//	/**
//	 * Initialize genes_ with the set of overlapping genes between the
//	 * functional data and the gene scores
//	 */
//	private void loadGenes(File functionalDataFile, HashMap<String, Gene> genesWithScores) {
//
//		// A hashmap with the overlapping genes and their index, in the order in
//		// which they occur in the funct data file
//		genes_ = new LinkedHashMap<String, Integer>();
//		ArrayList<Integer> colIndexes = new ArrayList<Integer>();
//
//		// Open the functional data file
//		FileParser parser = new FileParser(functionalDataFile);
//		// Read the header
//		int numCols = parser.readLine().length;
//		int count = 0;
//
//		while (true) {
//			// Read next line
//			String[] nextLine = parser.readLine();
//			if (nextLine == null)
//				break;
//
//			// Check number of columns
//			if (nextLine.length != numCols)
//				throw new RuntimeException("Line " + parser.getLineCounter() + " has " + nextLine.length + " columns (not same as header)");
//
//			// Add this gene if it has a score
//			String id = nextLine[0];
//			if (genesWithScores.containsKey(id)) {
//				// Check that the gene doesn't exist already
//				if (genes_.containsKey(id))
//					throw new RuntimeException("Gene '" + id + "' is listed twice");
//
//				genes_.put(id, count++);
//				colIndexes.add(parser.getLineCounter() - 1);
//			}
//		}
//		parser.close();
//
//		numGenes_ = genes_.size();
//		assert colIndexes.size() == numGenes_;
//		functDataNumGenes_ = parser.getLineCounter() - 2;
//
//		// Check if it's a gene x gene matrix (kernel)
//		if (numCols - 1 == functDataNumGenes_) {
//			isPairwiseData_ = true;
//			functDataColIndexes_ = colIndexes;
//		} else {
//			isPairwiseData_ = false;
//		}
//	}

	// ----------------------------------------------------------------------------

	/**
	 * Load the gene pairs that should be excluded from enrichment analysis, set
	 * corresponding data entries to NaN
	 */
	private void loadExcludedGenePairs(File excludedGenePairsFile) {

		if (!isPairwiseData_ || excludedGenePairsFile == null
				|| excludedGenePairsFile.length() == 0)
			return;

		GeneIdMapping idMapping = GeneIdMapping.getInstance();
		boolean translateToEntrez = Magnum.set.idTypeFunctionalData_.equalsIgnoreCase("entrez");

		// Open the file
		FileParser parser = new FileParser(excludedGenePairsFile);
		String[] header = parser.readLine();

		// Find the columns corresponding to the two gene ids
		int colGene1 = -1;
		int colGene2 = -1;
		for (int i = 0; i < header.length; i++) {
			if (header[i].equals("gene1_id"))
				colGene1 = i;
			else if (header[i].equals("gene2_id"))
				colGene2 = i;
		}

		// Check that the two columns are present
		if (colGene1 == -1)
			Magnum.log.error("Did not find mandatory column 'gene1_id'");
		if (colGene2 == -1)
			Magnum.log.error("Did not find mandatory column 'gene2_id'");

		int numExcluded = 0;

		while (true) {
			// Read next line
			String[] nextLine = parser.readLine();
			if (nextLine == null)
				break;

			// The two gene ids
			String id1 = nextLine[colGene1];
			String id2 = nextLine[colGene2];
			if (Magnum.set.idTypeFunctionalData_.equals("ensembl")) {
				id1 = idMapping.removeEnsemblVersion(id1);
				id2 = idMapping.removeEnsemblVersion(id2);
			}

			// All mapped ids (synonyms)
			HashSet<String> set1;
			HashSet<String> set2;

			if (translateToEntrez) {
				set1 = idMapping.ensembl2entrez(id1);
				set2 = idMapping.ensembl2entrez(id2);

			} else {
				set1 = new HashSet<String>(1);
				set2 = new HashSet<String>(1);
				set1.add(id1);
				set2.add(id2);
			}

			for (String mappedId1 : set1) {
				Integer index1 = genes_.get(mappedId1);
				if (index1 == null)
					continue;

				for (String mappedId2 : set2) {
					Integer index2 = genes_.get(mappedId2);
					if (index2 == null)
						continue;
					
					if ((index1 == 3979 && index2 == 9831) ||
							(index2 == 3979 && index1 == 9831))
						Magnum.log.println();

					// Set entry to NaN
					data_.set(index1, index2, Double.NaN);
					data_.set(index2, index1, Double.NaN);
					numExcluded++;
					
				}
			}
		}
		parser.close();
		Magnum.log.println("- " + numExcluded + " gene pairs excluded");
	}

	
	// ----------------------------------------------------------------------------

	/** 
	 * Exclude neighboring genes based on distance, set corresponding data entries to NaN.
	 * Note, genes have previously been expanded based on window size and neighborhood distance. 
	 */
	public void excludeNeighbors(ArrayList<Gene> genesScoreList) {
		
		for (int i=0; i<genesScoreList.size(); i++) {
			Gene g_i = genesScoreList.get(i);
			int index1 = genes_.get(g_i.id_);
			
			for (int j=i+1; j<genesScoreList.size(); j++) {
				Gene g_j = genesScoreList.get(j);
				assert !g_i.equals(g_j);
				
				if (g_i.overlaps(g_j)) {
					int index2 = genes_.get(g_j.id_);
					
					// Set entry to NaN
					data_.set(index1, index2, Double.NaN);
					data_.set(index2, index1, Double.NaN);
				}
			}
		}
	}

		
	// ----------------------------------------------------------------------------

	/**
	 * 
	 */
	private void initializeData() {

		// Read header -- this initializes colNames_, which are relevant for data that is not pairwise
		//parseGenePropertiesHeader(parser.readLine());

		// Initialize matrix
		if (isPairwiseData_) {
			data_ = new DenseDoubleMatrix2D(numGenes_, numGenes_);
		} else {
			throw new RuntimeException("Initialization of not pairwise FunctionalData not implemented without loading from file");
			//data_ = new DenseDoubleMatrix2D(numGenes_, colNames_.size());
		}

		assert genes_.size() == numGenes_;
		assert functDataColIndexes_.size() == numGenes_;
		assert colNames_.size() == numGenes_;
		assert data_.rows() == numGenes_;
		assert data_.columns() == numGenes_;
		
		for (int i=0; i<numGenes_; i++) {
			String id = colNames_.get(i);
			int kernelIndex = functDataColIndexes_.get(i);
			assert genes_.get(id) == i;
			
			// Copy the row to the filtered kernel
			for (int j=0; j<numGenes_; j++)
				// TODO make it quick
				data_.set(i, j, unfilteredData.get(kernelIndex, functDataColIndexes_.get(j)));
		}
	}

//	/**
//	 * Initialize genes_ with the set of overlapping genes between the
//	 * functional data and the gene scores
//	 */
//	private void loadData(File functionalDataFile) {
//
//		// Open the file
//		FileParser parser = new FileParser(functionalDataFile);
//		// Read header
//		parseGenePropertiesHeader(parser.readLine());
//
//		// Initialize matrix
//		if (isPairwiseData_)
//			data_ = new DenseDoubleMatrix2D(numGenes_, numGenes_);
//		else
//			data_ = new DenseDoubleMatrix2D(numGenes_, colNames_.size());
//
//		genesMissingScores_ = new ArrayList<String>();
//
//		while (true) {
//			// Read next line
//			String[] nextLine = parser.readLine();
//			if (nextLine == null)
//				break;
//
//			// The gene and its index
//			String id = nextLine[0];
//			Integer index = genes_.get(id);
//			// If the gene is not part of the gene score data, skip
//			if (index == null) {
//				genesMissingScores_.add(id);
//				continue;
//			}
//
//			// Parse properties
//			for (int i = 0; i < functDataColIndexes_.size(); i++)
//				data_.set(index, i, Double
//						.parseDouble(nextLine[functDataColIndexes_.get(i)]));
//		}
//		parser.close();
//	}

	
	// ----------------------------------------------------------------------------

	/** Parse the header of a gene property file, return number of columns */
	private void loadUnfilteredData(File functionalDataFile) {
		
		// Count the lines
		try {
			BufferedReader reader = new BufferedReader(new FileReader(functionalDataFile));
			int lines = 0;
			while (reader.readLine() != null) 
				lines++;

			reader.close();
			numGenes_ =  lines - 1;
		} catch (IOException e) {
			throw new RuntimeException(e);
		} 

		// Open the file
		FileParser parser = new FileParser(functionalDataFile);
		// Read header (sets unfilteredDataCols)
		parseGenePropertiesHeader(parser.readLine());

		// Initialize matrix
		unfilteredData = new DenseDoubleMatrix2D(numGenes_, unfilteredDataCols.size());
		unfilteredDataRows = new ArrayList<String>();
		
		int i = 0;
		while (true) {
			// Read next line
			String[] nextLine = parser.readLine();
			if (nextLine == null)
				break;

			// Add the gene id
			unfilteredDataRows.add(nextLine[0]);
			// Parse properties
			for (int j=1; j<unfilteredDataCols.size(); j++)
				// TODO make quick
				unfilteredData.set(i, j-1, Double.parseDouble(nextLine[j]));
			i++;
		}
		parser.close();

	}

	
	// ----------------------------------------------------------------------------

	/** Parse the header of a gene property file, return number of columns */
	private void parseGenePropertiesHeader(String[] header) {

		if (header.length < 2)
			throw new RuntimeException("Gene property files must have at least two tab-separated columns");

		// Is the data matrix square?
		isPairwiseData_ = (numGenes_ == header.length - 1);
		
		// Test if this is a header by trying to parse the second column as number, which should fail
		if (!isPairwiseData_) {
			boolean hasHeader = false;
			try {
				Double.valueOf(header[1]);
			} catch (NumberFormatException e) {
				hasHeader = true;
			}
			if (!hasHeader)
				throw new RuntimeException("File has no header (column names cannot be numbers)");
		}

		// Parse the gene property names
		unfilteredDataCols = new ArrayList<String>(header.length);
		for (int i = 1; i < header.length; i++)
			unfilteredDataCols.add(header[i]);
	}

	// ----------------------------------------------------------------------------

	/**
	 * Normalize by row/col sums to adjust for hubs. K'(i,j) =
	 * K(i,j)/sqrt(rowSums(K)[i] * colSums(K)[j])
	 */
	private void scaleKernel() {

		if (!isPairwiseData_)
			return;

		Magnum.log.println("- Scaling kernel (adjusting for hubs)...");

		// Note, I also did this using Blas following a stackoverflow comment
		// (multiply matrix by a vector of ones), but it was much slower
		// rowsums are equal to colsums because kernel is sysmmetric
		double[] rowsums = new double[numGenes_];
		for (int i = 0; i < numGenes_; i++) {
			double sum = 0;
			for (int j = 0; j < numGenes_; j++)
				sum += data_.get(i, j);
			rowsums[i] = sum;
		}

		for (int i = 0; i < numGenes_; i++)
			for (int j = 0; j < numGenes_; j++)
				data_.set(i, j,
						data_.get(i, j) / Math.sqrt(rowsums[i] * rowsums[j]));
	}

	// ============================================================================
	// GETTERS AND SETTERS

	public DoubleMatrix2D getData() {
		return data_;
	}

	public double get(int i, int j) {
		return data_.get(i, j);
	}

	public HashMap<String, Integer> getGenes() {
		return genes_;
	}

	public int getNumGenes() {
		return numGenes_;
	}

	public ArrayList<String> getColNames() {
		return colNames_;
	}

	public boolean getIsPairwiseData() {
		return isPairwiseData_;
	}

	public ArrayList<String> getGenesMissingScores() {
		return genesMissingScores_;
	}
}
