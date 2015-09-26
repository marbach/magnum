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
package edu.mit.magnum;

import java.io.File;
import java.util.ArrayList;
import java.util.Properties;


/** 
 * Abstract class offering functionality for the settings classes
 */
public class Settings {	
	
	/** Magnum instance */
	protected Magnum mag;
	/** The properties (settings file) */
	protected Properties prop = null;
	
	// ============================================================================
	// PUBLIC METHODS

	public Settings(Magnum mag) {
		this.mag = mag;
	}
	
	
	// ============================================================================
	// PROTECTED METHODS

	/** Get the string value of a parameter from the setting file */
	protected String getSetting(String param) {
		
		String value = prop.getProperty(param);
		if (value == null)
			throw new RuntimeException("Parameter not found in setting file: " + param);
		
		return value; 
	}

	
	// ----------------------------------------------------------------------------

	/** Get the integer value of a parameter from the setting file */
	protected int getSettingInt(String param) {
		return Integer.valueOf(getSetting(param)); 
	}

	/** Get the double value of a parameter from the setting file */
	protected double getSettingDouble(String param) {
		return Double.valueOf(getSetting(param)); 
	}

	
	// ----------------------------------------------------------------------------

	/** Parse a boolean property */
	protected boolean getSettingBoolean(String name) {
		
		String value = getSetting(name);
		if (value.equals("1") || value.equalsIgnoreCase("true") || value.equalsIgnoreCase("t"))
			return true;
		else if (value.equals("0") || value.equalsIgnoreCase("false") || value.equalsIgnoreCase("f"))
			return false;
		else
			throw new IllegalArgumentException("Parse error for boolean parameter '" + name + "': expected '1' or '0', found '" + value + "'");
	}

	
	// ----------------------------------------------------------------------------

	/** Parse an int array property */
	protected ArrayList<Integer> getSettingIntArray(String name, boolean positiveSorted) {
		
		String[] propStr = getSetting(name).split(",");
		ArrayList<Integer> prop = new ArrayList<Integer>();
		
		if (propStr.length == 1 && propStr[0].compareTo("") == 0)
			return prop;

		for (int i=0; i<propStr.length; i++)
			prop.add(Integer.valueOf(propStr[i]));
			
		if (positiveSorted && !mag.utils.posIntIncreasing(prop))
			throw new RuntimeException("Error parsing settings file, " + name + " has to be an ordered list of positive integers, given in increasing order");
		
		return prop;
	}


	// ----------------------------------------------------------------------------

	/** Parse an int array property */
	protected ArrayList<Double> getSettingDoubleArray(String name, boolean positiveSorted) {
		
		String[] propStr = getSetting(name).split(",");
		ArrayList<Double> prop = new ArrayList<Double>();
		
		if (propStr.length == 1 && propStr[0].compareTo("") == 0)
			return prop;

		for (int i=0; i<propStr.length; i++)
			prop.add(Double.valueOf(propStr[i]));
			
		if (positiveSorted && !mag.utils.posDoubleIncreasing(prop))
			throw new RuntimeException("Error parsing settings file, " + name + " has to be an ordered list of positive numbers, given in increasing order");
		
		return prop;
	}

	
	// ----------------------------------------------------------------------------

	/** Parse a string array property */
	protected String[] getStringArraySetting(Properties set, String name) {
		
		return set.getProperty(name).split(",");
	}
	
	
	// ----------------------------------------------------------------------------
    
    /** Get a file / directory saved as a string, throw exception if the name is empty */
    protected File getFileSetting(String param) {
    	
    	String filename = getSetting(param);
    	if (filename.isEmpty() || filename.equals(" "))
    		throw new RuntimeException(param + ": file/directory name is empty or has trailing whitespace");
    	
		return new File(filename);
    }


}
