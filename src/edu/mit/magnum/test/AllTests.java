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

package edu.mit.magnum.test;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

import edu.mit.magnum.enrich.test.EnrichMainTest;
import edu.mit.magnum.net.test.*;
import edu.mit.magnum.netops.test.PairwiseSumTest;
import edu.mit.magnum.netops.test.UnionTest;
import edu.mit.magnum.netprop.test.*;

@RunWith(Suite.class)
//@SuiteClasses({ NetworkTest.class, AnalyzerBasicPropertiesTest.class, AnalyzerShortestPathsTest.class, 
//	AnalyzerPstepKernelTest.class, HandlerTest.class, LinkCommunityTest.class })
@SuiteClasses({ 
	NetworkTest.class, 
	BasicPropertiesTest.class, 
	ShortestPathsTest.class, 
	PstepKernelTest.class,
	TanimotoCoefficientTest.class,
	UnionTest.class,
	PairwiseSumTest.class,
	EnrichMainTest.class
	})
public class AllTests {

}
