classdef Test_C11GlucoseKinetics < matlab.unittest.TestCase
	%% TEST_C11GLUCOSEKINETICS 

	%  Usage:  >> results = run(mlglucose_unittest.Test_C11GlucoseKinetics)
 	%          >> result  = run(mlglucose_unittest.Test_C11GlucoseKinetics, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 04-Dec-2017 11:15:32 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlglucose/test/+mlglucose_unittest.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.
 	
	properties
 		registry
 		testObj
 	end

	methods (Test)
		function test_afun(this)
 			import mlglucose.*;
 			this.assumeEqual(1,1);
 			this.verifyEqual(1,1);
 			this.assertEqual(1,1);
 		end
	end

 	methods (TestClassSetup)
		function setupC11GlucoseKinetics(this)
 			import mlglucose.*;
 			this.testObj_ = C11GlucoseKinetics;
 		end
	end

 	methods (TestMethodSetup)
		function setupC11GlucoseKineticsTest(this)
 			this.testObj = this.testObj_;
 			this.addTeardown(@this.cleanFiles);
 		end
	end

	properties (Access = private)
 		testObj_
 	end

	methods (Access = private)
		function cleanFiles(this)
 		end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

