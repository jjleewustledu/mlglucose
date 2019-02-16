classdef Test_GlucoseKineticsDirector < matlab.unittest.TestCase
	%% TEST_GLUCOSEKINETICSDIRECTOR 

	%  Usage:  >> results = run(mlraichle_unittest.Test_GlucoseKineticsDirector)
 	%          >> result  = run(mlraichle_unittest.Test_GlucoseKineticsDirector, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 27-Mar-2017 17:09:53 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlglucose/test/+mlglucose_unittest.
 	%% It was developed on Matlab 9.2.0.538062 (R2017a) for MACI64.  Copyright 2017 John Joowon Lee.
 	
	properties
        sessd
        sessionPath = '/data/nil-bluearc/raichle/PPGdata/jjlee2/HYGLY28'
 		testObj
        vnumber = 1
 	end

	methods (Test)
        function test_constructCmrglc_Stan(this)
            import mlglucose.*;
            fdg  = [];
            rois = [];
            bldr = FdgKineticsBuilderWithStan(fdg, rois);
            dtor = GlucoseKineticsDirector.createKinetics(bldr);
            dtor.constructCmrglc;
            this.verifyEqual(bldr.paramEstimate, nan);
        end
        function test_contructCmrglc_Bretthorst(this)
            import mlglucose.*;
            fdg  = [];
            rois = [];
            bldr = FdgKineticsBuilderWithBretthorst(fdg, rois);
            dtor = GlucoseKineticsDirector.createKinetics(bldr);
            dtor.constructCmrglc;
            this.verifyEqual(bldr.paramEstimate, nan);
        end
	end

 	methods (TestClassSetup)
		function setupGlucoseKineticsDirector(this)
            this.sessd = mlraichle.SessionData( ...
                'studyData', mlraichle.StudyData, ...
                'sessionPath', this.sessionPath, ...
                'ac', true);
            blindd = mlraichle.BlindedData('sessionData', this.sessd);
            roisb  = mlpet.BrainmaskBuilder('sessionData', this.sessd);
            scanb  = mlsiemens.BiographMMRBuilder('sessionData', this.sessd, ...
                'roisBuilder', roisb, ...
                'blindedData', blindd);
            wellb  = mlcapintec.CapracBuilder('sessionData', this.sessd);   
            solver = mlbayesian.BretthorstMcmc;
            
            % [18F]DG
 			import mlglucose.*;
            solver.model = F18DeoxyGlucoseModel( ...
                'scannerBuilder', scanb, 'aifBuilder', wellb, 'blindedData', blindd, ...
                'solverClass', 'mlbayesian.BretthorstMcmc', ...
                'sessionData', this.sessd);  
            this.testObj = GlucoseKineticsDirector.CreateKinetics;
 		end
	end

 	methods (TestMethodSetup)
		function setupF18DeoxyGlucoseKineticsTest(this)
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

