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
 		registry
        sessd
        sessionPath = '/data/nil-bluearc/raichle/PPGdata/jjlee2/HYGLY28'
 		testObj
        vnumber = 1
 	end

	methods (Test)
		function test_ctor(this)
 			glc = FDGKineticsWholebrain.godo2(this.sessd);
            glc = glc.aparcAsegBinarized_op_fdg;
            glc.plot;
        end
        function test_iterator(this)
            gkdIterator = this.testObj.iterator;
            iteration = 0;
            prods = {};
            while (gkdIterator.hasNext)
                aDirector = gkdIterator.next;
                aDirector = aDirector.constructRates;
                prods{iteration} = aDirector.product; %#ok<AGROW>
                iteration = iteration + 1;
            end
            this.verifyEqual(iteration, 1);
            this.verifyEqual(length(prods), 1);
        end
        function test_constructRates(this)
            this.testObj = this.testObj.constructRates;
        end
        function test_constructPhysiological(this)
            this.testObj = this.testObj.constructPhysiological;
        end
        function test_diagnose(this)
        end
        function test_plot(this)
        end
        function test_report(this)
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
            wellb  = mlcapintec.CapracBuilder('sessionData', this.sessd, ...
                'dtNyquist', this.DT_NYQUIST);  
            twilb  = mlswisstrace.TwiliteBuilder('sessionData', this.sessd, ...
                'dtNyquist', this.DT_NYQUIST);    
            solver = mlbayesian.BretthorstMcmc;
            
            % [15O]
            import mloxygen.*;
            solver.model = OyxgenModel( ...
                'scannerBuilder', scanb, 'aifBuilder', twilb, 'blindedData', blindd, ...
                'solverClass', 'mlbayesian.BretthorstMcmc', ...
                'sessionData', this.sessd);
            oxyd = OxygenKineticsDirector( ...
                'oxygenBldr', OxygenKineticsBuilder('solver', solver), ...
                'roisBldr', roisb);
            
            % [18F]DG
 			import mlglucose.*;
            solver.model = F18DeoxyGlucoseModel( ...
                'scannerBuilder', scanb, 'aifBuilder', wellb, 'blindedData', blindd, ...
                'solverClass', 'mlbayesian.BretthorstMcmc', ...
                'sessionData', this.sessd);  
            this.testObj = GlucoseKineticsDirector( ...
                'oxygenDir', oxyd, ...
                'glucoseBldr', F18DeoxyGlucoseKineticsBuilder('solver', solver), ...
                'roisBldr', roisb);
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

