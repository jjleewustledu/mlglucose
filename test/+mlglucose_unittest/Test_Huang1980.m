classdef Test_Huang1980 < matlab.unittest.TestCase
	%% TEST_HUANG1980 

	%  Usage:  >> results = run(mlglucose_unittest.Test_Huang1980)
 	%          >> result  = run(mlglucose_unittest.Test_Huang1980, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 10-Apr-2020 16:31:42 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/test/+mlglucose_unittest.
 	%% It was developed on Matlab 9.7.0.1319299 (R2019b) Update 5 for MACI64.  Copyright 2020 John Joowon Lee.
 	
	properties
        artery = [15.7231339417783 40.0670951423858 32.7896589986613 60.9033324985914 115.390827456407 20247.5870083001 138101.363561975 247624.893301051 219962.735811158 123800.355195601 84547.8333286971 74641.7780423665 61033.7705200579 51488.086840046 48137.398517549 45853.6272536275 43105.9076185098 40588.4399463862 39532.0074392734 37376.951543845 38029.7572601973 36896.5536596196 35837.3996944633 35239.9976922434 33563.3578174873 31172.6742405234 31661.2958727327 31726.6836574803 31190.0410893731 29736.9838758955 29345.5491772129 28612.6204028212 28087.9993320967 12488.6629201451 8478.08886746032 5324.32339588189 3353.30521211569 945.745212317361 665.76566052086]
        arteryTimesSampled = [0 3 7 10 14 18 22 25 30 35 39 44 48 52 55 58 63 67 70 74 78 81 84 90 93 95 98 103 106 109 113 117 120 417 597 897 1197 2997 3597]
        brainBinarized
        cbvfp = 'cbvdt20190523122016_on_T1001_decayUncorrect0'
        fdgfp = 'fdgdt20190523132832'
        fdgTausSampled = [10 13 14 16 17 19 20 22 23 25 26 28 29 31 32 34 35 37 38 40 41 43 44 46 47 49 50 52 53 56 57 59 60 62 63 65 66 68 69 71 72 74 76 78 79 81 82 84 85 87 88 91 92 94 95 97 98 100 101 104 105 108]
        home
        Laccumb % 683 voxels
        map
        map1
        pthRestricted = fullfile(getenv('SINGULARITY_HOME'), 'subjects', 'sub-s58163', 'resampling_restricted', '')
        pwd0
 		registry
        sesd
        sesf = 'CCIR_00559/ses-E03056/FDG_DT20190523132832.000000-Converted-AC'
        subf = 'subjects/sub-S58163/resampling_restricted'
        testImgObj
 		testObj
        v1 = 0.038
        v1ic
    end
    
    properties (Dependent)
        arteryInterpolated
        fdg
        fdgTimesSampled
        ks
    end
    
    methods 
        
        %% GET
        
        function g = get.arteryInterpolated(this)
            g = pchip(this.arteryTimesSampled, this.artery, 0:this.fdgTimesSampled(end));
        end
        function g = get.fdg(this)
            if ~isfile(fullfile(this.pthRestricted, [this.fdgfp '_on_T1001.4dfp.hdr']))
                pwd0 = pushd(this.pthRestricted);
                fv = mlfourdfp.FourdfpVisitor();
                t4 = [this.fdgfp '_to_T1001_t4'];
                fv.t4img_4dfp(t4, this.fdgfp, 'options', '-OT1001')            
                popd(pwd0)
            end            
            ic = mlfourd.ImagingContext2(fullfile(this.pthRestricted, [this.fdgfp '_on_T1001.4dfp.hdr']));
            g = ic.fourdfp.img;
        end
        function g = get.fdgTimesSampled(this)
            g = cumsum(this.fdgTausSampled);
            g = g - this.fdgTausSampled/2;
        end
        function g = get.ks(this)
            g = nan(1, this.map.Count);
            for i = 1:this.map.Count
                lbl = sprintf('k%i', i);
                g(i) = this.map(lbl).init; 
            end            
        end
    end

	methods (Test)
		function test_afun(this)
 			import mlglucose.*;
 			this.assumeEqual(1,1);
 			this.verifyEqual(1,1);
 			this.assertEqual(1,1);
        end
        function test_ctor(this)
            disp(this.testObj)
            disp(this.testObj.model)
            for k = this.testObj.model.map.keys
                fprintf('%s:\n', k{1})
                disp(this.testObj.model.map(k{1}))
            end
        end
        function test_imaging_createFromDeviceKit(this)
            devkit = mlsiemens.BiographMMRKit.createFromSession(this.sesd);
            o = mlglucose.ImagingHuang1980.createFromDeviceKit(devkit, 'cbv', [this.cbvfp '.4dfp.hdr'], 'roi', this.Laccumb);
            disp(o)
            disp(o.model)
            this.verifyEqual(o.glc, 104, 'RelTol', 1e-4)
            this.verifyEqual(o.hct, 0.398, 'RelTol', 1e-4)
            this.verifyEqual(o.LC, 0.81, 'RelTol', 1e-4)
            this.verifyEqual(o.times_sampled, this.fdgTimesSampled, 'RelTol', 1e-4)
            o.fdg.fsleyes()
            o.roi.fsleyes()
            o.v1.fsleyes()
        end
        function test_imaging_solve(this)
            devkit = mlsiemens.BiographMMRKit.createFromSession(this.sesd);
            o = mlglucose.ImagingHuang1980.createFromDeviceKit(devkit, 'cbv', [this.cbvfp '.4dfp.hdr'], 'roi', this.Laccumb);
            tic
            o = o.solve();
            toc
            save(o.ks)
            o.fdg.fsleyes()
            o.roi.fsleyes()
            o.v1.fsleyes()
            o.ks.fsleyes()
        end
        function test_numeric_createFromDeviceKit(this)
            devkit = mlsiemens.BiographMMRKit.createFromSession(this.sesd);
            o = mlglucose.Huang1980.createFromDeviceKit(devkit, 'cbv', 100*this.v1);
            disp(o)
            this.verifyEqual(o.measurement, ...
                single([0 3323.689453125 3802.98388671875 3954.17602539062 4092.505859375 4287.32470703125 4538.7744140625 4752.529296875 4938.9716796875 5151.2509765625 5268.09765625 5378.0703125 5481.6123046875 5642.17626953125 5710.3193359375 5803.84765625 5858.65771484375 5943.60791015625 5991.474609375 6044.99072265625 6091.6650390625 6107.958984375 6164.14306640625 6182.19091796875 6233.5361328125 6238.93359375 6265.95166015625 6284.30859375 6286.5458984375 6327.162109375 6336.34521484375 6337.15771484375 6337.69482421875 6334.58984375 6346.42431640625 6334.595703125 6336.93115234375 6335.9033203125 6305.41796875 6319.42431640625 6282.357421875 6269.5634765625 6241.36767578125 6246.48486328125 6217.9326171875 6226.42236328125 6257.197265625 6296.0625 6291.29443359375 6279.97705078125 6287.71533203125 6268.08251953125 6272.228515625 6262.6416015625 6258.60302734375 6236.263671875 6238.8759765625 6230.33203125 6209.56689453125 6193.8408203125 6185.36083984375 6120.12255859375]), ...
                'RelTol', 1e-4)
            disp(o.model)
            m = o.model;
            this.verifyEqual(m.glc, 104, 'RelTol', 1e-4)
            this.verifyEqual(m.hct, 0.398, 'RelTol', 1e-4)
            this.verifyEqual(m.LC, 0.81, 'RelTol', 1e-4)
            this.verifyEqual(m.times_sampled, this.fdgTimesSampled, 'RelTol', 1e-4)
            this.verifyEqual(double(m.v1), this.v1, 'RelTol', 1e-4)
        end
        function test_fdg(this)
            ic = mlfourd.ImagingContext2([this.fdgfp '_on_T1001.4dfp.hdr']);
            ic.fsleyes('T1001.4dfp.hdr')
        end
        function test_nest(this) % with synthetic data
            %diary(sprintf('Test_Huang1980_test_nest_%s.log', datestr(now, 'yyyymmddHHMMSS')))
 			o = mlglucose.Huang1980( ...
                'map', this.map, ...
                'times_sampled', this.fdgTimesSampled, ...
                'artery_interpolated', this.arteryInterpolated, ...
                'v1', this.v1, ...
                'solver', 'nest');
            o = o.simulated(this.ks);
            o = o.solve();
            %this.verifySimulated(o)
            %o = o.solve('STEP_Initial', [5e-4 5e-3 5-e2]);
            %o = o.solve('n', [10 50 100 500]);
            %o = o.solve('MCMC_Counter', [10 25 50 100 250 500]);
            %o = o.solve('MAX', [500 1000 2000 4000 8000]);
            %diary('off')
            
            %% BEST SO FAR:
            % -----------------------------------------------------------
            % # iterates ~ nest = 1000 <= MAX = 1000
            % # sampling particles ~ n = 10.000000
            % MCMC_Counter = 50.000000
            % STEP_Initial = 0.100000
            % Stopping criteria = 31.351198
            % Evidence:  ln(Z) = -3.959046 +/- 0.564772
            % Information:  H = 3.189671 nats = 4.601722 bits
            % Model:
            % 	k1 = 0.161748 +/- 0.024555
            % 	k2 = 0.362055 +/- 0.062365
            % 	k3 = 0.132983 +/- 0.044052
            % 	k4 = 0.000517 +/- 0.000101
            % 	sigma0 = 0.500000
            % 	map('k1') => min = 0.1, max = 0.4, init = 0.2
            % 	map('k2') => min = 0.15, max = 0.6, init = 0.3
            % 	map('k3') => min = 0.05, max = 0.2, init = 0.1
            % 	map('k4') => min = 0.00025, max = 0.001, init = 0.0005
            % 	sampled logL (k = 1000) = -0.000274
            % 	sampled logWt(k = 1000) = -102.252442
            % 
            % Elapsed time is 148.670101 seconds.
        end
        function test_simulAnneal(this) % with synthetic data
 			o = mlglucose.Huang1980( ...
                'map', this.map1, ...
                'times_sampled', this.fdgTimesSampled, ...
                'artery_interpolated', this.arteryInterpolated, ...
                'v1', this.v1, ...
                'solver', 'simulanneal');
            o = o.simulated(this.ks);
            o = o.solve();
            this.verifySimulated(o)

            %% BEST SO FAR:
            % runtests('mlglucose_unittest.Test_Huang1980','ProcedureName','test_simulAnneal')
            % Running mlglucose_unittest.Test_Huang1980
            % 
            % Local minimum possible. Constraints satisfied.
            % 
            % fmincon stopped because the size of the current step is less than
            % the value of the step size tolerance and constraints are 
            % satisfied to within the value of the constraint tolerance.
            % 
            % <stopping criteria details>
            % Maximum number of function evaluations exceeded: increase options.MaxFunctionEvaluations.
            % FMINCON: 
            % 
            % Local minimum possible. Constraints satisfied.
            % 
            % fmincon stopped because the size of the current step is less than
            % the value of the step size tolerance and constraints are 
            % satisfied to within the value of the constraint tolerance.
            % 
            % <stopping criteria details>
            % 
            % Optimization stopped because the relative changes in all elements of x are
            % less than options.StepTolerance = 1.000000e-10, and the relative maximum constraint
            % violation, 0.000000e+00, is less than options.ConstraintTolerance = 1.000000e-06.
            % 
            % 
            % 
            % mlglucose.Huang1980SimulAnneal         v1:    0.038000000000000
            % 
            % initial ks0:    0.127114565124791   0.126481312097762   0.053766144089710   0.000265205847209
            % 
            % est.     ks:    0.200214191190547   0.300638900430815   0.100068134839889   0.000499796756914
            % 
            %         loss:      1.609666696526757e-05
            % 
            %    exitflag:      0
            % 
            %      iterations: 11956
            %       funccount: 12441
            %         message: 'Maximum number of function evaluations exceeded: increase options.MaxFunctionEvaluations.?FMINCON: ??Local minimum possible. Constraints satisfied.??fmincon stopped because the size of the current step is less than?the value of the step size tolerance and constraints are ?satisfied to within the value of the constraint tolerance.??<stopping criteria details>??Optimization stopped because the relative changes in all elements of x are?less than options.StepTolerance = 1.000000e-10, and the relative maximum constraint?violation, 0.000000e+00, is less than options.ConstraintTolerance = 1.000000e-06.??'
            %        rngstate: [1×1 struct]
            %     problemtype: 'boundconstraints'
            %     temperature: [4×1 double]
            %       totaltime: 33.534028894999999
            % 
            %      Type: 'twister'
            %      Seed: 0
            %     State: [625×1 uint32]
            % 
            %    0.002029075495160
            %    0.001894381488540
            %    0.001877016902231
            %    0.002162491288036
            % 
            % Model:
            % 	k1 = 0.200214
            % 	k2 = 0.300639
            % 	k3 = 0.100068
            % 	k4 = 0.000500
            % 	sigma0 = 0.050000
            % 	map('k1') => min = 0.02, max = 2, init = 0.127114565124791
            % 	map('k2') => min = 0.03, max = 3, init = 0.126481312097762
            % 	map('k3') => min = 0.01, max = 1, init = 0.0537661440897098
            % 	map('k4') => min = 5e-05, max = 0.005, init = 0.00026520584720874
            % .
            % Done mlglucose_unittest.Test_Huang1980
            % __________
            % 
            % 
            % ans = 
            % 
            %   TestResult with properties:
            % 
            %           Name: 'mlglucose_unittest.Test_Huang1980/test_simulAnneal'
            %         Passed: 1
            %         Failed: 0
            %     Incomplete: 0
            %       Duration: 33.683903550000004
            %        Details: [1×1 struct]
            % 
            % Totals:
            %    1 Passed, 0 Failed, 0 Incomplete.
            %    33.6839 seconds testing time.
        end
        function test_hmc(this) % with synthetic data
            o = o.simulated(this.ks);
            o = o.solve();
            this.verifySimulated(o)
        end
        function test_bfgs(this) % with synthetic data
            o = o.simulated(this.ks);
            o = o.solve();
            this.verifySimulated(o)
        end
        function test_lm(this) % with synthetic data
            o = o.simulated(this.ks);
            o = o.solve();
            this.verifySimulated(o)
        end
        function test_cbv(this)
            cbv = this.testObj.buildCbv();
            assert(1 == numel(cbv))
            this.verifyEqual(cbv.fourdfp.img, 100*this.v1, 'RelTol', 1e-4)
        end
        function test_cmrglc(this)
            cmrglc = this.testObj.cmrglc();
            cmrglc.fsleyes()
            cmrglc = cmrglc.volumeAveraged(this.brainBinarized);
            this.verifyEqual(cmrglc.img, [], 'RelTol', 1e-4)
        end
        function test_K1(this)
            K1 = this.testObj.K1();
            K1.fsleyes()
            K1 = K1.volumeAveraged(this.brainBinarized);
            this.verifyEqual(K1.img, [], 'RelTol', 1e-4)
        end
        function test_numeric_ks(this)
            devkit = mlsiemens.BiographMMRKit.createFromSession(this.sesd);
            o = mlglucose.NumericHuang1980.createFromDeviceKit(devkit, 'cbv', 100*this.v1, 'roi', this.Laccumb);
            
            tic
            o = o.solve();
            toc
            this.verifyEqual([o.k1() o.k2() o.k3() o.k4()], [0.013943826385222   0.000449068811131   0.000586885759740   0.000000000000035], 'RelTol', 1e-2)
            
            % Model:
            % 	k1 = 0.013944
            % 	k2 = 0.000449
            % 	k3 = 0.000587
            % 	k4 = 0.000000
            % 	sigma0 = 0.050000
            % 	map('k1') => min = 0.001, max = 0.5, init = 0.03, sigma = 0.02
            % 	map('k2') => min = 2.22044604925031e-16, max = 0.1, init = 0.003, sigma = 0.01
            % 	map('k3') => min = 2.22044604925031e-16, max = 0.05, init = 0.003, sigma = 0.002
            % 	map('k4') => min = 2.22044604925031e-16, max = 0.01, init = 0.0003, sigma = 0.001
            % Elapsed time is 21.416997 seconds.
        end
        function test_v1ic(this)
            this.verifyEqual(size(this.v1ic), [256 256 256])
            this.verifyEqual(dipsum(double(this.v1ic)), 7.6588078e+04, 'RelTol', 1e-4)
            this.verifyEqual(dipmean(double(this.v1ic)), 0.0045650, 'RelTol', 1e-4)
            %this.v1ic.fsleyes()
        end
	end

 	methods (TestClassSetup)
		function setupHuang1980(this)
 			import mlglucose.*
            this.home = fullfile(getenv('SINGULARITY_HOME'), this.subf);
            this.pwd0 = pushd(this.home);
            this.sesd = mlraichle.SessionData.create(this.sesf);
            this.map = containers.Map;
            this.map('k1') = struct('min',    0.02, 'max', 2,    'init', 0.2);
            this.map('k2') = struct('min',    0.03, 'max', 3,    'init', 0.3);
            this.map('k3') = struct('min',    0.01, 'max', 1,    'init', 0.1);
            this.map('k4') = struct('min',    5e-5, 'max', 5e-3, 'init', 5e-4);
            this.map1 = containers.Map;
            this.map1('k1') = struct('min',    0.02, 'max', 2,    'init', 0.2*rand());
            this.map1('k2') = struct('min',    0.03, 'max', 3,    'init', 0.3*rand());
            this.map1('k3') = struct('min',    0.01, 'max', 1,    'init', 0.1*rand());
            this.map1('k4') = struct('min',    5e-5, 'max', 5e-3, 'init', 5e-4*rand());
            this.v1ic = mlfourd.ImagingContext2([this.cbvfp '.4dfp.hdr'])/100;
            this.Laccumb = mlfourd.ImagingFormatContext('wmparc_222.4dfp.hdr');
            this.Laccumb.img = this.Laccumb.img == 26;
            this.Laccumb = mlfourd.ImagingContext2(this.Laccumb, 'fileprefix', 'wmparc_222_26');
 			
%            devkit = mlsiemens.BiographMMRKit.createFromSession(this.sesd);
%            this.testObj_ = mlglucose.Huang1980.createFromDeviceKit(devkit, 'cbv', 100*this.v1);
 		end
	end

 	methods (TestMethodSetup)
		function setupHuang1980Test(this)            
% 			this.testObj = copy(this.testObj_);
            pushd(this.home);
 			this.addTeardown(@this.cleanTestMethod);
 		end
    end

    %% PRIVATE
    
	properties (Access = private)
        fdg_
        testImgObj_
 		testObj_
 	end

	methods (Access = private)
		function cleanTestMethod(this)
            popd(this.pwd0)
        end
        function verifySimulated(this, o)
            this.verifyEqual(o.K1(), this.v1*this.ks(1), 'RelTol', 0.1)
            this.verifyEqual(o.k1(), this.ks(1), 'RelTol', 0.1)
            this.verifyEqual(o.k2(), this.ks(2), 'RelTol', 0.1)
            this.verifyEqual(o.k3(), this.ks(3), 'RelTol', 0.1)
            this.verifyEqual(o.k4(), this.ks(4), 'RelTol', 0.1)
            %this.verifyEqual(o.t0(), [], 'RelTol', 0.1)
            %this.verifyEqual(o.cmrglc(), [], 'RelTol', 0.1)
        end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

