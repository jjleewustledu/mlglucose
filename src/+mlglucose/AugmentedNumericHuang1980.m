classdef AugmentedNumericHuang1980 < handle & mlglucose.Huang1980
	%% AUGMENTEDNUMERICHUANG1980  

	%  $Revision$
 	%  was created 04-Jan-2021 20:38:13 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.9.0.1538559 (R2020b) Update 3 for MACI64.  Copyright 2021 John Joowon Lee.
 	 	
    properties (Constant)
        LENK = 5
    end
    
    properties
        Dt_aif
    end
    
    properties (Dependent)        
        artery_sampled
    end
    
    methods (Static)
        function [this,tac_,aif_] = createFromDeviceKit(devkit, devkit2, varargin)
            %% adjusts AIF timings for coincidence of inflow with tissue activity from scanner
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param required devkit2 is mlpet.IDeviceKit.
            %  @param required scanner is an mlpet.AbstractDevice.
            %  @param required scanner2 is an mlpet.AbstractDevice.
            %  @param required counter is an mlpet.AbstractDevice.
            %  @param required counter2 is an mlpet.AbstractDevice.
            %  @param solver is in {'nest' 'simulanneal' 'hmc' 'lm' 'bfgs'}, default := 'simulanneal'.
            %  @param roi is mlfourd.ImagingContext2.
            %  @param roi2 is mlfourd.ImagingContext2.
            %  @param cbv is mlfourd.ImagingContext2.
            %  @param cbv2 is mlfourd.ImagingContext2.
            %  @param blurFdg := {[], 0, 4.3, ...}   
            %  @param Dt_aif isscalar.
            %  @param times_sampled non-uniformly scheduled by the time-resolved PET reconstruction.
            %  @param artery_interpolated, default from devkit.
            %  @param glc is numeric, default from devkit.
            %  @param hct is numeric, default from devkit.
            %  @param LC is numeric, default from mloxygen.Raichle1983Model
            %  @param sigma0, default from mloptimization.SimulatedAnnealing.
            %  @param fileprefix, default from devkit.
            %  @param fracMixing in [0 1] for mixing tacs and aifs.
            %  @return this.
            %  @return tac_, blurred by ipr.blurFdg.
            %  @return aif_.
            
            import mlglucose.AugmentedNumericHuang1980
            import mlglucose.AugmentedNumericHuang1980.mix
            import mlglucose.Huang1980
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addRequired(ip, 'devkit2', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'scanner', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'scanner2', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'counter', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'counter2', [], @(x) isa(x, 'mlpet.AbstractDevice'))            
            addParameter(ip, 'roi', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'roi2', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'cbv', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'cbv2', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'blurFdg', 4.3, @isnumeric)
            addParameter(ip, 'Dt_aif', 0, @isscalar)
            addParameter(ip, 'fracMixing', 0.5, @isscalar)
            parse(ip, devkit, devkit2, varargin{:})
            ipr = ip.Results;
            this.Dt_aif = ipr.Dt_aif;
            
            % scanners provide calibrations, ancillary data            
            
            scanner = ipr.scanner.volumeAveraged(ipr.roi);
            tac = scanner.activityDensity();
            scanner2 = ipr.scanner2.volumeAveraged(ipr.roi2);            
            tac2 = scanner2.activityDensity();
            
            % counters also have calibrations
            
            aif = ipr.counter.activityDensity();         
            aif2 = ipr.counter2.activityDensity();
            
            % reconcile timings             
              
            if this.Dt_aif > 0 % counting2 is behind counting               
                Dt = AugmentedNumericHuang1980.DTimeToShift(ipr.counter, scanner);
                aif = makima(ipr.counter.times + Dt, aif, 0:scanner.times(end));
                tac = makima(scanner.times, tac, scanner.times); 
                aif2 = makima(ipr.counter2.times + Dt - this.Dt_aif, aif2, 0:scanner.times(end));
                tac2 = makima(scanner2.times - this.Dt_aif, tac2, scanner.times); 
                timesMid_ = scanner.timesMid;
            else              
                Dt = AugmentedNumericHuang1980.DTimeToShift(ipr.counter2, scanner2);
                aif2 = makima(ipr.counter2.times + Dt, aif2, 0:scanner2.times(end));
                tac2 = makima(scanner2.times, tac2, scanner2.times); 
                aif = makima(ipr.counter.times + Dt + this.Dt_aif, aif, 0:scanner2.times(end));
                tac = makima(scanner.times + this.Dt_aif, tac, scanner2.times);  
                timesMid_ = scanner2.timesMid;
            end 
            aif(aif < 0) = 0;
            tac(tac < 0) = 0;
            aif2(aif2 < 0) = 0;
            tac2(tac2 < 0) = 0;
            
            tac_ = mix(tac, tac2, ipr.fracMixing); % calibrated, decaying
            aif_ = mix(aif, aif2, ipr.fracMixing);
            
            % v1              

            v1 = ipr.cbv * 0.0105;
            v1 = v1.volumeAveraged(ipr.roi);
            v1 = v1.fourdfp.img;
            v12 = ipr.cbv2 * 0.0105;
            v12 = v12.volumeAveraged(ipr.roi2);
            v12 = v12.fourdfp.img;
            v1_ = mix(v1, v12, ipr.fracMixing);
            
            %
            
            glc = Huang1980.glcFromRadMeasurements(ipr.counter.radMeasurements);
            glc2 = Huang1980.glcFromRadMeasurements(ipr.counter2.radMeasurements);
            glc_ = mix(glc, glc2, ipr.fracMixing);
            hct = Huang1980.hctFromRadMeasurements(ipr.counter.radMeasurements);
            hct2 = Huang1980.hctFromRadMeasurements(ipr.counter2.radMeasurements);
            hct_ = mix(hct, hct2, ipr.fracMixing);
            fp = sprintf('mlglucose_AugmentedNumericHuang1980_createFromDeviceKit_dt%s', datestr(now, 'yyyymmddHHMMSS'));    
            
            this = AugmentedNumericHuang1980( ...
                'fdg', tac_, ...
                'solver', 'simulanneal', ...
                'devkit', ipr.devkit, ...
                'Dt', Dt, ...
                'times_sampled', timesMid_, ...
                'artery_interpolated', aif_, ...
                'v1', v1_, ...
                'glc', glc_, ...
                'hct', hct_, ...
                'fileprefix', fp, ...
                varargin{:});
        end
        function Dt = DTimeToShift(varargin)
            %% @return Dt is the delay of the tac behind the aif, but typically Dt ~ -20.
            
            ip = inputParser;
            addRequired(ip, 'counter')
            addRequired(ip, 'scanner')
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            t_c        = asrow(ipr.counter.times);
            activity_c = asrow(ipr.counter.activityDensity());
            t_s        = asrow(ipr.scanner.timesMid);
            activity_s = asrow(ipr.scanner.imagingContext.fourdfp.img);            
            
            unif_t = 0:max([t_c t_s]);
            unif_activity_c = makima(t_c, activity_c, unif_t);
            unif_activity_c(unif_activity_c < 0) = 0;
            unif_activity_s = makima(t_s, activity_s, unif_t);
            unif_activity_s(unif_activity_s < 0) = 0;
            d_activity_s = diff(unif_activity_s); % uniformly sampled time-derivative
            
            % shift dcv in time to match inflow with dtac    
            % use 0.1 of max since counting SNR >> 10 and idx_scanner ~ 1
            [~,idx_c] = max(unif_activity_c > 0.1*max(unif_activity_c));
            [~,idx_s] = max(d_activity_s > 0.1*max(d_activity_s));
            Dt = unif_t(idx_s) - unif_t(idx_c); % Dt ~ -20
            if Dt < -unif_t(idx_c) || Dt > 0
                warning('mlglucose:ValueError', ...
                    'AugmentedNumericHuang1980.DTimeToShift.Dt -> %g; forcing -> %g', Dt, -t_c(idx_c))
                Dt = -unif_t(idx_c);
            end
        end
        function Dt = DTimeToShiftAifs(varargin)
            %% @return Dt is the delay of counting2 behind counting.
            
            ip = inputParser;
            addRequired(ip, 'counting')
            addRequired(ip, 'counting2')
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            t         = asrow(ipr.counting.times);
            activity  = asrow(ipr.counting.activityDensity());
            t2        = asrow(ipr.counting2.times);
            activity2 = asrow(ipr.counting2.activityDensity());
            
            unif_t = 0:max([t t2]);
            unif_activity = makima(t, activity, unif_t);
            unif_activity2 = makima(t2, activity2, unif_t);
              
            % shift activity in time to match inflow with activity2
            % use 0.1 of max since counting SNR >> 10 
            [~,idx] = max(unif_activity > 0.1*max(unif_activity));
            [~,idx2] = max(unif_activity2 > 0.1*max(unif_activity2));
            Dt = unif_t(idx2) - unif_t(idx);
        end
        function mixed = mix(obj, obj2, f)
            assert(f > 0)
            assert(f < 1)
            if isnumeric(obj) && isnumeric(obj2)
                mixed = f*obj + (1 - f)*obj2;
                return
            end
            obj = mlfourd.ImagingContext2(obj);
            obj2 = mlfourd.ImagingContext2(obj2);
            mixed = obj * f + obj2 * (1 - f);
        end
    end
    
    methods
        
        %% GET
        
        function g = get.artery_sampled(this)
            g = this.model.sampleOnScannerFrames(this.artery_interpolated, this.times_sampled);
        end
        
        %%
        
        function ks = buildKs(this, varargin)
            this = solve(this, varargin{:});
            ks = [k1(this) k2(this) k3(this) k4(this) k5(this)];
        end
        function [k,sk] = k5(this, varargin)
            [k,sk] = k5(this.strategy_, varargin{:});
        end
        function [K,sK] = Ks(this, varargin)
            K = zeros(1,this.LENK);
            sK = zeros(1,this.LENK);
            [K(1),sK(1)] = K1(this.strategy_, varargin{:});
            [K(2),sK(2)] = k2(this.strategy_, varargin{:});
            [K(3),sK(3)] = k3(this.strategy_, varargin{:});
            [K(4),sK(4)] = k4(this.strategy_, varargin{:});
            [K(5),sK(5)] = k5(this.strategy_, varargin{:});
        end
        function [k,sk] = ks(this, varargin)
            k = zeros(1,this.LENK);
            sk = zeros(1,this.LENK);
            [k(1),sk(1)] = k1(this.strategy_, varargin{:});
            [k(2),sk(2)] = k2(this.strategy_, varargin{:});
            [k(3),sk(3)] = k3(this.strategy_, varargin{:});
            [k(4),sk(4)] = k4(this.strategy_, varargin{:});
            [k(5),sk(5)] = k5(this.strategy_, varargin{:});
        end
        function fdg = checkSimulated(this, varargin)
            %% CHECKSIMULATED simulates tissue activity with passed and internal parameters without changing state.
            %  @param required ks is [k1 k2 k3 k4 k5 Dt].
            %  @param v1 is CBV < 1 and dimensionless; default is this.v1.
            %  @param aif is numeric; default is this.artery_interpolated for model state.
            %  @return fdg simulation is numeric.
        
            ip = inputParser;
            addOptional(ip, 'ks', this.ks(), @isnumeric)
            addParameter(ip, 'v1', this.v1, @isscalar)
            addParameter(ip, 'aif', this.artery_interpolated, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;  
            
            fdg = this.model.simulated(ipr.ks, 'v1', ipr.v1, 'aif', ipr.aif, 'Dt', this.Dt);
        end
    end
    
    %% PROTECTED

	methods (Access = protected)
 		function this = AugmentedNumericHuang1980(varargin)
 			%% AUGMENTEDNUMERICHUANG1980
            %  @param fdg is numeric.
            %  @param solver is in {'nest' 'simulanneal' 'hmc' 'lm' 'bfgs'}. 
            %  @param devkit is mlpet.IDeviceKit.          
            %  @param Dt is numeric, s of time-shifting for AIF.
            %  @param times_sampled is numeric.
            %  @param artery_interpolated is numeric.     
            %  
            %  for mlglucose.DispersedHuang1980Model: 
            %  @param v1
            %  @param glc
            %  @param hct
            %  @param LC
            %  @param map is a containers.Map.  Default := DispersedHuang1980Model.preferredMap.
            %  @param times_sampled for scanner is typically not uniform.
            %  @param artery_interpolated must be uniformly interpolated.
            %
            %  for mlglucose.DispersedHuang1980SimulAnneal:
            %  @param context is mlglucose.Huang1980.
            %  @param sigma0.
            %  @param fileprefix.

 			this = this@mlglucose.Huang1980( ...
                'model', mlglucose.DispersedHuang1980Model(varargin{:}), ...
                varargin{:});	
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'fdg', [], @(x) isnumeric(x))
            addParameter(ip, 'solver', 'simulanneal', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
                        
            this.measurement = ipr.fdg;
            switch lower(ipr.solver)
                case 'simulanneal'
                    this.strategy_ = mlglucose.DispersedHuang1980SimulAnneal( ...
                        'context', this, varargin{:});
                otherwise
                    error('mlglucose:NotImplementedError', 'Huang1980.ipr.solver->%s', ipr.solver)
            end
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

