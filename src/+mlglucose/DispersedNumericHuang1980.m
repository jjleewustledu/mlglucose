classdef DispersedNumericHuang1980 < handle & mlglucose.Huang1980
	%% DISPERSEDNUMERICHUANG1980 

	%  $Revision$
 	%  was created 29-Apr-2020 23:31:01 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.7.0.1319299 (R2019b) Update 5 for MACI64.  Copyright 2020 John Joowon Lee.
 	
    properties (Constant)
        LENK = 5
    end
    
    methods (Static)
        function [this,fdg,aif] = createFromDeviceKit(devkit, varargin)
            %% adjusts AIF timings for coincidence of inflow with tissue activity from scanner
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param fdg is numeric, default from devkit.
            %  @param roi ...
            %  @param solver is in {'nest' 'simulanneal' 'hmc' 'lm' 'bfgs'}, default := 'simulanneal'.
            %  @param cbv is numeric or understood by mlfourd.ImagingContext2.
            %  @param blurFdg := {[], 0, 4.3, ...}                       
            %  @param map, default := mlglucose.Huang1980Model.preferredMap().
            %  @param times_sampled non-uniformly scheduled by the time-resolved PET reconstruction.
            %  @param artery_interpolated, default from devkit.
            %  @param glc is numeric, default from devkit.
            %  @param hct is numeric, default from devkit.
            %  @param LC is numeric, default from mloxygen.Raichle1983Model
            %  @param sigma0, default from mloptimization.SimulatedAnnealing.
            %  @param fileprefix, default from devkit.
            %  @return this.
            %  @return fdg, blurred by ipr.blurFdg.
            %  @return aif.
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'fdg', [], @isnumeric)
            addParameter(ip, 'roi', 'brain_222.4dfp.hdr')
            addParameter(ip, 'cbv', [])
            addParameter(ip, 'blurFdg', 4.3, @isnumeric)
            parse(ip, devkit, varargin{:})
            ipr = ip.Results;
            
            % prepare atlas data
            
            sesd = ipr.devkit.sessionData;
            sesd.jitOn222(sesd.wmparc1OnAtlas())
            sesd.jitOn222(sesd.fdgOnAtlas())
            
            % scanner provides calibrations, ancillary data
            
            ipr.roi = mlfourd.ImagingContext2(ipr.roi);
            roibin = ipr.roi.binarized();            
            scanner = ipr.devkit.buildScannerDevice();
            scanner = scanner.blurred(ipr.blurFdg);
            scanner = scanner.volumeAveraged(roibin);
            fdg = scanner.activityDensity(); % calibrated, decaying
            
            % v1
            
            fp = sprintf('mlglucose_DispersedHuang1980_createFromDeviceKit_dt%s', datestr(now, 'yyyymmddHHMMSS'));  
            if isnumeric(ipr.cbv)
                v1 = 0.0105*ipr.cbv;
            else
                v1 = mlfourd.ImagingContext2(ipr.cbv) .* 0.0105;
                v1 = v1.volumeAveraged(roibin);
            end
            
            % AIF            
            % Dt shifts the AIF in time:  Dt < 0 shifts left; Dt > 0 shifts right.
            
            counting = ipr.devkit.buildCountingDevice();
            Dt = mlglucose.DispersedNumericHuang1980.DTimeToShift(counting, scanner);
            aif = pchip(counting.times + Dt, counting.activityDensity(), 0:scanner.times(end));
            radm = counting.radMeasurements;
            
            this = mlglucose.DispersedNumericHuang1980( ...
                'devkit', devkit, ...
                'fdg', fdg, ...
                'solver', 'simulanneal', ...
                'v1', v1, ...
                'times_sampled', scanner.timesMid, ...
                'artery_interpolated', aif, ...
                'Dt', Dt, ...
                'glc', mlglucose.Huang1980.glcFromRadMeasurements(radm), ...
                'hct', mlglucose.Huang1980.hctFromRadMeasurements(radm), ...
                'fileprefix', fp, ...
                'roi', ipr.roi, ...
                varargin{:});
        end
        function Dt = DTimeToShift(varargin)
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
            unif_activity_s = pchip(t_s, activity_s, unif_t);
            d_activity_s = diff(unif_activity_s); % uniformly sampled time-derivative
            
            % shift dcv in time to match inflow with dtac    
            % use 0.1 of max since counting SNR >> 10 and idx_scanner ~ 1
            [~,idx_c] = max(activity_c > 0.1*max(activity_c));
            idx_c = idx_c - 1;
            [~,idx_s] = max(d_activity_s > 0.1*max(d_activity_s));
            Dt = unif_t(idx_s) - t_c(idx_c); % Dt ~ -20
            if Dt < -t_c(idx_c) || Dt > 0
                warning('mlglucose:ValueError', ...
                    'DispersedNumericHuang1980.DTimeToShift.Dt -> %g; forcing -> %g', Dt, -t_c(idx_c))
                Dt = -t_c(idx_c);
            end
        end
    end

	methods 
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
        
 		function this = DispersedNumericHuang1980(varargin)
 			%% DISPERSEDNUMERICHUANG1980   
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
            %  @param map is a containers.Map.  Default := this.preferredMap.
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

