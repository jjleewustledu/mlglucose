classdef NumericHuang1980 < handle & mlglucose.Huang1980
	%% NUMERICHUANG1980 

	%  $Revision$
 	%  was created 29-Apr-2020 23:31:01 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.7.0.1319299 (R2019b) Update 5 for MACI64.  Copyright 2020 John Joowon Lee.
 	
    properties
    end
    
    methods (Static)
        function this = createFromDeviceKit(devkit, varargin)
            %% adjusts AIF timings for coincidence of inflow with tissue activity from scanner
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param
            %  @param fdg is numeric, default from devkit.
            %  @param solver is in {'nest' 'simulanneal' 'hmc' 'lm' 'bfgs'}, default := 'simulanneal'.
            %  @param map, default := mlglucose.Huang1980Model.preferredMap().
            %  @param cbv is numeric or understood by mlfourd.ImagingContext2.
            %  @param times_sampled non-uniformly scheduled by the time-resolved PET reconstruction.
            %  @param artery_interpolated, default from devkit.
            %  @param glc is numeric, default from devkit.
            %  @param hct is numeric, default from devkit.
            %  @param LC is numeric, default from mlglucose.Huang1980Model
            %  @param sigma0, default from mloptimization.SimulatedAnnealing.
            %  @param fileprefix, default from devkit.
            %  @param blurFdg := {[], 0, 4.3, ...}
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'fdg', [], @isnumeric)
            addParameter(ip, 'roi', 'brain.4dfp.hdr', @(x) isa(x, 'mlfourd.ImagingContext2'))
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
            
            fp = sprintf('mlglucose_Huang1980_createFromDeviceKit_dt%s', datestr(now, 'yyyymmddHHMMSS'));  
            if isnumeric(ipr.cbv)
                v1 = 0.0105*ipr.cbv;
            else
                v1 = mlfourd.ImagingContext2(ipr.cbv) .* 0.0105;
                v1 = v1.volumeAveraged(roibin);
            end
            
            % AIF            
            % Dt shifts the AIF in time:  Dt < 0 shifts left; Dt > 0 shifts right.
            
            counting = ipr.devkit.buildCountingDevice();
            Dt = mlglucose.NumericHuang1980.DTimeToShift(counting, scanner);
            aif = pchip(counting.times + Dt, counting.activityDensity(), 0:scanner.times(end));
            radm = counting.radMeasurements;
            
            this = mlglucose.NumericHuang1980( ...
                devkit, ...
                'fdg', fdg, ...
                'solver', 'simulanneal', ...
                'v1', v1, ...
                'times_sampled', scanner.timesMid, ...
                'artery_interpolated', aif, ...
                'Dt', Dt, ...
                'glc', mlglucose.Huang1980.glcFromRadMeasurements(radm), ...
                'hct', mlglucose.Huang1980.hctFromRadMeasurements(radm), ...
                'fileprefix', fp, ...
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
            [~,idx_s] = max(d_activity_s > 0.1*max(d_activity_s));
            Dt = unif_t(idx_s) - t_c(idx_c); % Dt ~ -20
            if Dt < -t_c(idx_c) || Dt > 0
                warning('mlglucose:ValueError', ...
                    'NumericHuang1980.DTimeToShift.Dt -> %g; forcing -> %g', Dt, -t_c(idx_c))
                Dt = -t_c(idx_c);
            end
        end
    end

	methods 		  
 		function this = NumericHuang1980(devkit, varargin)
 			%% NUMERICHUANG1980
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param fdg is numeric.
            %  @param solver is in {'nest' 'simulanneal' 'hmc' 'lm' 'bfgs'}.
            %  @param map
            %  @param v1
            %  @param times_sampled non-uniformly scheduled by the time-resolved PET reconstruction.
            %  @param artery_interpolated
            %  @param Dt is numeric, s of time-shifting for AIF.
            %  @param glc is numeric, mmol/L.
            %  @param hct is numeric, percent.
            %  @param LC is numeric.
            %  @param sigma0.
            %  @param fileprefix.

 			this = this@mlglucose.Huang1980(devkit, varargin{:});	
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'fdg', [], @(x) isnumeric(x))
            addParameter(ip, 'solver', 'simulanneal', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
                        
            this.measurement = ipr.fdg;
            this.model = mlglucose.Huang1980Model(varargin{:});
            switch lower(ipr.solver)
                case 'nest'
                    this.strategy_ = mlglucose.Huang1980Nest( ...
                        'context', this, varargin{:});
                case 'simulanneal'
                    this.strategy_ = mlglucose.Huang1980SimulAnneal( ...
                        'context', this, varargin{:});
                case 'hmc'
                    this.strategy_ = mlglucose.Huang1980HMC( ...
                        'context', this, varargin{:});
                case 'lm'
                    this.strategy_ = mlglucose.Huang1980LM( ...
                        'context', this, varargin{:});
                case 'bfgs'
                    this.strategy_ = mlglucose.Huang1980BFGS( ...
                        'context', this, varargin{:});
                otherwise
                    error('mlglucose:NotImplementedError', 'Huang1980.ipr.solver->%s', ipr.solver)
            end
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

