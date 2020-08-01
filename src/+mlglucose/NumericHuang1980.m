classdef NumericHuang1980 < handle & mlglucose.Huang1980
	%% NUMERICHUANG1980  

	%  $Revision$
 	%  was created 29-Apr-2020 23:31:01 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.7.0.1319299 (R2019b) Update 5 for MACI64.  Copyright 2020 John Joowon Lee.
 	
    methods (Static)
        function this = createFromDeviceKit(devkit, varargin)
            %% 
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
            
            import mlfourd.ImagingContext2
            import mlglucose.Huang1980.glcFromRadMeasurements
            import mlglucose.Huang1980.hctFromRadMeasurements
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'fdg', [], @isnumeric)
            addParameter(ip, 'roi', 'brain.4dfp.hdr', @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'cbv', [])
            addParameter(ip, 'blurFdg', 4.3, @isnumeric)
            parse(ip, devkit, varargin{:})
            ipr = ip.Results;
            
            % scanner provides calibrations, ancillary data
            
            ipr.roi = mlfourd.ImagingContext2(ipr.roi);
            roibin = ipr.roi.binarized();            
            scanner = ipr.devkit.buildScannerDevice();
            if ~isempty(ipr.blurFdg) && ipr.blurFdg > 0
                scanner = scanner.blurred(ipr.blurFdg);
            end
            scanner = scanner.volumeAveraged(roibin);
            fdg = scanner.activityDensity(); % calibrated
            
            % v1
            fp = sprintf('mlglucose_Huang1980_createFromDeviceKit_dt%s', datestr(now, 'yyyymmddHHMMSS'));  
            if isnumeric(ipr.cbv)
                v1 = 0.0105*ipr.cbv;
            else
                v1 = 0.0105*mlfourd.ImagingContext2(ipr.cbv);
                v1 = v1.volumeAveraged(roibin);
            end
            
            % AIF            
            % Dt shifts the AIF in time:  Dt < 0 shifts left; Dt > 0 shifts right.
            counting = ipr.devkit.buildCountingDevice();
            Dt = mlglucose.NumericHuang1980.shiftAif(counting, scanner);
            aif = pchip(counting.times + Dt, counting.activityDensity(), 0:scanner.times(end));
            radm = counting.radMeasurements;
            
            this = mlglucose.NumericHuang1980( ...
                'fdg', fdg, ...
                'solver', 'simulanneal', ...
                'v1', v1, ...
                'times_sampled', scanner.timesMid, ...
                'artery_interpolated', aif, ...
                'glc', glcFromRadMeasurements(radm), ...
                'hct', hctFromRadMeasurements(radm), ...
                'fileprefix', fp, ...
                varargin{:});
        end
        function Dt = shiftAif(varargin)
            ip = inputParser;
            addRequired(ip, 'counter')
            addRequired(ip, 'scanner')
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            t_c        = asrow(ipr.counter.times);
            activity_c = asrow(ipr.counter.activityDensity());
            t_s        = asrow(ipr.scanner.times);
            activity_s = asrow(ipr.scanner.imagingContext.fourdfp.img);
            
            unif_t = 0:max([t_c t_s]);
            unif_activity_s = pchip(t_s, activity_s, unif_t);
            d_activity_s = diff(unif_activity_s); % uniformly sampled time-derivative
            
            % shift dcv in time to match inflow with dtac            
            [~,idx_c] = max(activity_c > max(activity_c)/2);
            [~,idx_s] = max(d_activity_s > max(d_activity_s)/2);
            Dt = unif_t(idx_s) - t_c(idx_c); % Dt ~ -10
            if Dt > 0
                warning('mlglucose:ValueError', 'NumericHuang1980.shiftAif.Dt -> %g', Dt)
                Dt = -20;
            end
            if Dt < -300
                warning('mlglucose:ValueError', 'NumericHuang1980.shiftAif.Dt -> %g', Dt)
                Dt = -20;
            end
        end 
    end

	methods 		  
 		function this = NumericHuang1980(varargin)
 			%% NUMERICHUANG1980
            %  @param fdg is numeric.
            %  @param solver is in {'nest' 'simulanneal' 'hmc' 'lm' 'bfgs'}.
            %  @param map
            %  @param v1
            %  @param times_sampled non-uniformly scheduled by the time-resolved PET reconstruction.
            %  @param artery_interpolated
            %  @param glc is numeric, mmol/L.
            %  @param hct is numeric, percent.
            %  @param LC is numeric.
            %  @param sigma0.
            %  @param fileprefix.

 			this = this@mlglucose.Huang1980(varargin{:});	
            
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

