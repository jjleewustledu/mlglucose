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
            %  @param convert_wb2plasma, default from mlglucose.Huang1980Model.
            %  @param sigma0, default from mloptimization.SimulatedAnnealing.
            %  @param fileprefix, default from devkit.
            
            import mlfourd.ImagingContext2
            import mlglucose.Huang1980.glcFromRadMeasurements
            import mlglucose.Huang1980.hctFromRadMeasurements
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'fdg', [], @isnumeric)
            addParameter(ip, 'roi', 'brain.4dfp.hdr', @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'cbv', [])
            parse(ip, devkit, varargin{:})
            ipr = ip.Results;
            
            ipr.roi = mlfourd.ImagingContext2(ipr.roi);
            roibin = ipr.roi.binarized();            
            scanner = ipr.devkit.buildScannerDevice();
            scanner = scanner.blurred(4.3);
            scanner = scanner.volumeAveraged(roibin);
            fdg = asrow(scanner.imagingContext.fourdfp.img);
            fp = sprintf('mlglucose_Huang1980_createFromDeviceKit_dt%s', datestr(now, 'yyyymmddHHMMSS'));  
            if isnumeric(ipr.cbv)
                v1 = ipr.cbv/100;
            else
                v1 = mlfourd.ImagingContext2(ipr.cbv)/100;
                v1 = v1.volumeAveraged(roibin);
            end
            counting = ipr.devkit.buildCountingDevice();
            aif = pchip(counting.times, counting.activityDensity(), 0:scanner.times(end));
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
            %  @param convert_wb2plasma is logical.
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

