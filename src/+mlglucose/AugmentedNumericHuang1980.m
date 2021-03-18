classdef AugmentedNumericHuang1980 < handle & mlglucose.DispersedNumericHuang1980
	%% AUGMENTEDNUMERICHUANG1980  

	%  $Revision$
 	%  was created 04-Jan-2021 20:38:13 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.9.0.1538559 (R2020b) Update 3 for MACI64.  Copyright 2021 John Joowon Lee.

    
    methods (Static)
        function [this,tac_,aif_] = createFromDeviceKit(devkit, devkit2, varargin)
            %% adjusts AIF timings for coincidence of inflow with tissue activity from scanner
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param required devkit2 is mlpet.IDeviceKit.
            %  @param required scanner is an mlpet.AbstractDevice.
            %  @param required scanner2 is an mlpet.AbstractDevice.
            %  @param required arterial is an mlpet.AbstractDevice.
            %  @param required arterial2 is an mlpet.AbstractDevice.
            %  @param solver is in {'nest' 'simulanneal' 'hmc' 'lm' 'bfgs'}, default := 'simulanneal'.
            %  @param roi is mlfourd.ImagingContext2.
            %  @param roi2 is mlfourd.ImagingContext2.
            %  @param cbv is mlfourd.ImagingContext2.
            %  @param cbv2 is mlfourd.ImagingContext2. 
            %  @param DtMixing isscalar.
            %  @param LC is numeric, default from mlglucose.DispersedHuang1980Model.
            %  @param sigma0, default from mloptimization.SimulatedAnnealing.
            %  @param fileprefix, default from devkit.
            %  @param fracMixing in [0 1] for mixing tacs and aifs.
            %  @return this.
            %  @return tac_, blurred by ipr.blurFdg.
            %  @return aif_.
            
            import mlglucose.AugmentedNumericHuang1980
            import mlpet.AugmentedData.mixTacsAifs
            import mlpet.AugmentedData.mix
            import mlglucose.Huang1980
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addRequired(ip, 'devkit2', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'scanner', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'scanner2', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'arterial', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'arterial2', [], @(x) isa(x, 'mlpet.AbstractDevice'))            
            addParameter(ip, 'roi', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'roi2', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'cbv', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'cbv2', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'DtMixing', 0, @isscalar)
            addParameter(ip, 'fracMixing', 0.9, @isscalar)
            parse(ip, devkit, devkit2, varargin{:})
            ipr = ip.Results;
            
            % mix components for augmentation
            [tac_,timesMid_,aif_,Dt] = mixTacsAifs(devkit, devkit2, varargin{:});
            
            % v1              

            v1 = ipr.cbv * 0.0105;
            v1 = v1.volumeAveraged(ipr.roi);
            v1 = v1.fourdfp.img;
            v12 = ipr.cbv2 * 0.0105;
            v12 = v12.volumeAveraged(ipr.roi2);
            v12 = v12.fourdfp.img;
            v1_ = mix(v1, v12, ipr.fracMixing);
            
            % glc, hct
            
            glc = Huang1980.glcFromRadMeasurements(ipr.arterial.radMeasurements);
            glc2 = Huang1980.glcFromRadMeasurements(ipr.arterial2.radMeasurements);
            glc_ = mix(glc, glc2, ipr.fracMixing);
            hct = Huang1980.hctFromRadMeasurements(ipr.arterial.radMeasurements);
            hct2 = Huang1980.hctFromRadMeasurements(ipr.arterial2.radMeasurements);
            hct_ = mix(hct, hct2, ipr.fracMixing);
            
            %
            
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
            this.DtMixing = ipr.DtMixing;
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

 			this = this@mlglucose.DispersedNumericHuang1980(varargin{:});
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

