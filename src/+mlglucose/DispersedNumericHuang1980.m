classdef DispersedNumericHuang1980 < handle & mlpet.AugmentedData & mlglucose.Huang1980
	%% DISPERSEDNUMERICHUANG1980  

	%  $Revision$
 	%  was created 04-Jan-2021 20:38:13 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.9.0.1538559 (R2020b) Update 3 for MACI64.  Copyright 2021 John Joowon Lee.
 	 	
    
    properties (Constant)
        LENK = 5        
    end
    
    properties (Dependent)
        artery_sampled
        Delta
    end
    
    properties 
        Dt % time-shift for AIF; Dt < 0 shifts backwards in time.
    end
    
    methods (Static)
        function [this,tac_,aif_] = createFromDeviceKit(devkit, varargin)
            %% adjusts AIF timings for coincidence of inflow with tissue activity from scanner
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param required scanner is an mlpet.AbstractDevice.
            %  @param required arterial is an mlpet.AbstractDevice.
            %  @param solver is in {'nest' 'simulanneal' 'hmc' 'lm' 'bfgs'}, default := 'simulanneal'.
            %  @param roi is mlfourd.ImagingContext2.
            %  @param cbv is mlfourd.ImagingContext2.
            %  @param LC is numeric, default from mlglucose.DispersedHuang1980Model.
            %  @param sigma0, default from mloptimization.SimulatedAnnealing.
            %  @param fileprefix, default from devkit.
            %  @return this.
            %  @return tac_, blurred by ipr.blurFdg.
            %  @return aif_.
            
            import mlglucose.DispersedNumericHuang1980
            import mlpet.AugmentedData.mixTacAif
            import mlpet.AugmentedData.mix
            import mlglucose.Huang1980
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'scanner', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'arterial', [], @(x) isa(x, 'mlpet.AbstractDevice'))          
            addParameter(ip, 'roi', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'cbv', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            parse(ip, devkit, varargin{:})
            ipr = ip.Results;
            
            [tac_,timesMid_,aif_,Dt] = mixTacAif(devkit, varargin{:});
            
            % v1              

            v1 = ipr.cbv * 0.0105;
            v1 = v1.volumeAveraged(ipr.roi);
            v1_ = v1.fourdfp.img;
            
            %
            
            glc_ = Huang1980.glcFromRadMeasurements(ipr.arterial.radMeasurements);
            hct_ = Huang1980.hctFromRadMeasurements(ipr.arterial.radMeasurements);
            fp = sprintf('mlglucose_DispersedNumericHuang1980_createFromDeviceKit_dt%s', datestr(now, 'yyyymmddHHMMSS'));    
            
            this = DispersedNumericHuang1980( ...
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
    end
    
    methods
        
        %% GET
        
        function g = get.artery_sampled(this)
            g = this.model.solutionOnScannerFrames( ...
                this.artery_interpolated(this.tBuffer+1:end), this.times_sampled);
        end
        function g = get.Delta(this)
            g = this.k5();
        end
        
        %%
        
        function ks = buildKs(this, varargin)
            this = solve(this, varargin{:});
            ks = [k1(this) k2(this) k3(this) k4(this) k5(this) k6(this)];
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
        function [k,sk] = k5(this, varargin)
            [k,sk] = k5(this.strategy_, varargin{:});
        end
        function [k,sk] = k6(this, varargin)
            k = this.Dt;
            sk = nan;
        end
        function ks_ = ks(this, varargin)
            %% ks == [k1 k2 k3 k4 Delta Dt]
            %  @param 'typ' is char, understood by imagingType.            
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'typ', 'mlfourd.ImagingContext2', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            k(1) = k1(this.strategy_, varargin{:});
            k(2) = k2(this.strategy_, varargin{:});
            k(3) = k3(this.strategy_, varargin{:});
            k(4) = k4(this.strategy_, varargin{:});
            k(5) = k5(this.strategy_, varargin{:});
            k(6) = k6(this);
             
            roibin = logical(this.roi);
            ks_ = copy(this.roi.fourdfp);
            ks_.img = zeros([size(this.roi) length(k)]);
            for t = 1:length(k)
                img = zeros(size(this.roi), 'single');
                img(roibin) = k(t);
                ks_.img(:,:,:,t) = img;
            end
            ks_.fileprefix = this.sessionData.ksOnAtlas('typ', 'fp', 'tags', [this.blurTag this.regionTag]);
            ks_ = imagingType(ipr.typ, ks_);
        end
    end
    
    %% PROTECTED

	methods (Access = protected)
 		function this = DispersedNumericHuang1980(varargin)
 			%% DISPERSEDNUMERICHUANG1980
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

 			this = this@mlglucose.Huang1980( ...
                'model', mlglucose.DispersedHuang1980Model(varargin{:}), ...
                varargin{:});	
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'fdg', [], @(x) isnumeric(x))
            addParameter(ip, 'solver', 'simulanneal', @ischar)
            addParameter(ip, 'Dt', 0, @isscalar)
            parse(ip, varargin{:})
            ipr = ip.Results;
                        
            this.measurement = ipr.fdg;
            switch lower(ipr.solver)
                case 'simulanneal'
                    this.strategy_ = mlglucose.DispersedHuang1980SimulAnneal( ...
                        'context', this, varargin{:});
                otherwise
                    error('mlglucose:NotImplementedError', 'DispersedNumericHuang1980.ipr.solver->%s', ipr.solver)
            end
            this.Dt = ipr.Dt;
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

