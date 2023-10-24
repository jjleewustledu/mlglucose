classdef DispersedNumericHuang1980 < handle & mlpet.AugmentedData & mlglucose.Huang1980
	%% DISPERSEDNUMERICHUANG1980  

	%  $Revision$
 	%  was created 04-Jan-2021 20:38:13 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.9.0.1538559 (R2020b) Update 3 for MACI64.  Copyright 2021 John Joowon Lee.
 	 	
    
    properties (Constant)
        LENK = 7        
    end
    
    properties (Dependent)
        artery_sampled
        Delta
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
            %  @return tac_.
            %  @return aif_.
            
            import mlglucose.DispersedNumericHuang1980
            import mlkinetics.ScannerKit.mixTacAif
            import mlglucose.Huang1980
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'scanner', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'arterial', [], @(x) isa(x, 'mlpet.AbstractDevice') || iscell(x))          
            addParameter(ip, 'roi', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'cbv', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            parse(ip, devkit, varargin{:})
            ipr = ip.Results;

            if iscell(ipr.arterial)
                [this,tac_,aif_] = DispersedNumericHuang1980.createFromDualDeviceKit(devkit, varargin{:});
                return
            end
            
            [tac_,timesMid_,aif_,Dt] = mixTacAif(devkit, varargin{:});
            
            % v1              

            v1 = ipr.cbv * 0.0105;
            v1 = v1.volumeAveraged(ipr.roi);
            v1_ = v1.imagingFormat.img;
            
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
        function [this,tac_,aif_] = createFromDualDeviceKit(devkit, varargin)
            %% manage Twilite and Caprac

            import mlglucose.DispersedNumericHuang1980
            import mlglucose.DispersedNumericHuang1980.mixTacAifHybrid
            import mlglucose.Huang1980
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'scanner', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'arterial', [], @iscell)
            addParameter(ip, 'roi', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'cbv', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            parse(ip, devkit, varargin{:})
            ipr = ip.Results;

            [tac_,timesMid_,aif_,Dt] = mixTacAifHybrid(devkit, varargin{:});
            
            % v1              

            v1 = ipr.cbv * 0.0105;
            v1 = v1.volumeAveraged(ipr.roi);
            v1_ = v1.imagingFormat.img;
            
            %
            
            glc_ = Huang1980.glcFromRadMeasurements(ipr.arterial{1}.radMeasurements);
            hct_ = Huang1980.hctFromRadMeasurements(ipr.arterial{1}.radMeasurements);
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
        function [tac__,timesMid__,aif__,Dt,datetimePeak] = mixTacAifHybrid(devkit, varargin)
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'scanner', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'arterial', [], @iscell)
            addParameter(ip, 'roi', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            parse(ip, devkit, varargin{:})
            ipr = ip.Results;      
            ad = mlaif.AifData.instance();
            
            % scannerDevs provide calibrations & ROI-volume averaging            
            s = ipr.scanner.volumeAveraged(ipr.roi);
            tac = s.activityDensity();
            tac(tac < 0) = 0;                       
            tac = ad.normalizationFactor*tac; % empirical normalization
            tac__ = tac;
            timesMid__ = s.timesMid;
            Nt = ceil(timesMid__(end));
            
            % arterialDevs calibrate & align arterial times-series to localized scanner time-series  
            [a1, datetimePeak] = devkit.alignArterialToScanner( ...
                ipr.arterial{1}, s, 'sameWorldline', false);
            Dt = a1.Dt;
            aif1 = a1.activityDensity();
            aif2 = ipr.arterial{2}.activityDensityInterp1(); % 1 Hz interp1
            daif = aif1(end) - aif2(1);
            aif = [aif1, aif2+daif];
            if length(aif) > Nt
                aif = aif(1:Nt); % truncate late aif
            elseif length(aif) < Nt
                Nremain = Nt - length(aif);
                aif = [aif, aif(end)*ones(1,Nremain)]; % extrapolate late aif
            end
            t = (0:Nt-1) - seconds(s.datetime0 - a1.datetime0);
            
            % use tBuffer to increase fidelity of kinetic model
            while any(-ad.tBuffer == t)
                ad.T = ad.T + 1;
            end
            aif = interp1([-ad.tBuffer t], [0 aif], -ad.tBuffer:s.timesMid(end), 'linear', 0);
            aif(aif < 0) = 0;   
            aif__ = aif;
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
        function g = get.Dt(this)
            g = this.strategy_.Dt;
        end
        
        %%
        
		function a = artery_local_mediated(this, varargin)
            %% ARTERY_LOCAL returns artery activities mapped into R^(3+1), space-times,
            %  shifted by this.Dt and disperses by this.Delta
            %  @return a is an ImagingFormatContext2, the artery activities sampled on scanner space-times.
            %  See also ml*.Dispersed*Model.

            n = length(this.artery_interpolated);
            times = 0:1:n-1;            
            auc0 = trapz(this.artery_interpolated);
            artery_interpolated_ = conv(this.artery_interpolated, exp(-this.Delta*times));
            if this.Dt < 0 % shift back to right
                artery_interpolated1 = zeros(1, n);
                artery_interpolated1(-this.Dt+1:end) = artery_interpolated_(1:n+this.Dt);
            elseif this.Dt > 0 % shift back to left
                artery_interpolated1 = artery_interpolated_(this.Dt+1:this.Dt+n);
            else
                artery_interpolated1 = artery_interpolated_(1:n);
            end 
            artery_interpolated1 = artery_interpolated1*auc0/trapz(artery_interpolated1);
            artery_interpolated1 = artery_interpolated1(this.tBuffer+1:end);
            avec = this.model.solutionOnScannerFrames(artery_interpolated1, this.times_sampled);
            avec = single(asrow(avec));

            Nel = numel(this.roi);
            roi_vec = reshape(logical(this.roi), [Nel, 1]);
            Nroi = sum(roi_vec);
            Nt = length(avec);
            avec_mat = repmat(avec, Nroi, 1);
            
            img = zeros(Nel, Nt, 'single');
            img(roi_vec,:) = avec_mat;

            a = this.roi.imagingFormat;
            a.img = reshape(img, [size(this.roi) Nt]);
            a.fileprefix = strcat(stackstr(), '_a');
        end 
        function ks = buildKs(this, varargin)
            this = solve(this, varargin{:});
            ks = [k1(this) k2(this) k3(this) k4(this) k5(this) k6(this) loss(this)];
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
        function ks_ = ks_mediated(this, varargin)
            %% ks == [k1 k2 k3 k4 Delta Dt loss]
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
            k(7) = loss(this);
            
            Nel = numel(this.roi);
            roi_vec = reshape(logical(this.roi), [Nel, 1]);
            Nroi = sum(roi_vec);
            Nt = length(k);
            k_mat = repmat(k, Nroi, 1);
            
            img = zeros(Nel, Nt, 'single');
            img(roi_vec,:) = k_mat;

            ks_ = this.roi.imagingFormat;
            ks_.img = reshape(img, [size(this.roi) Nt]);
            ks_.fileprefix = this.ksOnAtlas('typ', 'fp', 'tags', [this.blurTag this.regionTag]);
        end
        function ks_ = ks(this, varargin)
            %% ks == [k1 k2 k3 k4 Delta Dt loss]
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
            k(7) = loss(this);
             
            roibin = logical(this.roi);
            ks_ = copy(this.roi.imagingFormat);
            ks_.img = zeros([size(this.roi) length(k)]);
            for t = 1:length(k)
                img = zeros(size(this.roi), 'single');
                img(roibin) = k(t);
                ks_.img(:,:,:,t) = img;
            end
            ks_.fileprefix = this.ksOnAtlas('typ', 'fp', 'tags', [this.blurTag this.regionTag]);
            ks_ = imagingType(ipr.typ, ks_);
        end
        function fp = ksOnAtlas(this, varargin)
            if isa (this.sessionData, 'mlnipet.SessionData')
                fp = this.sessionData.ksOnAtlas('typ', 'fp', 'tags', [this.blurTag this.regionTag]);
                return
            end
            if isa (this.sessionData, 'mlpipeline.ImagingMediator')
                ic = this.sessionData.metricOnAtlas('ks', [this.blurTag this.regionTag]);
                fp = ic.fileprefix;
                return
            end
            error('mloxygen:RuntimeError', stackstr())
        end
    end
    
    %% PROTECTED

	methods (Access = protected)
 		function this = DispersedNumericHuang1980(varargin)
 			%% DISPERSEDNUMERICHUANG1980
            %  @param fdg is numeric.
            %  @param solver is in {'nest' 'simulanneal' 'hmc' 'lm' 'bfgs'}. 
            %  @param devkit is mlpet.IDeviceKit.            
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
                    error('mlglucose:NotImplementedError', 'DispersedNumericHuang1980.ipr.solver->%s', ipr.solver)
            end
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

