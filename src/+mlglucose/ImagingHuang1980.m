classdef ImagingHuang1980 < handle & matlab.mixin.Copyable 
	%% IMAGINGHUANG1980  

	%  $Revision$
 	%  was created 28-Apr-2020 23:53:00 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.7.0.1319299 (R2019b) Update 5 for MACI64.  Copyright 2020 John Joowon Lee.
 	 	
	properties (Constant)
        JITTER = 0.01 % > 0 to aid deep learning
        MAX_NORMAL_BACKGROUND = 20 % Bq/mL
        MIN_V1 = 0.001; % fraction
    end
    
	properties
        artery_plasma_interpolated
        blur
        fdg
        glc
        hct
        ks
        LC
        meanAbsError
        measurement % expose for performance when used by mlglucose.Huang1980Strategy
        model       %
        normMeanAbsError
        Nroi
        prediction
        regionTag
        residual
        roi
        roibin
        taus
 		times_sampled
        v1
    end
    
    methods (Static)
        function this = createFromDeviceKit(devkit, varargin)
            %% makes no adjustments of AIF timing
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param cbv is understood by mlfourd.ImagingContext2.
            %  @param roi is understood by mlfourd.ImagingContext2.
            %  @param blurFdg := {[], 0, 4.3, ...}
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'cbv', [], @(x) ~isempty(x))
            addParameter(ip, 'roi', [], @(x) ~isempty(x))
            addParameter(ip, 'blurFdg', 4.3, @isnumeric)
            parse(ip, devkit, varargin{:})
            ipr = ip.Results;
            
            scanner  = ipr.devkit.buildScannerDevice();
            fdg = scanner.imagingContext;
            fdg = fdg.blurred(ipr.blurFdg);
            
            counting = ipr.devkit.buildCountingDevice();
            aif = pchip(counting.times, counting.activityDensity(), 0:counting.times(end));
            radm = counting.radMeasurements;
            
            this = mlglucose.ImagingHuang1980( ...
                'fdg', fdg, ...
                'taus', scanner.taus, ...
                'times_sampled', scanner.timesMid, ...
                'artery_sampled', aif, ...
                'glc', mlglucose.Huang1980.glcFromRadMeasurements(radm), ...
                'hct', mlglucose.Huang1980.hctFromRadMeasurements(radm), ...
                varargin{:});
        end
    end
    
    methods
        function this = ensureModel(this)
            if isempty(this.model)                               
                map = mlglucose.Huang1980Model.preferredMap();
                this.model = mlglucose.Huang1980Model( ...
                    'map', map, ...
                    'times_sampled', this.times_sampled, ...
                    'artery_interpolated', this.artery_plasma_interpolated, ...
                    'glc', this.glc, ...
                    'hct', this.hct, ...
                    'LC', this.LC);
            end            
        end
        function [ic,nic] = buildMeanAbsError(this)
            % @return \Sigma_i activity_i \frac{tau_i}/{T}
            
            if isempty(this.residual)
                this.buildResidual()
            end
            assert(~isempty(this.taus))
            
            % MAE
            this.meanAbsError = abs(copy(this.residual));
            this.meanAbsError = this.meanAbsError.timeAveraged('taus', this.taus);
            this.meanAbsError.fileprefix = [this.fdg.fileprefix this.regionTag '_MAE'];
            ic = this.meanAbsError;
            ic.fileprefix = [fdgDCorr.fileprefix this.regionTag '_MAE'];
            
            % NMAE
            fdgTimeAverage = copy(this.fdg);
            fdgTimeAverage = fdgTimeAverage.timeAveraged();
            this.normMeanAbsError = copy(this.meanAbsError);            
            this.normMeanAbsError = this.meanAbsError ./ fdgTimeAverage;
            this.normMeanAbsError.fileprefix = [this.fdg.fileprefix this.regionTag '_NMAE'];
            nic = this.normMeanAbsError;
            nic.fileprefix = [fdgDCorr.fileprefix this.regionTag '_NMAE'];
        end
        function ic = buildPrediction(this, varargin)
            %% @param reuseExisting is logical; default is false.
            
            ip = inputParser;
            addParameter(ip, 'reuseExisting', false, @islogical)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            assert(~isempty(this.ks))
            
            % initialize ifc from this.fdg
            % return existing if reuseExisting
            fdg_ifc = copy(this.fdg.fourdfp);
            fdg_ifc.fileprefix = [fdg_ifc.fileprefix this.regionTag '_predicted'];
            if ipr.reuseExisting && isfile(fdg_ifc.fqfilename)
                ic = mlfourd.ImagingContext2(fdg_ifc.fqfilename);
            fileprefix0 = [working_ifc.fileprefix this.regionTag '_predicted'];
                this.prediction = ic;
                return
            end
            fdg_sz = size(fdg_ifc);
            fdg_ifc.img = zeros(fdg_sz);
            
            % represent prediction as voxels (x) times
            img2d = zeros(this.Nroi, fdg_sz(4));
            ks_ = zeros(this.Nroi,5);
            for ik = 1:5     
                rate = this.ks.fourdfp.img(:,:,:,ik);                
                ks_(:, ik) = rate(this.roibin);
            end
            v1_ = this.v1.fourdfp.img(this.roibin);
            this.ensureModel() % wihtout voxelwise adjustments of AIF timings
            for vxl = 1:this.Nroi                
                img2d(vxl,:) = this.model.simulated(ks_(vxl,:), 'v1', v1_(vxl));
            end
            
            % embed prediction in R^3 (x) times
            for t = 1:fdg_sz(4)
                frame = zeros(fdg_sz(1:3));
                frame(this.roibin) = img2d(:,t);
                fdg_ifc.img(:,:,:,t) = frame;
            end
            this.prediction = mlfourd.ImagingContext2(fdg_ifc);
            this.prediction = this.prediction.blurred(4.3);
            ic = this.prediction;
        end
        function ic = buildResidual(this)
            if isempty(this.prediction)
                this.buildPrediction()
            end
            this.residual = copy(this.fdg);
            this.residual = this.fdg - this.prediction;
            this.residual.fileprefix = [this.fdg.fileprefix this.regionTag '_residual'];
            ic = this.residual;
            ic.fileprefix = [fdgDCorr.fileprefix this.regionTag  '_residual'];
        end
        function this = solve(this)
            fdg_img_2d = this.projectedFdgArray();
            v1_img_1d = this.v1.fourdfp.img(this.roibin);
            ks_img_2d = zeros(dipsum(this.roibin), 4);
            map = mlglucose.Huang1980Model.preferredMap();
            
            for vxl = 1:this.Nroi
                try
                    tic
                    fprintf('mlglucose.ImagingHuang1980.solve():  vxl->%i this.Nroi->%i\n', vxl, this.Nroi)
                    this.measurement = fdg_img_2d(vxl, :);
                    the_v1 = v1_img_1d(vxl);
                    if this.sufficientData(the_v1, this.measurement)
                        this.model = mlglucose.Huang1980Model( ...
                            'map', map, ...
                            'times_sampled', this.times_sampled, ...
                            'v1', the_v1, ...
                            'artery_interpolated', this.artery_plasma_interpolated, ...
                            'glc', this.glc, ...
                            'hct', this.hct, ...
                            'LC', this.LC);
                        strategy = mlglucose.Huang1980SimulAnneal('context', this);
                        strategy = solve(strategy);
                        
                        % store latest solutions
                        ks_img_2d(vxl, 1) = k1(strategy);
                        ks_img_2d(vxl, 2) = k2(strategy);
                        ks_img_2d(vxl, 3) = k3(strategy);
                        ks_img_2d(vxl, 4) = k4(strategy);
                        
                        % use latest solutions for initial conditions for solving neighboring voxels
                        map_k1 = map('k1'); % cache mapped struct
                        map_k2 = map('k2');
                        map_k3 = map('k3');
                        map_k4 = map('k4');                        
                        map_k1.init = ks_img_2d(vxl, 1)*(1 + this.JITTER*randn()); % cached struct.init := latest solutions with jitter
                        map_k2.init = ks_img_2d(vxl, 2)*(1 + this.JITTER*randn());
                        map_k3.init = ks_img_2d(vxl, 3)*(1 + this.JITTER*randn());
                        map_k4.init = ks_img_2d(vxl, 4)*(1 + this.JITTER*randn());
                        map('k1') = map_k1; % update mapped struct with adjusted cache
                        map('k2') = map_k2;
                        map('k3') = map_k3;
                        map('k4') = map_k4; 
                    
                        % else ks_img_2d retains zeros                        
                    end                    
                    toc
                catch ME
                    %handexcept(ME)
                    handwarning(ME)
                end
            end
            
            this.ks = this.invProjectedKsArray(ks_img_2d);
        end
        function tf = sufficientData(this, v1, measurement)
            tf = v1 > this.MIN_V1 && mean(measurement) > this.MAX_NORMAL_BACKGROUND;
        end
    end

    %% PROTECTED
    
	methods (Access = protected)
        function this = ImagingHuang1980(varargin)
            
            %% 
            %  @param fdg is understood by mlfourd.ImagingContext2.
            %  @param times_sampled from scanner is numeric.
            %  @param artery_sampled from counter is numeric.
            %  @param cbv is understood by mlfourd.ImagingContext2.
            %  @param roi is understood by mlfourd.ImagingContext2.
            %  @param glc
            %  @param hct
            %  @param LC, default := 0.81.   
            %  @param regionTag, default := '_brain'
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'fdg', [])
            addParameter(ip, 'taus', [], @isnumeric)
            addParameter(ip, 'times_sampled', [], @isnumeric)
            addParameter(ip, 'artery_sampled', [], @isnumeric)
            addParameter(ip, 'cbv', [], @(x) ~isempty(x))
            addParameter(ip, 'roi', [], @(x) ~isempty(x))
            addParameter(ip, 'glc', @isnumeric)
            addParameter(ip, 'hct', @isnumeric)
            addParameter(ip, 'LC', 0.81, @isnumeric)
            addParameter(ip, 'blur', 4.3, @isnumeric)
            addParameter(ip, 'regionTag', '_brain', @ischar)
            parse(ip, varargin{:})
            addParameter(ip, 'regionTag', '', @ischar)
            ipr = ip.Results;
            
            this.blur = ipr.blur;
            this.fdg = mlfourd.ImagingContext2(ipr.fdg);
            assert(4 == ndims(this.fdg))
            this.taus = ipr.taus;
            this.times_sampled = ipr.times_sampled;  
            this.hct = ipr.hct;
            
            % artery management            
            t = 0:this.times_sampled(end);
            this.artery_plasma_interpolated = pchip(0:length(ipr.artery_sampled)-1, ipr.artery_sampled, t);
            
            cbvic = mlfourd.ImagingContext2(ipr.cbv);
            if dipisnan(cbvic)
                error('mlglucose:RuntimeError', 'ImagingHuang1980 found dipisnan(cbvic)')
            end
            this.v1 = cbvic/100;
            this.roi = mlfourd.ImagingContext2(ipr.roi);
            this.roibin = this.roi.fourdfp.img > 0;
            this.fdg = this.masked(this.fdg);
            this.Nroi = dipsum(this.roibin);
            this.glc = ipr.glc;
            this.LC = ipr.LC;
            this.regionTag = ipr.regionTag;
        end
        function that = copyElement(this)
            %%  See also web(fullfile(docroot, 'matlab/ref/matlab.mixin.copyable-class.html'))
            
            that = copyElement@matlab.mixin.Copyable(this);
        end
        function ic = invProjectedKsArray(this, arr)
            img = zeros([size(this.roibin) 4]);
            for ik = 1:4
                cube = img(:,:,:,ik);
                cube(this.roibin) = arr(:,ik); 
                img(:,:,:,ik) = cube;
            end
            fp = strrep(this.fdg.fileprefix, 'fdg', 'ks'); % sprintf('ks_%s_', this.roi.fileprefix));
            ic = mlfourd.ImagingContext2(img, 'filename', [fp '_' this.roi.fileprefix '.4dfp.hdr'], 'mmppix', [2 2 2]);
        end
        function ic = masked(this, ic)
            %% retains original fileprefix.
            
            fp = ic.fileprefix;
            ic = ic.masked(this.roi);
            ic.fileprefix = fp;
        end
        function arr = projectedFdgArray(this)
            Nt = size(this.fdg, 4);            
            arr = zeros(dipsum(this.roibin), Nt);
            for it = 1:Nt
                cube = this.fdg.fourdfp.img(:,:,:,it);
                arr(:,it) = cube(this.roibin);
            end
        end
  	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

