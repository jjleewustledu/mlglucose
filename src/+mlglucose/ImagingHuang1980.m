classdef ImagingHuang1980 < handle & matlab.mixin.Copyable 
	%% IMAGINGHUANG1980
    %  builds Huang models for imaging voxels.
    %  For performance considerations, see also https://blogs.mathworks.com/loren/2012/03/26/considering-performance-in-object-oriented-matlab-code/

	%  $Revision$
 	%  was created 28-Apr-2020 23:53:00 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.7.0.1319299 (R2019b) Update 5 for MACI64.  Copyright 2020 John Joowon Lee.
 	 	
	properties (Constant)
        HALFLIFE = 6586.272; % s
        JITTER = 0.0 % > 0 to aid deep learning
        LENK = 4
        MAX_NORMAL_BACKGROUND = 20 % Bq/mL
        MIN_V1 = 0.001; % fraction
    end
    
	properties
        artery_plasma_interpolated
        devkit
        fdg % scanner.activityDensity(), decaying & calibrated
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
        residual
        roi
        roibin % is logical
        taus
 		times_sampled
        v1
    end
    
    properties (Dependent)        
        regionTag
        sessionData
    end
    
    methods (Static)
        function this = createFromDeviceKit(devkit, varargin)
            %% makes no adjustments of AIF timing
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param fdg is numeric, default from devkit.
            %  @param roi is understood by mlfourd.ImagingContext2.
            %  @param cbv is understood by mlfourd.ImagingContext2.
            %  @param blurFdg := {[], 0, 4.3, ...}
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'fdg', [], @isnumeric)
            addParameter(ip, 'roi', [], @(x) ~isempty(x))
            addParameter(ip, 'cbv', [], @(x) ~isempty(x))
            addParameter(ip, 'blurFdg', 4.3, @isnumeric)
            parse(ip, devkit, varargin{:})
            ipr = ip.Results;
            
            % scanner provides calibrations, ancillary data
            
            scanner = ipr.devkit.buildScannerDevice();
            scanner = scanner.blurred(ipr.blurFdg);
            fdg = scanner.activityDensity(); % calibrated, decaying
            
            % AIF  
            
            counting = ipr.devkit.buildCountingDevice();
            aif = pchip(counting.times, counting.activityDensity(), 0:scanner.times(end));
            radm = counting.radMeasurements;
            
            this = mlglucose.ImagingHuang1980( ...
                'devkit', devkit, ...
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
        
        %% GET
        
        function g = get.regionTag(this)
            g = this.sessionData.regionTag;
        end
        function g = get.sessionData(this)
            g = this.devkit.sessionData;
        end
        
        %%
        
        function [ic,nic] = buildMeanAbsError(this)
            % @return \Sigma_i activity_i \frac{tau_i}/{T}
            
            if isempty(this.residual)
                this.buildResidual()
            end
            assert(~isempty(this.taus))
            sesd = this.devkit.sessionData;
            fdgDCorr = sesd.fdgOnAtlas('typ', 'mlfourd.ImagingContext2');
            
            % MAE
            ic = abs(copy(this.residual));
            ic = ic.timeAveraged('taus', this.taus);
            ic.fileprefix = [fdgDCorr.fileprefix this.regionTag '_MAE'];
            this.meanAbsError = ic;
            
            % NMAE
            fdgAvgt = copy(fdgDCorr);
            fdgAvgt = fdgAvgt.timeAveraged('taus', this.taus);         
            nic = copy(this.meanAbsError) ./ fdgAvgt;
            nic.fileprefix = [fdgDCorr.fileprefix this.regionTag '_NMAE'];
            this.normMeanAbsError = nic;
        end
        function ic = buildPrediction(this, varargin)
            %% @param reuseExisting is logical; default is false.
            
            ip = inputParser;
            addParameter(ip, 'reuseExisting', false, @islogical)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            assert(~isempty(this.ks))
            sesd = this.devkit.sessionData;
            scanner = this.devkit.buildScannerDevice();
            
            % init working_ifc from fdgOnAtlas
            % return existing if reuseExisting
            working_ifc = sesd.fdgOnAtlas('typ', 'mlfourd.ImagingFormatContext');
            fileprefix0 = [working_ifc.fileprefix this.regionTag '_predicted'];
            working_ifc.fileprefix = fileprefix0;
            sz = size(working_ifc);
            if ipr.reuseExisting && isfile(working_ifc.fqfilename)
                ic = mlfourd.ImagingContext2(working_ifc.fqfilename);
                this.prediction = ic;
                return
            end            
            
            % represent prediction in R^2:  voxels^1 (x) times            
            ks_ = zeros(this.Nroi,this.LENK+1);
            for ik = 1:this.LENK+1     
                rate = this.ks.fourdfp.img(:,:,:,ik);                
                ks_(:, ik) = rate(this.roibin);
            end
            v1_ = this.v1.blurred(4.3);
            v1_ = v1_.fourdfp.img(this.roibin);            
            this.ensureModel() % without voxelwise adjustments of AIF timings  
            img2d = zeros(this.Nroi, sz(4));          
            for vxl = 1:this.Nroi                
                img2d(vxl,:) = this.model.simulated(ks_(vxl,:), 'v1', v1_(vxl)); % adjust AIF timings
            end
            
            % embed prediction in R^4:  voxels^3 (x) times
            working_ifc.img = zeros(sz);
            for t = 1:sz(4)
                frame = zeros(sz(1:3));
                frame(this.roibin) = img2d(:,t);
                working_ifc.img(:,:,:,t) = frame;
            end
            
            % decay-correct & remove calibrations
            ic = scanner.decayCorrectLike(working_ifc);
            ic = ic ./ scanner.invEfficiencyf(sesd);
            ic.fileprefix = fileprefix0;
            this.prediction = ic;
        end
        function ic = buildResidual(this)
            if isempty(this.prediction)
                this.buildPrediction()
            end            
            sesd = this.devkit.sessionData;
            fdgDCorr = sesd.fdgOnAtlas('typ', 'mlfourd.ImagingContext2');
            ic = this.prediction - fdgDCorr;
            ic.fileprefix = [fdgDCorr.fileprefix this.regionTag  '_residual'];
            this.residual = ic;
        end
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
            %% IMAGINGHUANG1980
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param fdg is understood by mlfourd.ImagingContext2.
            %  @param times_sampled from scanner is numeric.
            %  @param artery_sampled from counter is numeric.
            %  @param cbv is understood by mlfourd.ImagingContext2.
            %  @param roi is understood by mlfourd.ImagingContext2.
            %  @param glc
            %  @param hct
            %  @param LC, default := 0.81.  
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'fdg', [])
            addParameter(ip, 'taus', [], @isnumeric)
            addParameter(ip, 'times_sampled', [], @isnumeric)
            addParameter(ip, 'artery_sampled', [], @isnumeric)
            addParameter(ip, 'cbv', [], @(x) ~isempty(x))
            addParameter(ip, 'roi', [], @(x) ~isempty(x))
            addParameter(ip, 'glc', @isnumeric)
            addParameter(ip, 'hct', @isnumeric)
            addParameter(ip, 'LC', 0.81, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            this.devkit = ipr.devkit;
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
            this.v1 = cbvic ./ 105;
            this.roi = mlfourd.ImagingContext2(ipr.roi);
            this.roibin = this.roi.fourdfp.img > 0;
            this.fdg = this.masked(this.fdg);
            this.Nroi = dipsum(this.roibin);
            this.glc = ipr.glc;
            this.LC = ipr.LC;
        end
        function that = copyElement(this)
            %%  See also web(fullfile(docroot, 'matlab/ref/matlab.mixin.copyable-class.html'))
            
            that = copyElement@matlab.mixin.Copyable(this);
        end
        function ic = invProjectedKsArray(this, arr)
            img = zeros([size(this.roibin) this.LENK]);
            for ik = 1:this.LENK
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

