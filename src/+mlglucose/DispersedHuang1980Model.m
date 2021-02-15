classdef DispersedHuang1980Model < mlpet.TracerKineticsModel
	%% DISPERSEDHUANG1980MODEL provides model data and methods to the strategy design pattern comprising
    %  mlglucose.{Huang1980, DispersedHuang1980SimulAnneal}.
    %  It operates on single voxels or regions.  It includes a dispersion parameter $\Delta$:
    %  aif_\text{disp}(t) = aif(t) \otimes e^{-\Delta t}.

	%  $Revision$
 	%  was created 10-Apr-2020 17:55:02 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.7.0.1319299 (R2019b) Update 5 for MACI64.  Copyright 2020 John Joowon Lee.
 		
	properties
        glc
        hct
        LC
        v1
    end
    
    methods (Static) 
        function loss = loss_function(ks, v1, artery_interpolated, times_sampled, measurement, ~)
            import mlglucose.DispersedHuang1980Model.sampled            
            estimation  = sampled(ks, v1, artery_interpolated, times_sampled);
            measurement = measurement(1:length(estimation));
            positive    = measurement > 0.05*max(measurement);
            eoverm      = estimation(positive)./measurement(positive);            
            Q           = mean(abs(1 - eoverm));
            %Q           = sum((1 - eoverm).^2);
            loss        = 0.5*Q; %/sigma0^2; % + sum(log(sigma0*measurement)); % sigma ~ sigma0*measurement
        end
        function m    = preferredMap()
            %% init from Huang's table 1
            m = containers.Map;
            m('k1') = struct('min', eps,  'max',  0.5,   'init', 0.048,   'sigma', 0.0048);
            m('k2') = struct('min', eps,  'max',  0.02,  'init', 0.0022,  'sigma', 0.0022);
            m('k3') = struct('min', eps,  'max',  0.01,  'init', 0.001,   'sigma', 0.0001);
            m('k4') = struct('min', eps,  'max',  0.001, 'init', 0.00011, 'sigma', 0.00011);
            m('k5') = struct('min', 0.01, 'max',  1,     'init', 1,       'sigma', 0.05);
        end
        function qs   = sampled(ks, v1, artery_interpolated, times_sampled)
            %  @param artery_interpolated is uniformly sampled at high sampling freq.
            %  @param times_sampled are samples scheduled by the time-resolved PET reconstruction
            
            import mlglucose.DispersedHuang1980Model.solution 
            import mlglucose.DispersedHuang1980Model.solutionOnScannerFrames  
            qs = solution(ks, v1, artery_interpolated);
            qs = solutionOnScannerFrames(qs, times_sampled);
        end
        function loss = simulanneal_objective(ks, v1, artery_interpolated, times_sampled, qs0, sigma0)
            import mlglucose.DispersedHuang1980Model.sampled          
            qs = sampled(ks, v1, artery_interpolated, times_sampled);            
            loss = 0.5 * sum((1 - qs ./ qs0).^2) / sigma0^2; % + sum(log(sigma0*qs0)); % sigma ~ sigma0 * qs0
        end 
        function qs   = solution(ks, v1, artery_interpolated)
            %  @param artery_interpolated is uniformly sampled at high sampling freq. starting at time = 0.
            
            k1 = ks(1);
            k2 = ks(2);
            k3 = ks(3);
            k4 = ks(4);
            Delta = ks(5);
            scale = 1;            
            n = length(artery_interpolated);
            times = 0:1:n-1;
            
            % use Delta
            auc0 = trapz(artery_interpolated);
            artery_interpolated1 = conv(artery_interpolated, exp(-Delta*times));
            artery_interpolated1 = artery_interpolated1(1:n);
            artery_interpolated1 = artery_interpolated1*auc0/trapz(artery_interpolated1);
            
            % use k1:k4
            k234 = k2 + k3 + k4;         
            bminusa = sqrt(k234^2 - 4 * k2 * k4);
            alpha = 0.5 * (k234 - bminusa);
            beta  = 0.5 * (k234 + bminusa);   
            conva = conv(exp(-alpha .* times), artery_interpolated1);
            convb = conv(exp(-beta .* times), artery_interpolated1);
            conva = conva(1:n);
            convb = convb(1:n);
            conv2 = (k4 - alpha) .* conva + (beta - k4) .* convb;
            conv3 =                 conva -                convb;
            q2 = (k1 / bminusa)      * conv2;
            q3 = (k3 * k1 / bminusa) * conv3;
            qs = v1 * (artery_interpolated1 + scale * (q2 + q3));            
        end 
    end

	methods		  
 		function this = DispersedHuang1980Model(varargin)
            %  @param v1
            %  @param glc
            %  @param hct
            %  @param LC
                
            this = this@mlpet.TracerKineticsModel(varargin{:});
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'v1', 0.038, @isnumeric)
            addParameter(ip, 'glc', 100, @isnumeric)
            addParameter(ip, 'hct', 0.4, @isnumeric)
            addParameter(ip, 'LC', 0.81, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            this.v1 = ipr.v1;
            this.glc = ipr.glc;
            this.hct = ipr.hct;            
            if (this.hct > 1)
                this.hct = this.hct/100;
            end
            this.LC = ipr.LC;         
        end
        function fdg  = simulated(this, varargin)
            %% SIMULATED simulates tissue activity with passed and internal parameters.
            %  @param required ks is [k1 k2 k3 k4 k5 Dt].
            %  @param v1 is CBV < 1 and dimensionless; default is this.v1.
            %  @param aif is numeric; default is this.artery_interpolated for model state.
            %  @param Dt is numeric, in sec.
        
            ip = inputParser;
            addRequired(ip, 'ks', @isnumeric)
            addParameter(ip, 'v1', this.v1, @isscalar)
            addParameter(ip, 'aif', this.artery_interpolated, @isnumeric)
            addParameter(ip, 'Dt', 0, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;            
            
            ks = ipr.ks(1:5);
            if length(ipr.ks) > 5
                ipr.Dt = ipr.ks(6);
            end
            if ipr.Dt ~= 0
                times = 0:length(ipr.aif)-1;
                aif = pchip(times - ipr.Dt, ipr.aif, times); % remove the delay Dt found by model
            else
                aif = ipr.aif;
            end
            fdg = mlglucose.DispersedHuang1980Model.sampled(ks, ipr.v1, aif, this.times_sampled);
        end
    end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

