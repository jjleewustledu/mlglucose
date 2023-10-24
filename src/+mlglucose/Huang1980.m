classdef Huang1980 < handle & mlpet.TracerKinetics
	%% HUANG1980 is the context to a strategy design pattern which implements:
    %  mlglucose.{Huang1980Nest, Huang1980SimulAnneal, Huang1980HMC, Huang1980LM, Huang1980BFGS}.
    %  For performance considerations, see also https://blogs.mathworks.com/loren/2012/03/26/considering-performance-in-object-oriented-matlab-code/

	%  $Revision$
 	%  was created 10-Apr-2020 15:14:24 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.7.0.1319299 (R2019b) Update 5 for MACI64.  Copyright 2020 John Joowon Lee.
    
    properties (Dependent)
        v1
    end
    
    methods (Static)
        function this = createFromDeviceKit(devkit, varargin)
            %% 
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param Dt is numeric, s of time-shifting for AIF.
            %  @param fdg is numeric, default from devkit.
            %  @param cbv is numeric.
            %  @param solver is in {'nest' 'simulanneal' 'hmc' 'lm' 'bfgs'}, default := 'simulanneal'.
            %  @param roi ...
            %  @param map, default := mlglucose.Huang1980Model.preferredMap().
            %  @param times_sampled non-uniformly scheduled by the time-resolved PET reconstruction.
            %  @param artery_interpolated, default from devkit.
            %  @param glc is numeric, default from devkit.
            %  @param hct is numeric, default from devkit.
            %  @param LC is numeric, default from mlglucose.Huang1980Model.
            %  @param sigma0, default from mloptimization.SimulatedAnnealing.
            %  @param fileprefix, default from devkit.
            %  @return this.
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'fdg', [])
            addParameter(ip, 'cbv', [])
            parse(ip, devkit, varargin{:})
            ipr = ip.Results;
            
            if isnumeric(ipr.fdg) && isnumeric(ipr.cbv)
                this = mlglucose.NumericHuang1980.createFromDeviceKit(devkit, varargin{:});
            else
                this = mlglucose.ImagingHuang1980.createFromDeviceKit(devkit, varargin{:});
            end
        end
        function g = glcFromRadMeasurements(radm)
            %  @return mg/dL
            
            tbl = radm.laboratory;
            rows = tbl.Properties.RowNames;
            select = contains(rows, 'glc');
            g = tbl.measurement(select);
            g = mean(g(find(g)), 'omitnan'); %#ok<FNDSB>
        end
        function g = glcConversion(g, unitsIn, unitsOut)
            %  @param required g is numeric
            %  @param required unitsIn, unitsOut in {'mg/dL' 'mmol/L' 'umol/hg'}
            
            assert(isnumeric(g))
            assert(ischar(unitsIn))
            assert(ischar(unitsOut))
            if strcmp(unitsIn, unitsOut)
                return
            end
            
            switch unitsIn % to SI
                case 'mg/dL'
                    g = g * 0.0555;
                case 'mmol/L'
                case 'umol/hg'
                    % [mmol/L] == [umol/hg] [mmol/umol] [hg/g] [g/mL] [mL/L] 
                    g = g * 1e-3 * 1e-2 * 1.05 * 1e3;
                otherwise
                    error('mlglucose:ValueError', 'Huang1980.gclConversion')
            end
            
            switch unitsOut % SI to desired
                case 'mg/dL'
                    g = g / 0.0555;
                case 'mmol/L'
                case 'umol/hg'
                    % [umol/hg] == [mmol/L] [umol/mmol] [L/mL] [mL/g] [g/hg] 
                    g = g * 1e3 * 1e-3 * (1/1.05) * 1e2;
                otherwise
                    error('mlglucose:ValueError', 'Huang1980.gclConversion')
            end
        end
        function h = hctFromRadMeasurements(radm)
            h = radm.laboratory{'Hct', 'measurement'};
            if h > 1
                h = h/100;
            end
        end
    end
    
	methods  
        
        %% GET
        
        function g = get.v1(this)
            g = this.strategy_.v1;
        end
        
        %%
        
        function ks = buildKs(this, varargin)
            this = solve(this, varargin{:});
            ks = [k1(this) k2(this) k3(this) k4(this)];
        end        
        function c = chi(this, varargin)
            %  @return 1/s
            
            c = this.model.v1*k1(this, varargin{:})*k3(this, varargin{:})/ ...
                (k2(this, varargin{:}) + k3(this, varargin{:}));
        end
        function r = cmrglc(this, varargin)
            %  @return umol/hg/min
            
            % [umol/mmol] [(mmol/L) / (mg/dL)] [L/dL] [dL/mL] [g/hg] [mL/g] == [umol/hg]
            glc_ = this.glcConversion(this.model.glc, 'mg/dL', 'umol/hg'); 
            r = 60*this.chi(varargin{:})*glc_/this.model.LC;
        end
        function t = ctxglc(~, varargin)
            t = nan;
        end
        function f = freeglc(~, varargin)
            f = nan;
        end
        function [K,sK] = K1(this, varargin)
            [K,sK] = K1(this.strategy_, varargin{:});
        end
        function [k,sk] = k1(this, varargin)
            [k,sk] = k1(this.strategy_, varargin{:});
        end
        function [k,sk] = k2(this, varargin)
            [k,sk] = k2(this.strategy_, varargin{:});
        end
        function [k,sk] = k3(this, varargin)
            [k,sk] = k3(this.strategy_, varargin{:});
        end
        function [k,sk] = k4(this, varargin)
            [k,sk] = k4(this.strategy_, varargin{:});
        end
        function [K,sK] = Ks(this, varargin)
            K = zeros(1,4);
            sK = zeros(1,4);
            [K(1),sK(1)] = K1(this.strategy_, varargin{:});
            [K(2),sK(2)] = k2(this.strategy_, varargin{:});
            [K(3),sK(3)] = k3(this.strategy_, varargin{:});
            [K(4),sK(4)] = k4(this.strategy_, varargin{:});
        end
        function [k,sk] = ks(this, varargin)
            k = zeros(1,4);
            sk = zeros(1,4);
            [k(1),sk(1)] = k1(this.strategy_, varargin{:});
            [k(2),sk(2)] = k2(this.strategy_, varargin{:});
            [k(3),sk(3)] = k3(this.strategy_, varargin{:});
            [k(4),sk(4)] = k4(this.strategy_, varargin{:});
        end
    end
    
    %% PROTECTED
    
    methods (Access = protected)
        function this = Huang1980(varargin)  
            this = this@mlpet.TracerKinetics(varargin{:});
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

