classdef Huang1980SimulAnneal < mlpet.TracerSimulAnneal & mlglucose.Huang1980Strategy
	%% HUANG1980SIMULANNEAL operates on single voxels/regions.

	%  $Revision$
 	%  was created 03-Jan-2020 17:45:56 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.7.0.1261785 (R2019b) Update 3 for MACI64.  Copyright 2020 John Joowon Lee.
 	
    properties
        v1 % blood volume fraction
    end

	methods        
 		function this = Huang1980SimulAnneal(varargin)
 			%% HUANG1980SIMULANNEAL
            %  @param context is mlglucose.Huang1980.
            %  @param sigma0.
            %  @param fileprefix.
            
            this = this@mlpet.TracerSimulAnneal(varargin{:});
                      
            [this.ks_lower,this.ks_upper,this.ks0] = remapper(this);
            this.artery_interpolated = this.model.artery_interpolated;
            this.v1 = this.model.v1;
        end
        
        function fprintfModel(this)
            fprintf('Simulated Annealing:\n');
            for ky = 1:length(this.ks)
                fprintf('\tk%i = %f\n', ky, this.ks(ky));
            end
            fprintf('\tsigma0 = %f\n', this.sigma0);
            fprintf('\tv1 = %f\n', this.v1);
            for ky = this.map.keys
                fprintf('\tmap(''%s'') => %s\n', ky{1}, struct2str(this.map(ky{1})));
            end
        end
        function [k,sk] = K1(this, varargin)
            [k,sk] = k1(this, varargin{:});
            k = k*this.v1;
            sk = sk*this.v1;
        end
        function [k,sk] = k1(this, varargin)
            [k,sk] = find_result(this, 'k1');
        end
        function [k,sk] = k2(this, varargin)
            [k,sk] = find_result(this, 'k2');
        end
        function [k,sk] = k3(this, varargin)
            [k,sk] = find_result(this, 'k3');
        end    
        function [k,sk] = k4(this, varargin)
            [k,sk] = find_result(this, 'k4');
        end          
        function [k,sk] = k5(this, varargin)
            [k,sk] = find_result(this, 'k5');
        end  
        function h = plot(this, varargin)
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'showAif', true, @islogical)
            addParameter(ip, 'xlim', [-20 1800], @isnumeric)            
            addParameter(ip, 'ylim', [], @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            aif = this.artery_interpolated;
            h = figure;
            times = this.times_sampled;
            if ipr.showAif
                plot(times, this.Measurement, ':o', ...
                    times, this.model.sampled(this.ks, this.v1, aif, times), '-', ...
                    0:length(aif)-1, this.v1*aif, '--')                
                legend('measurement', 'estimation', 'v1*aif')
            else
                plot(times, this.Measurement, 'o', ...
                    times, this.model.sampled(this.ks, this.v1, aif, times), '-')                
                legend('measurement', 'estimation')
            end
            if ~isempty(ipr.xlim); xlim(ipr.xlim); end
            if ~isempty(ipr.ylim); ylim(ipr.ylim); end
            xlabel('times / s')
            ylabel('activity / (Bq/mL)')
            annotation('textbox', [.175 .25 .3 .3], 'String', sprintfModel(this), 'FitBoxToText', 'on', 'FontSize', 7, 'LineStyle', 'none')
            title('Huang1980SimulAnneal.plot()')
        end 
        function this = solve(this, varargin)
            ip = inputParser;
            addRequired(ip, 'loss_function', @(x) isa(x, 'function_handle'))
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            options_fmincon = optimoptions('fmincon', ...
                'FunctionTolerance', 1e-9, ...
                'OptimalityTolerance', 1e-9);
            if this.visualize_anneal
                options = optimoptions('simulannealbnd', ...
                    'AnnealingFcn', 'annealingboltz', ...
                    'FunctionTolerance', eps, ...
                    'HybridFcn', {@fmincon, options_fmincon}, ...
                    'InitialTemperature', 20, ...
                    'ReannealInterval', 200, ...
                    'TemperatureFcn', 'temperatureexp', ...
                    'Display', 'diagnose', ...
                    'PlotFcns', {@saplotbestx,@saplotbestf,@saplotx,@saplotf,@saplotstopping,@saplottemperature});
            else
                options = optimoptions('simulannealbnd', ...
                    'AnnealingFcn', 'annealingboltz', ...
                    'FunctionTolerance', eps, ...
                    'HybridFcn', {@fmincon, options_fmincon}, ...
                    'InitialTemperature', 20, ...
                    'ReannealInterval', 200, ...
                    'TemperatureFcn', 'temperatureexp');
            end
 			[ks_,sse,exitflag,output] = simulannealbnd( ...
                @(ks__) ipr.loss_function( ...
                       ks__, double(this.v1), this.artery_interpolated, this.times_sampled, double(this.Measurement), this.sigma0), ...
                this.ks0, this.ks_lower, this.ks_upper, options); 
            
            this.results_ = struct('ks0', this.ks0, 'ks', ks_, 'sse', sse, 'exitflag', exitflag, 'output', output); 
            if ~this.quiet
                fprintfModel(this)
            end
            if this.visualize
                plot(this)
            end
        end         
        function s    = sprintfModel(this)
            s = sprintf('Simulated Annealing:\n');
            for ky = 1:length(this.ks)
                s = [s sprintf('\tk%i = %f\n', ky, this.ks(ky))]; %#ok<AGROW>
            end
            s = [s sprintf('\tsigma0 = %f\n', this.sigma0)];
            s = [s sprintf('\tv1 = %f\n', this.v1)];
            for ky = this.map.keys
                s = [s sprintf('\tmap(''%s'') => %s\n', ky{1}, struct2str(this.map(ky{1})))]; %#ok<AGROW>
            end
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

