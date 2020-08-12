classdef Huang1980SimulAnneal < mloptimization.SimulatedAnnealing & mlglucose.Huang1980Strategy
	%% HUANG1980SIMULANNEAL operates on single voxels.

	%  $Revision$
 	%  was created 03-Jan-2020 17:45:56 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.7.0.1261785 (R2019b) Update 3 for MACI64.  Copyright 2020 John Joowon Lee.
 	
    properties
        artery_interpolated
        ks0
        ks_lower
        ks_upper
        quiet = false
        v1                   % blood volume fraction
        visualize = false
        visualize_anneal = false
    end
    
	properties (Dependent)   
        ks
        results
 	end
    
    methods (Static)
        function loss = loss_function(ks, v1, artery_interpolated, times_sampled, measurement, sigma0)
            import mlglucose.Huang1980Model.sampled            
            estimation  = sampled(ks, v1, artery_interpolated, times_sampled);
            positive    = measurement > 0;
            measurement = measurement(positive);
            eoverm      = estimation(positive)./measurement;
            Q           = sum((1 - eoverm).^2);
            loss        = 0.5*Q/sigma0^2; % + sum(log(sigma0*measurement)); % sigma ~ sigma0*measurement
        end   
        function conc = slide_fast(conc, Dt)
            %% SLIDE_FAST slides discretized function conc(t) to conc(t - Dt);
            %  @param conc is row vector without NaN.
            %  @param t is row vector with same size as conc.
            %  @param Dt is scalar rounded to integer.
            %
            %  Dt > 0 will slide conc(t) towards later times t.
            %  Dt < 0 will slide conc(t) towards earlier times t.
            
            Dt = round(Dt);
            if Dt == 0
                return
            end
            if Dt < 0
                T = length(conc);
               conc_ = conc(end)*ones(1, length(conc));
               conc_(1:T+Dt) = conc(1-Dt:end);
               conc = conc_;
               return
            end
            conc_ = zeros(size(conc));
            conc_(1+Dt:end) = conc(1:end-Dt);
            conc = conc_;
        end
    end

	methods 
        
        %% GET
        
        function g = get.ks(this)
            g = this.results_.ks;
        end
        function g = get.results(this)
            g = this.results_;
        end
        
        %%
        
 		function this = Huang1980SimulAnneal(varargin)
 			%% HUANG1980SIMULANNEAL
            %  @param context is mlglucose.Huang1980.
            %  @param sigma0.
            %  @param fileprefix.
            
            this = this@mloptimization.SimulatedAnnealing(varargin{:});
                      
            [this.ks_lower,this.ks_upper,this.ks0] = remapper(this);
            this.artery_interpolated = this.model.artery_interpolated;
            this.v1 = this.model.v1;
        end
        
        function disp(this)
            fprintf('\n')
            fprintf(class(this))
            fprintf('         v1: '); disp(this.v1)
            if isempty(this.results_)
                return
            end
            fprintf('initial ks0: '); disp(this.results_.ks0)
            fprintf('est.     ks: '); disp(this.results_.ks)
            fprintf('        sse: '); disp(this.results_.sse)
            fprintf('   exitflag: '); disp(this.results_.exitflag)
            disp(this.results_.output)
            disp(this.results_.output.rngstate)
            disp(this.results_.output.temperature)
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
        function h = plot(this, varargin)
            ip = inputParser;
            addParameter(ip, 'showAif', true, @islogical)
            addParameter(ip, 'xlim', [], @isnumeric)            
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
                legend('measurement', 'estimation', 'aif')
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
        function        save(this)
            save([this.fileprefix '.mat'], this);
        end
        function        saveas(this, fn)
            save(fn, this);
        end
        function this = solve(this, varargin)
            import mlglucose.Huang1980SimulAnneal.loss_function   
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
                @(ks__) loss_function( ...
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
    
    %% PROTECTED
    
    properties (Access = protected)        
        results_
    end
    
    methods (Access = protected)
        function [m,sd] = find_result(this, lbl)
            ks_ = this.ks;
            assert(strcmp(lbl(1), 'k'))
            ik = str2double(lbl(2));
            m = ks_(ik);
            sd = 0;
        end
        function [lb,ub,ks0] = remapper(this)
            for i = 1:this.map.Count
                lbl = sprintf('k%i', i);
                lb(i)  = this.map(lbl).min; %#ok<AGROW>
                ub(i)  = this.map(lbl).max; %#ok<AGROW>
                ks0(i) = this.map(lbl).init; %#ok<AGROW>
            end
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

