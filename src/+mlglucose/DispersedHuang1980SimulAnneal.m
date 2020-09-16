classdef DispersedHuang1980SimulAnneal < mlglucose.Huang1980SimulAnneal
	%% DISPERSEDHUANG1980SIMULANNEAL  

	%  $Revision$
 	%  was created 11-Aug-2020 20:44:37 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.7.0.1434023 (R2019b) Update 6 for MACI64.  Copyright 2020 John Joowon Lee.
    
    methods  (Static)
        function loss = loss_function(ks, v1, artery_interpolated, times_sampled, measurement, sigma0)
            import mlglucose.DispersedHuang1980Model.sampled            
            estimation  = sampled(ks, v1, artery_interpolated, times_sampled);
            positive    = measurement > 0.05*max(measurement);
            eoverm      = estimation(positive)./measurement(positive);
            Q           = sum((1 - eoverm).^2);
            loss        = 0.5*Q/sigma0^2; % + sum(log(sigma0*measurement)); % sigma ~ sigma0*measurement
        end 
    end

	methods		  
 		function this = DispersedHuang1980SimulAnneal(varargin)
 			this = this@mlglucose.Huang1980SimulAnneal(varargin{:});
        end
        
        function [k,sk] = k5(this, varargin)
            [k,sk] = find_result(this, 'k5');
        end 
        function this = solve(this, varargin)
            import mlglucose.DispersedHuang1980SimulAnneal.loss_function   
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
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

