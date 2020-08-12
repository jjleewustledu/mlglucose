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
            positive    = measurement > 0;
            measurement = measurement(positive);
            eoverm      = estimation(positive)./measurement;
            Q           = sum((1 - eoverm).^2);
            loss        = 0.5*Q/sigma0^2; % + sum(log(sigma0*measurement)); % sigma ~ sigma0*measurement
        end 
    end

	methods		  
 		function this = DispersedHuang1980SimulAnneal(varargin)
 			%% DISPERSEDHUANG1980SIMULANNEAL
 			%  @param .

 			this = this@mlglucose.Huang1980SimulAnneal(varargin{:});
        end
        
        function [k,sk] = k5(this, varargin)
            [k,sk] = find_result(this, 'k5');
        end   
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

