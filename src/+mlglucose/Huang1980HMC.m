classdef Huang1980HMC
	%% HUANG1980HMC  

	%  $Revision$
 	%  was created 10-Apr-2020 18:03:37 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.7.0.1319299 (R2019b) Update 5 for MACI64.  Copyright 2020 John Joowon Lee.
 	
	properties
 		
 	end

	methods 
		  
 		function this = Huang1980HMC(varargin)
 			%% HUANG1980HMC
 			%  @param .

 			this = this@mlglucose.Huang1980MCMC(varargin{:});
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

