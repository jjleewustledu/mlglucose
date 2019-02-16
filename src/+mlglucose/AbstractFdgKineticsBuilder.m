classdef (Abstract) AbstractFdgKineticsBuilder < handle & mlglucose.AbstractGlucoseKineticsBuilder
	%% ABSTRACTFDGKINETICSBUILDER  

	%  $Revision$
 	%  was created 18-Aug-2018 22:17:57 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	methods 
		  
 		function this = AbstractFdgKineticsBuilder(varargin)
 			%% ABSTRACTFDGKINETICSBUILDER
 			%  @param .

 			this = this@mlglucose.AbstractGlucoseKineticsBuilder(varargin{:});
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

