classdef (Abstract) AbstractGlucoseKineticsBuilder < handle & mlkinetics.AbstractKineticsBuilder
	%% ABSTRACTGLUCOSEKINETICSBUILDER  

	%  $Revision$
 	%  was created 17-Aug-2018 21:21:30 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	methods
		  
 		function this = AbstractGlucoseKineticsBuilder(varargin)
 			%% ABSTRACTGLUCOSEKINETICSBUILDER
            
            this = this@mlkinetics.AbstractKineticsBuilder(varargin{:});
 		end
    end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

