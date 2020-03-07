classdef (Abstract) MetabBuilder < handle & mlkinetics.MetabBuilder
	%% METABBUILDER  

	%  $Revision$
 	%  was created 17-Dec-2018 00:15:36 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties
 		
 	end

	methods 
		  
 		function this = MetabBuilder(varargin)
 			%% METABBUILDER
 			%  @param .
 			
            this = this@mlkinetics.MetabBuilder(varargin{:});
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

