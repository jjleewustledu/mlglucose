classdef F18DeoxyGlucoseKineticsBuilder < mlglucose.GlucoseKineticsBuilder
	%% F18DEOXYGLUCOSEKINETICSBUILDER  

	%  $Revision$
 	%  was created 04-Dec-2017 11:27:40 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.
 	
	properties 		
 	end

	methods 
		  
 		function this = F18DeoxyGlucoseKineticsBuilder(varargin)
 			%% F18DEOXYGLUCOSEKINETICSBUILDER
 			%  Usage:  this = F18DeoxyGlucoseKineticsBuilder()

 			this = this@mlglucose.GlucoseKineticsBuilder(varargin{:});
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

