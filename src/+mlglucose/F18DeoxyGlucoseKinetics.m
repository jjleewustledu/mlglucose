classdef F18DeoxyGlucoseKinetics < mlglucose.AbstractGlucoseKinetics
	%% F18DEOXYGLUCOSEKINETICS  

	%  $Revision$
 	%  was created 04-Dec-2017 11:15:43 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.
 	
	properties
 		
 	end

	methods 
		  
 		function this = F18DeoxyGlucoseKinetics(varargin)
 			%% F18DEOXYGLUCOSEKINETICS
 			%  Usage:  this = F18DeoxyGlucoseKinetics()

 			this = this@mlglucose.AbstractGlucoseKinetics(varargin{:});
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

