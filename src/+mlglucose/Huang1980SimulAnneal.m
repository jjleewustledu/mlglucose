classdef Huang1980SimulAnneal < mlpet.TCSimulAnneal
	%% HUANG1980SIMULANNEAL operates on single voxels/regions.

	%  $Revision$
 	%  was created 03-Jan-2020 17:45:56 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.7.0.1261785 (R2019b) Update 3 for MACI64.  Copyright 2020 John Joowon Lee.
 	
    properties
    end

	methods        
 		function this = Huang1980SimulAnneal(varargin)
 			%% HUANG1980SIMULANNEAL
            %  @param context is mlglucose.Huang1980.
            %  @param sigma0.
            %  @param fileprefix.
            
            this = this@mlpet.TCSimulAnneal(varargin{:});
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

