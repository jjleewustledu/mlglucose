classdef (Abstract) Huang1980Strategy 
	%% HUANG1980STRATEGY is the interface for concrete strategies used with mlglucose.Huang1980.

	%  $Revision$
 	%  was created 10-Apr-2020 17:38:15 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.7.0.1319299 (R2019b) Update 5 for MACI64.  Copyright 2020 John Joowon Lee.
    
	properties (Abstract)
        artery_interpolated
        context
        fileprefix
        map                  % containers.Map containing model params as structs with fields:  min, max, init
        model                %
        Measurement          % external data
        results
        sigma0               % fraction of Measurement < 1
        times_sampled        % numeric times for Measurement; midpoints of frames for PET
        v1                   % blood volume fraction
        visualize 		
 	end

	methods (Abstract)
        K1(this)
        k1(this)
        k2(this)
        k3(this)
        k4(this)
        solve(this)
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

