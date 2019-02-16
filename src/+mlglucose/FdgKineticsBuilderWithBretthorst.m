classdef FdgKineticsBuilderWithBretthorst < handle & mlglucose.FdgKineticsBuilderWithBayes
	%% FDGKINETICSBUILDERWITHBRETTHORST  

	%  $Revision$
 	%  was created 17-Aug-2018 20:51:13 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties
 		
 	end

	methods 
        
        function plot(this, varargin)
            this.mcmc_.plotAnnealing;
            this.mcmc_.plotParameterCovariances;
            this.mcmc_.plotLogProbabilityQC;
            this.mcmc_.histStdOfError;
        end
		  
 		function this = FdgKineticsBuilderWithBretthorst(varargin)
 			%% FDGKINETICSBUILDERWITHBRETTHORST
 			%  @param .

 			this = this@mlglucose.FdgKineticsBuilderWithBayes(varargin{:});
 		end
    end 
    
    %% PRIVATE
    
    properties (Access = private)
        mcmc_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

