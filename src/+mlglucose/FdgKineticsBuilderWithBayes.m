classdef FdgKineticsBuilderWithBayes < handle & mlglucose.AbstractFdgKineticsBuilder
	%% FDGKINETICSBUILDERWITHBAYES  

	%  $Revision$
 	%  was created 19-Aug-2018 01:31:46 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties (Dependent)
        fit
 		paramEstimate
 	end

	methods 
        
        %% GET
        
        function g = get.fit(this)
            g = this.fit_;
        end
        function g = get.paramEstimate(this)
            g = this.product;
        end
		  
        %%
        
        function buildKs(this)            
        end
        function buildCtxglc(this)
        end
        function buildCmrglc(this)
        end   
        function buildFreeglc(this)
        end
        function buildEnet(this)
        end
        function plot(this)
        end
        function print(this)
            print(this.fit_);
        end
        
 		function this = FdgKineticsBuilderWithBayes(varargin)
 			%% FDGKINETICSBUILDERWITHBAYES
 			%  @param .

 			this = this@mlglucose.AbstractFdgKineticsBuilder(varargin{:});
 		end
    end 
    
    %% PROTECTED
    
    properties (Access = protected)
        fit_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

