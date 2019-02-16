classdef GlucoseKineticsDirector < handle & mlkinetics.AbstractKineticsDirector
	%% GLUCOSEKINETICSDIRECTOR  

	%  $Revision$
 	%  was created 04-Dec-2017 12:15:52 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee. 	
    
    
    methods (Static)
        function this = createKinetics(bldr, varargin)
            this = mlglucose.GlucoseKineticsDirector('kineticsBldr', bldr, varargin{:});
        end
    end
    
	methods 
        
        %%
        
        function constructKs(this)
            this.kineticsBuilder_.buildKs;
        end
        function constructBrainmaskEstimate(this)
            if (~isa(this.roisBuilder_, 'mlrois.BrainmaskBuilder'))
                this.roisBuilder_ = mlrois.BrainmaskBuilder();
            end
            this.constructKs;
        end
        function constructAvisRegionalEstimates(this)
        end
        function constructAparcAsegEstimates(this)
        end
        function constructAparcA2009sAsegEstimates(this)
        end     
        function constructVoxelEstimates(this)
        end     
        function constructCtxglc(this)
        end
        function constructCmrglc(this)
        end   
        function constructFreeglc(this)
        end
        function constructEnet(this)
        end
        
    end
    
    methods (Access = protected)
                
 		function this = GlucoseKineticsDirector(varargin)
 			%% GLUCOSEKINETICSDIRECTOR

            this = this@mlkinetics.AbstractKineticsDirector(varargin{:});
 		end
    end 
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end
