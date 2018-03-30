classdef (Abstract) AbstractSolvedResults 
	%% ABSTRACTSOLVEDRESULTS  

	%  $Revision$
 	%  was created 12-Dec-2017 18:52:20 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlanalysis/src/+mlanalysis.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.
 	
    properties (Constant)
        DT_NYQUIST = 0.5
    end
    
	properties (Dependent)
 		useSynthetic
    end
    
    methods (Abstract)        
        this = constructSelfTest(this)
    end

	methods 
        
        %% GET/SET
        
        function g = get.useSynthetic(this)
            if (~isempty(his.oxygenDirector_))
                g = this.oxygenDirector_.useSynthetic;
            end            
            if (~isempty(his.glucoseDirector_))
                g = g || this.glucoseDirector_.useSynthetic;
            end
        end
        function this = set.useSynthetic(this, tf)
            assert(islogical(tf));          
            if (~isempty(this.oxygenDirector_))
                this.oxygenDirector_.useSynthetic = tf;
            end            
            if (~isempty(this.glucoseDirector_))
                this.glucoseDirector_.useSynthetic = tf;
            end            
        end
		  
        %%
            
 		function this = AbstractSolvedResults(varargin)
 		end
    end 
    
    %% PROTECTED
    
    properties (Access = protected)
        glucoseDirector_
        oxygenDirector_
    end
    
    

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

