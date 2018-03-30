classdef C11GlucoseModel < mlglucose.GlucoseModel
	%% C11GLUCOSEMODEL  

	%  $Revision$
 	%  was created 12-Dec-2017 16:28:38 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.
 	
    
	properties (Constant)
 		LC = 1
        LAMBDA_BLOOD = 0.7 % Dillon RS. Importance of the hematocrit in interpretation of blood sugar. 
                           % Diabetes 1965;14:672-674
    end
    
	properties
        fu
        k12
        k21
        k32
        k43
        t0 
        CBF
        CBV
    end
    
    methods (Static)        
        function Cwb = plasma2wb(Cp, hct, ~)
            if (hct > 1)
                hct = hct/100;
            end
            import mlglucose.*; 
            Cwb = Cp.*(1 + hct*(C11GlucoseModel.LAMBDA_BLOOD - 1));
        end
        function Cp  = wb2plasma(Cwb, hct, ~)
            if (hct > 1)
                hct = hct/100;
            end
            import mlglucose.*; 
            Cp = Cwb./(1 + hct*(C11GlucoseModel.LAMBDA_BLOOD - 1));
        end
    end

	methods     
        function prms = solverParameters(this)
            import mlglucose.*;
            switch (class(this.solver_))
                case 'mlanalysis.LevenbergMarquardt'
                    prms = C11GlucoseLMParameters;
                case 'mlbayesian.BretthorstMcmc'
                    prms = C11GlucoseMcmcParameters;
                case 'mlnest.NestedSamplingMain'
                    prms = C11GlucoseNestParameters;
                case 'mlstan.FlatHMC'
                    prms = C11GlucoseFlatHMCParameters;
                case 'mlstan.HierarchicalHMC'
                    prms = C11GlucoseHierarchicalHMCParameters;
                otherwise
                    error('mlraichle:unsupportedSwitchStrategy', ...
                        'C11GlucoseModel.solverParameters:class(solver)->%s', this.solver_);
            end
        end 
		  
 		function this = C11GlucoseModel(varargin)
 			%% C11GLUCOSEMODEL

 			this = this@mlglucose.GlucoseModel(varargin{:}); 	
            this = thoseParameters2this(this, this.solverParameters);		
 		end
 	end 
    
    %% PRIVATE
    
    methods (Access = private)
        function those = thisParameters2those(this, those)
            those.LC  = this.LC;
            those.fu  = this.fu;
            those.k12 = this.k12;
            those.k21 = this.k21;
            those.k32 = this.k32;
            those.k43 = this.k43;
            those.t0  = this.t0;
            those.CBF = this.CBF;
            those.CBV = this.CBV;            
        end
        function this = thoseParameters2this(this, those)
            this.LC  = those.LC;
            this.fu  = those.fu;
            this.k12 = those.k12;
            this.k21 = those.k21;
            this.k32 = those.k32;
            this.k43 = those.k43;
            this.t0  = those.t0;
            this.CBF = those.CBF;
            this.CBV = those.CBV;            
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

