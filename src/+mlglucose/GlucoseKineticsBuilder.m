classdef (Abstract) GlucoseKineticsBuilder < mlkinetics.AbstractKineticsBuilder
	%% GLUCOSEKINETICSBUILDER specifies an abstract interface for creating parts of an IGlucoseKinetics object.

	%  $Revision$
 	%  was created 04-Dec-2017 11:28:24 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.
 	
    properties
    end

	methods     
        
        %%
		
        function rates = buildRates(this)
            this.solver_ = this.solver_.estimateParameters;
            assert(this.solver_.isfinished);
            this.plot;
            this.saveFigures;
            this.save;
            this.writetable;
            rates = this.solver_.model.glcRates;
        end
        function phys  = buildPhysiologicals(this, rates)
            assert(~isempty(rates));
            this.save;
            this.writetable;
            phys = this.solver_.model.glcRates2physiologicals(rates);
        end
        
 		function this = GlucoseKineticsBuilder(varargin)
 			%% GLUCOSEKINETICSBUILDER
            
            this = this@mlkinetics.AbstractKineticsBuilder(varargin{:});
 		end
    end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

