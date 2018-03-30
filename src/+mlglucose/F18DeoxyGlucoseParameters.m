classdef F18DeoxyGlucoseParameters 
	%% F18DEOXYGLUCOSEPARAMETERS  

	%  $Revision$
 	%  was created 12-Dec-2017 23:30:17 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.
 	    
	properties
        LC = 0.81 % Wu, et al., Molecular Imaging and Biology, 5(1), 32-41, 2003.
        fu = 1 % fudge
        
        %% Mean values from Powers xlsx "Final Normals WB PET PVC & ETS"
        k1 = 3.946/60 % 1/s
        k2 = 0.3093/60
        k3 = 0.1862/60
        k4 = 0.01382/60
        t0 = 0 % s
        CBF = nan % mL/hg/min
        CBV = 0.0383 % mL/hg
        
        %% Std values from Powers
        sigmak1 = 1.254/60
        sigmak2 = 0.4505/60
        sigmak3 = 0.1093/60
        sigmak4 = 0.004525/60 
        
        %% plot
        xLabel = 'times/s'
        yLabel = 'activity/(Bq/mL)'
        notes  = ''
    end
    
    properties (Dependent)
        detailedTitle
    end

	methods 
        
        %% GET
        
        function dt   = get.detailedTitle(this)
            dt = sprintf('F18DeoxyGlucoseParameters\nLC %g, fu %g, k1 %g, k2 %g, k3 %g, k4 %g, t0 %g, CBV %g\n%S', ...
                this.LC, this.fu, this.k1, this.k2, this.k3, this.k4, this.t0, this.CBV, this.notes);
        end
        
        %%
		  
 		function this = F18DeoxyGlucoseParameters(varargin)
 			%% F18DEOXYGLUCOSEPARAMETERS
 			%  Usage:  this = F18DeoxyGlucoseParameters()

 			
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

