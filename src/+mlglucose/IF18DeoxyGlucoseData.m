classdef (Abstract) IF18DeoxyGlucoseData 
	%% IF18DEOXYGLUCOSEDATA  

	%  $Revision$
 	%  was created 12-Feb-2018 03:11:07 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties (Abstract)  
        LC % = 0.81 % Wu, et al., Molecular Imaging and Biology, 5(1), 32-41, 2003.
        fu % = 1
        k1 % = 3.946/60 % 1/s
        k2 % = 0.3093/60
        k3 % = 0.1862/60
        k4 % = 0.01382/60
        t0 % = 0 % s
        CBF % = 55 % mL/hg/min
        CBV % = 0.0383 % mL/hg
        
        sk1 % = 1.254/60 % 1/s
        sk2 % = 0.4505/60
        sk3 % = 0.1093/60
        sk4 % = 0.004525/60
        st0 % = 30 % s
        
        fixed % = false(5,1)
        fixedValue % = nan(5,1)
        
        % Joanne Markham used the notation K_1 = V_B*k_{21}, rate from compartment 1 to 2.
        % Mean values from Powers xlsx "Final Normals WB PET PVC & ETS"
    end
    
    methods (Abstract)
        save(this, filename)
        plot(this, varargin)
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

