classdef GlucoseRates 
	%% GLUCOSERATES  

	%  $Revision$
 	%  was created 04-Dec-2017 14:59:20 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.
 	
	properties              
        LC
        fu
        k1
        k2
        k3
        k4
        t0
        CBF
        CBV
    end
    
	properties (Dependent)
 		k04
        k12
        k21
        k32
        k43
    end
    
    methods
        
        %% GET/SET
        
        function g = get.k04(this)
            g = this.CBF/this.CBV;
        end
        function g = get.k12(this)
            g = this.k2;
        end
        function g = get.k21(this)
            g = this.k1;
        end
        function g = get.k32(this)
            g = this.k3;
        end
        function g = get.k43(this)
            g = this.k4;
        end
        
        function this = set.k12(this, s)
            this.k2 = s;
        end
        function this = set.k21(this, s)
            this.k1 = s;
        end
        function this = set.k32(this, s)
            this.k3 = s;
        end
        function this = set.k43(this, s)
            this.k4 = s;
        end        
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

