classdef F18DeoxyGlucoseKernel < mlkinetics.KineticsKernel
	%% F18DEOXYGLUCOSEKERNEL is designed for recoding for compiled objects. 

	%  $Revision$
 	%  was created 22-Dec-2017 14:30:15 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.
 	
	properties
        CBV
 	end

    methods (Static)
        function this = main(varargin)
            this = mlglucose.F18DeoxyGlucoseKernel(varargin{:});
        end
        function ps   = adjustParams(ps)
            if (ps(4) > ps(3))
                tmp   = ps(3);
                ps(3) = ps(4);
                ps(4) = tmp;
            end
        end
        
        function alpha_ = a(k2, k3, k4)
            k234   = k2 + k3 + k4;
            alpha_ = k234 - sqrt(k234^2 - 4*k2*k4);
            alpha_ = alpha_/2;
        end
        function beta_  = b(k2, k3, k4)
            k234  = k2 + k3 + k4;
            beta_ = k234 + sqrt(k234^2 - 4*k2*k4);
            beta_ = beta_/2;
        end
        function q      = q2(Aa, k1, a, b, k4, t)
            scale = k1/(b - a);
            q = scale * conv((k4 - a)*exp(-a*t) + (b - k4)*exp(-b*t), Aa);
            q = q(1:length(t));
        end
        function q      = q3(Aa, k1, a, b, k3, t)
            scale = k3*k1/(b - a);
            q = scale * conv(exp(-a*t) - exp(-b*t), Aa);
            q = q(1:length(t));
        end
        function q      = glucose(Aa, k1, k2, k3, k4, t, cbv)
            import mlglucose.*;
            Aa = cbv*Aa;
            a  = F18DeoxyGlucoseKernel.a(k2, k3, k4);
            b  = F18DeoxyGlucoseKernel.b(k2, k3, k4);
            q  = F18DeoxyGlucoseKernel.q2(Aa, k1, a, b, k4, t) + ...
                 F18DeoxyGlucoseKernel.q3(Aa, k1, a, b, k3, t) + Aa;
        end
    end
    
    methods
        function edf = estimateData(this, ps)
            %% ESTIMATEDATAFAST needs small data structures that will all fit in L1 caches.
            
            qNyquist = this.glucose( ...
                this.aifSpecificActivityInterp, ps(1), ps(2), ps(3), ps(4), this.aifTimesInterp, this.CBV);
            edf = pchipDt(this.aifTimesInterp, qNyquist, this.independentData, ps(5));
        end
		  
 		function this = F18DeoxyGlucoseKernel(cbv, varargin)
 			%% F18DEOXYGLUCOSEKERNEL
            
            this = this@mlkinetics.KineticsKernel(varargin{:});
            assert(isnumeric(cbv));	
            this.CBV = cbv;
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

