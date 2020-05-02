classdef Huang1980Nest < mlnest.Dynamics & mlglucose.Huang1980Strategy
	%% HUANG1980NEST  

	%  $Revision$
 	%  was created 10-Apr-2020 18:03:26 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.7.0.1319299 (R2019b) Update 5 for MACI64.  Copyright 2020 John Joowon Lee.
 	
	properties
        artery_interpolated
        main
        v1
 	end
    
	methods
 		function this = Huang1980Nest(varargin)            
            this = this@mlnest.Dynamics(varargin{:});
            
            this.artery_interpolated = this.model.artery_interpolated;
            this.v1 = this.model.v1;
        end
        
        function est  = Estimation(this, Obj)            
            Obn = Obj2native(this, Obj);
            %t0_ = Obn.t0;
            ks_ = [Obj.k1 Obn.k2 Obn.k3 Obn.k4]; 
            artery_ = this.artery_interpolated;
            times_smpl_ = this.times_sampled;   

            % adjust artery_ by t0_
            %artery_ = mlnest.Dynamics.slide_fast(artery_, t0_ - times_smpl_(1));

            % delegate to this.model
            est = this.model.sampled(ks_, this.v1, artery_, times_smpl_);
        end 
        function h = H(this, varargin)
            if isempty(this.main)
                h = [];
                return
            end
            h = this.main.H;
        end
        function [k,sk] = K1(this, varargin)
            [k,sk] = k1(this, varargin{:});
            k = k*this.v1;
            sk = sk*this.v1;
        end
        function [k,sk] = k1(this, varargin)
            [k,sk] = find_result(this, 'k1');
        end
        function [k,sk] = k2(this, varargin)
            [k,sk] = find_result(this, 'k2');
        end
        function [k,sk] = k3(this, varargin)
            [k,sk] = find_result(this, 'k3');
        end
        function [k,sk] = k4(this, varargin)
            [k,sk] = find_result(this, 'k4');
        end
        function lz = logZ(this, varargin)
            if isempty(this.main)
                lz = [];
                return
            end
            lz = this.main.logZ;
        end
        function this = solve(this, varargin)
            this.main = mlglucose.Huang1980Nest.run(this, varargin{:});
        end
        function [t,st] = t0(this, varargin)
            [t,st] = find_result(this, 't0');
        end
    end
    
    %% PROTECTED
    
    methods (Access = protected)
        function [m,sd] = find_result(this, lbl)
            r = this.main.results;
            [~,idx] = max(strcmp(r.flds, lbl));
            m = r.moment1(idx);
            sd = sqrt(r.moment2(idx) - m^2);
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

