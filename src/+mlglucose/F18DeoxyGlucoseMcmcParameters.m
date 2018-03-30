classdef F18DeoxyGlucoseMcmcParameters < mlglucose.F18DeoxyGlucoseParameters 
	%% F18DEOXYGLUCOSEMCMCPARAMETERS  

	%  $Revision$
 	%  was created 12-Dec-2017 23:33:15 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.
 	
	properties (Dependent)
 		mapParams
 	end

	methods 
        
        %% GET/SET
        
        function this = set.mapParams(this, m)
            assert(isa(m, 'containers.Map'));
            this.mapParams_ = m;
        end
        function m    = get.mapParams(this)
            if (~isempty(this.mapParams_))
                m = this.mapParams_;
                return
            end            
            m = containers.Map;
            N = 20;
            
            % From Powers xlsx "Final Normals WB PET PVC & ETS"
            m('fu') = struct('fixed', 1, 'min', 0.01,                              'mean', this.fu, 'max',   1);  
            m('k1') = struct('fixed', 0, 'min', 0.05/60,                           'mean', this.k1, 'max',  20/60);
            m('k2') = struct('fixed', 0, 'min', max(0.04517/60   - N*this.sk2, 0), 'mean', this.k2, 'max',   1.7332/60   + N*this.sk2);
            m('k3') = struct('fixed', 0, 'min', max(0.05827/60   - N*this.sk3, 0), 'mean', this.k3, 'max',   0.41084/60  + N*this.sk3);
            m('k4') = struct('fixed', 0, 'min', max(0.0040048/60 - N*this.sk4, 0), 'mean', this.k4, 'max',   0.017819/60 + N*this.sk4);
            m('u0') = struct('fixed', 0, 'min', -100,                              'mean', this.u0, 'max', 100);  
            m('v1') = struct('fixed', 1, 'min', 0.01,                              'mean', this.v1, 'max',   0.1);  
        end
        
        %%
		  
 		function this = F18DeoxyGlucoseMcmcParameters(varargin)
 			%% F18DEOXYGLUCOSEMCMCPARAMETERS
 			%  Usage:  this = F18DeoxyGlucoseMcmcParameters()

 			
 		end
    end 
    
    %% PRIVATE
    
    properties (Access = private)
        mapParams_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

