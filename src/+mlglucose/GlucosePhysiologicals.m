classdef GlucosePhysiologicals 
	%% GLUCOSEPHYSIOLOGICALS  

	%  $Revision$
 	%  was created 04-Dec-2017 14:59:29 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.
 	
	properties
    end
    
    methods
        
        %%
        
        function p = CMRglc(this)
        end
        function p = CTXglc(this)
        end
        function p = Enet(this)
        end
        function p = freeglc(this)
        end
        function p = waterMetab(this)
        end
        
        function this = GlucosePhysiologicals(varargin)
            %% GLUCOSEPHYSIOLOGICALS
            %  @param rates is an mlglucose.GlucoseRates
            
            ip = inputParser;
            addRequired(ip, 'rates', @(x) isa(x, 'mlglucose.GlucoseRates'));
            parse(ip, varargin{:});
            
            this.rates_ = ip.Results.rates;
        end
    end
    
    %% PRIVATE
    
    properties (Access = private)
        rates_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

