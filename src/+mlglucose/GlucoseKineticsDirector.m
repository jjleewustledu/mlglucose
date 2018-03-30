classdef GlucoseKineticsDirector < mlkinetics.AbstractIterableKineticsDirector
	%% GLUCOSEKINETICSDIRECTOR  

	%  $Revision$
 	%  was created 04-Dec-2017 12:15:52 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee. 	
    
    properties (Dependent)
        oxygenDirector
        product
    end
    
	methods 
        
        %% GET/SET
        
        function g = get.oxygenDirector(this)
            g = this.oxygenDirector_;
        end
        function g = get.product(this)
            g = { mlglucose.GlucoseView('rates',          this.rates_, ...
                                        'physiologicals', this.physiologicals_) ...
                  this.oxygenDirector_.product };
        end
        
        function this = set.oxygenDirector(this, s)
            assert(isa(s, 'mloxygen.OxygenDirector'))
            this.oxygenDirector_ = s;
        end
        
        %%
        
        function this = constructRates(this)
            this.oxygenDirector_ = this.oxygenDirector_.constructRates;
            this.rates_ = this.kineticsBuilder_.buildRates;
        end
        function this = constructPhysiological(this)
            assert(~isempty(this.rates_));
            this.oxygenDirector_ = this.oxygenDirector_.constructPhysiological;
            this.physiologicals_ = this.kineticsBuilder_.buildPhysiologicals(this.rates_);
        end
        
 		function this = GlucoseKineticsDirector(varargin)
 			%% GLUCOSEKINETICSDIRECTOR
 			%  @param glucoseBldr is an mlkinetics.IKineticsBuilder.
            %  @param named model is an mlglucose.GlucoseModel.

            this = this@mlkinetics.AbstractIterableKineticsDirector(varargin{:});
 			ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'oxygenDir', mlkinetics.NullKineticsDirector(), @(x) isa(x, 'mlkinetics.IKineticsDirector'));
            addParameter(ip, 'glucoseBldr', @(x) isa(x, 'mlglucose.GlucoseKineticsBuilder'));
            addParameter(ip, 'roisBldr', @(x) isa(x, 'mlrois.IRoisBuilder'));
            parse(ip, varargin{:});
            
            this.oxygenDirector_ = ip.Results.oxygenDir;
            this.kineticsBuilder_ = ip.Results.glucoseBldr;
            this.roisBuilder_ = ip.Results.roisBldr;
 		end
    end 

    %% PRIVATE
    
    properties (Access = private)
        oxygenDirector_
 	end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end
