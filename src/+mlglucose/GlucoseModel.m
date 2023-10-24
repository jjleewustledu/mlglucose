classdef (Abstract) GlucoseModel < mlkinetics.KineticsModel_20190307
	%% GLUCOSEMODEL  

	%  $Revision$
 	%  was created 05-Dec-2017 21:20:18 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.
 	
    properties      
    end
    
    methods (Abstract)
        those = thisParameters2those(this)
        this  = thoseParameters2this(this, those)
    end
    
	methods 
		  
        function r = glcRates(this)
            r = thisParameters2those(this, mlglucose.GlucoseRates);
        end
        function phys = glcRates2physiologicals(this, varargin)
            phys = mlglucose.GlucosePhysiologicals(this, varargin{:}); % plasma2wb, wb2plasma will be needed
        end 
        
 		function this = GlucoseModel(varargin)
 			%% GLUCOSEMODEL
 			%  @param named blindedData is an mlglucose.BlindedData
            
            this = this@mlkinetics.KineticsModel_20190307(varargin{:});            
            import mlpet.*;
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'blindedData', mlglucose.BlindedData(), @(x) isa(x, 'mlglucose.BlindedData'));
            parse(ip, varargin{:}); 
            this.blindedData_ = ip.Results.blindedData;
 		end
 	end 
    
    %% PROTECTED
    
    properties (Access = protected)
        blindedData_
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end
