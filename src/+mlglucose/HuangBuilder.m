classdef HuangBuilder < mlglucose.MetabBuilder
	%% HUANGBUILDER  

	%  $Revision$
 	%  was created 17-Dec-2018 00:10:37 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties
 		
 	end

	methods 
        
        function obj  = createAifData(this)
        end
        function obj  = createTacData(this)
        end
        function obj  = createCalData(this)
        end
        function obj  = createRoiData(this)
        end
        function obj  = createParameters(this)
        end
		  
 		function this = HuangBuilder(varargin)
 			%% HUANGBUILDER
 			%  @param .

            this = this@mlglucose.MetabBuilder(varargin{:});

            this.aifData_ = ip.Results.referenceTissueData;
            this.tacData_ = ip.Results.tissueData;
            this.calData_ = ip.Results.calibrationData;
            this.solver_  = ip.Results.solver;
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

