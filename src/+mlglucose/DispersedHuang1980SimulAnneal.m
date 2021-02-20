classdef DispersedHuang1980SimulAnneal < mlglucose.Huang1980SimulAnneal
	%% DISPERSEDHUANG1980SIMULANNEAL  

	%  $Revision$
 	%  was created 11-Aug-2020 20:44:37 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.7.0.1434023 (R2019b) Update 6 for MACI64.  Copyright 2020 John Joowon Lee.
    
    properties
        Dt
        registry
    end
    
	methods		  
 		function this = DispersedHuang1980SimulAnneal(varargin)
 			this = this@mlglucose.Huang1980SimulAnneal(varargin{:});
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            ip.PartialMatching = false;
            addParameter(ip, 'Dt', [], @isnumeric)
            parse(ip, varargin{:})
            this.Dt = ip.Results.Dt;
            this.registry = mlraichle.RaichleRegistry.instance();
        end  
        
        function fprintfModel(this)
            fprintf('Simulated Annealing:\n');
            for ky = 1:length(this.ks)
                fprintf('\tk%i = %g\n', ky, this.ks(ky));
            end
            fprintf('\tDt = %g\n', this.Dt);
            fprintf('\ttBuffer = %g\n', this.registry.tBuffer)   
            fprintf('\tsigma0 = %g\n', this.sigma0);
            fprintf('\tv1 = %g\n', this.v1);
            for ky = this.map.keys
                fprintf('\tmap(''%s'') => %s\n', ky{1}, struct2str(this.map(ky{1})));
            end
        end
        function s = sprintfModel(this)
            s = sprintf('Simulated Annealing:\n');
            for ky = 1:length(this.ks)
                s = [s sprintf('\tk%i = %g\n', ky, this.ks(ky))]; %#ok<AGROW>
            end
            s = [s sprintf('\tDt = %g\n', this.Dt)];
            s = [s sprintf('\ttBuffer = %g\n', this.registry.tBuffer)];   
            s = [s sprintf('\tsigma0 = %g\n', this.sigma0)];
            s = [s sprintf('\tv1 = %g\n', this.v1)];
            for ky = this.map.keys
                s = [s sprintf('\tmap(''%s'') => %s\n', ky{1}, struct2str(this.map(ky{1})))]; %#ok<AGROW>
            end
        end        
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

