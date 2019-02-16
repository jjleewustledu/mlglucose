classdef FdgKineticsBuilderWithStan < handle & mlglucose.FdgKineticsBuilderWithBayes
	%% FDGKINETICSBUILDERWITHSTAN  

	%  $Revision$
 	%  was created 17-Aug-2018 20:51:06 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties (Dependent)
        data
        model_code
        stan_filename
 	end

	methods        
        
        %% GET/SET
        
        function g = get.data(this)
            g = this.data_;
        end
        function g = get.model_code(this)
            g = this.model_code_;
        end
        function g = get.stan_filename(this)
            g = this.stan_filename_;
        end
        
        %%
        
        function buildKs(this)
            this.fit_ = stan('file', this.stan_filename, 'data', this.data, 'verbose', this.verbose);
            this.fit_.block(); % See also https://github.com/brian-lau/MatlabStan/wiki/Blocking-the-Matlab-command-line
            print(this.fit_);
            this.product_ = this.fit_.extract('permuted',true);
        end        
        function e = extract(this, varargin)
            %  @param named 'permuted' is boolean.
            %  @param named 'pars' names an element of this.model_code->parameters
            %  @returns individual chains per array element if 'permuted' == true.
            %  @returns struct with all parameters if 'permuted' == false.
            
            e = this.fit_.extract(varargin{:});
        end
        function peek(this)
            peek(this.fit_);
        end
        function traceplot(this)
            this.fit_.traceplot;
        end
		  
 		function this = FdgKineticsBuilderWithStan(varargin)
 			%% FDGKINETICSBUILDERWITHSTAN
 			%  @param .

 			this = this@mlglucose.FdgKineticsBuilderWithBayes(varargin{:});
            
            y = [151, 145, 147, 155, 135, 159, 141, 159, 177, 134, ...
                160, 143, 154, 171, 163, 160, 142, 156, 157, 152, 154, 139, 146, ...
                157, 132, 160, 169, 157, 137, 153, 199, 199, 214, 200, 188, 210, ...
                189, 201, 236, 182, 208, 188, 200, 221, 216, 207, 187, 203, 212, ...
                203, 205, 190, 191, 211, 185, 207, 216, 205, 180, 200, 246, 249, ...
                263, 237, 230, 252, 231, 248, 285, 220, 261, 220, 244, 270, 242, ...
                248, 234, 243, 259, 246, 253, 225, 229, 250, 237, 257, 261, 248, ...
                219, 244, 283, 293, 312, 272, 280, 298, 275, 297, 350, 260, 313, ...
                273, 289, 326, 281, 288, 280, 283, 307, 286, 298, 267, 272, 285, ...
                286, 303, 295, 289, 258, 286, 320, 354, 328, 297, 323, 331, 305, ...
                338, 376, 296, 352, 314, 325, 358, 312, 324, 316, 317, 336, 321, ...
                334, 302, 302, 323, 331, 345, 333, 316, 291, 324];
            y = reshape(y,30,5);
            x = [8 15 22 29 36];

            this.data_ = struct('N',size(y,1),'TT',size(y,2),'x',x,'y',y,'xbar',mean(x));
            this.model_code_ = {''};
            this.stan_filename_ = '/Users/jjlee/Tmp/rats.stan';
 		end
    end 
    
    %% PROTECTED
    
    properties (Access = protected)
        data_
        model_code_
        stan_filename_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

