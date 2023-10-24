classdef F18DeoxyGlucoseModel < mlglucose.GlucoseModel
	%% F18DEOXYGLUCOSEMODEL  
    %  TO DO:  implement strategy design pattern for different solvers.

	%  $Revision$
 	%  was created 05-Dec-2017 21:20:07 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlraichle/src/+mlraichle.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.
 	
	properties    
        LC = 0.81 % Wu, et al., Molecular Imaging and Biology, 5(1), 32-41, 2003.
        fu = 1
        k1 = 3.946/60 % 1/s
        k2 = 0.3093/60
        k3 = 0.1862/60
        k4 = 0.01382/60
        t0 = 0 % s
        CBF = 55 % mL/hg/min
        CBV = 0.0383 % mL/hg
        
        sk1 = 1.254/60 % 1/s
        sk2 = 0.4505/60
        sk3 = 0.1093/60
        sk4 = 0.004525/60
        st0 = 30 % s
        
        fixed = false(5,1)
        fixedValue = nan(5,1)
        
        % Joanne Markham used the notation K_1 = V_B*k_{21}, rate from compartment 1 to 2.
        % Mean values from Powers xlsx "Final Normals WB PET PVC & ETS"
    end
    
    methods (Static)
        function name = parameterIndexToName(idx)
            names = {'k1' 'k2' 'k3' 'k4' 't0'};
            name = names{idx};
        end
        function Cwb  = plasma2wb(varargin)
            Cwb = mlraichle.RBCPartition.plasma2wb(varargin{:});
        end        
        function Cp   = wb2plasma(varargin)
            Cp = mlraichle.RBCPartition.wb2plasma(varargin{:});
        end
    end
    
	methods
        function a     = aifSpecificActivity(this)
            a = this.wb2plasma( ...
                this.aifSpecificActivity@mlglucose.GlucoseModel, ...
                this.blindedData_.hct, this.aifTimes);
        end
        function a     = aifSpecificActivityInterp(this)
            a = this.wb2plasma( ...
                this.aifSpecificActivityInterp@mlglucose.GlucoseModel, ...
                this.blindedData_.hct, this.aifTimeInterpolants);          
        end        
        function ps    = modelParameters(this)
            ps = [this.k1 this.k2 this.k3 this.k4 this.t0]';
        end
        function sps   = modelStdParameters(this)
            sps = [this.sk1 this.sk2 this.sk3 this.sk4 this.st0]';
        end
        function this  = constructSyntheticKernel(this)
            import mlglucose.*;
            this.kernel_ = F18DeoxyGlucoseKernel( ...
                this.CBV, ...
                this.aifTimesInterp, ...
                this.aifSpecificActivityInterp, ...
                this.scannerTimes, ...
                F18DeoxyGlucoseKernel.glucose( ...
                    asai, this.k1, this.k2, this.k3, this.k4, this.aifTimesInterp, this.CBV));
        end
        function this  = constructKernelWithData(this)
            this.kernel_ = mlglucose.F18DeoxyGlucoseKernel( ...
                this.CBV, ...
                this.aifTimesInterp, this.aifSpecificActivityInterp, this.scannerTimes, this.scannerSpecificActivity);
        end
        function those = thisParameters2those(this, those)
            those.LC  = this.LC;
            those.fu  = this.fu;
            those.k1  = this.k1;
            those.k2  = this.k2;
            those.k3  = this.k3;
            those.k4  = this.k4;
            those.t0  = this.t0;
            those.CBF = this.CBF;
            those.CBV = this.CBV;   
            
            those.sk1  = this.sk1;
            those.sk2  = this.sk2;
            those.sk3  = this.sk3;
            those.sk4  = this.sk4;
            those.st0  = this.st0;
            
            those.fixed = this.fixed;
            those.fixedValue = this.fixedValue;
        end
        function this  = thoseParameters2this(this, those)
            this.LC  = those.LC;
            this.fu  = those.fu;
            this.k1  = those.k1;
            this.k2  = those.k2;
            this.k3  = those.k3;
            this.k4  = those.k4;
            this.t0  = those.t0;
            this.CBF = those.CBF;
            this.CBV = those.CBV;   
            
            this.sk1  = those.sk1;
            this.sk2  = those.sk2;
            this.sk3  = those.sk3;
            this.sk4  = those.sk4;
            this.st0  = those.st0;  
            
            this.fixed = those.fixed;
            this.fixedValue = those.fixedValue;       
        end
        
 		function this  = F18DeoxyGlucoseModel(varargin)
 			%% F18DEOXYGLUCOSEMODEL

 			this = this@mlglucose.GlucoseModel(varargin{:});
            
            this = this.setupKernel;                       
            this = this.setupFilesystem;
            this = this.checkModel;
 		end
    end 
    
    %% PROTECTED
    
    methods (Access = protected)        
        function ps   = mcmcParameters(this)
            %% MCMCPARAMETERS must be in heap memory for speed
            %  @return struct containing:
            %  fixed      is logical, length := length(this.modelParameters)
            %  fixedValue is numeric, "
            %  min        is numeric, "
            %  mean       is numeric, "
            %  max        is numeric, "
            %  std        is numeric, "
            %  nAnneal    =  20, number of loops per annealing temp
            %  nBeta      =  50, number of temperature steps
            %  nPop       =  50, number of population for annealing/burn-in and proposal/sampling
            %  nProposals = 100, number of proposals for importance sampling
            %  nSamples   is numeric, numel of independentData
            
            ps = struct( ...
                'fixed',      this.fixed, ...
                'fixedValue', this.fixedValue, ...
                'min',        [0.05/60 ...
                               this.positive(0.04517/60   - this.M*this.sk2), ...
                               this.positive(0.05827/60   - this.M*this.sk3), ...
                               this.positive(0.0040048/60 - this.M*this.sk4), ...
                                                           -this.M*this.st0]', ...
                'mean',       this.modelParameters, ...
                'max',        [20/60 ...
                               1.7332/60   + this.M*this.sk2 ...
                               0.41084/60  + this.M*this.sk3 ...
                               0.017819/60 + this.M*this.sk4 ...
                                             this.M*this.st0]', ...
                'std',        this.modelStdParameters, ...
                'nAnneal',    this.nAnneal, ...
                'nBeta',      this.nBeta, ...
                'nPop',       50, ...
                'nProposals', this.nProposals, ...
                'nSamples',   numel(this.independentData));
        end
        function this = mcmcParameters2model(this, solvr)
            assert(isa(solvr, 'mlbayesian.IMcmcSolver'));
            bfp = ensureRowVector(solvr.kernel.bestFitParams);
            sp  = ensureRowVector(solvr.kernel.stdParams);
            this.k1  = bfp(1);
            this.k2  = bfp(2);
            this.k3  = bfp(3);
            this.k4  = bfp(4);
            this.t0  = bfp(5);
            this.sk1 = sp(1);
            this.sk2 = sp(2);
            this.sk3 = sp(3);
            this.sk4 = sp(4);
            this.st0 = sp(5);
        end
    end
        
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

