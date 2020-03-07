classdef MetabDirector < handle & mlkinetics.MetabDirector
	%% METABDIRECTOR is the director of a builder design pattern.

	%  $Revision$
 	%  was created 16-Dec-2018 23:43:11 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties
 		
    end
    
    methods (Static)
        function this = constructHuang(varargin)
            %% is the client to the builder design pattern.  See also GoF, Design Patterns, Builder for notation.
            
            ip = inputParser;
            addParameter(ip, 'sessionData', [], @(x) isa(x, 'mlpipeline.ISessionData'));
            parse(ip, varargin{:});
            sessd = ip.Results.sessionData;
            
            this = mlglucose.MetabDirector( ...
                'calBuilder', mlpet.CalibrationBuilder('sessionData', sessd), ...
                'aifBuilder', mlcapintec.CapracBuilder('sessionData', sessd), ...
                'roiBuilder', mlrois.RoisBuilder('sessionData', sessd), ...
                'tacBuilder', mlsiemens.BiographMMRBuilder('sessionData', sessd), ...
                'modelBuilder', mlglucose.HuangBuilder('sessionData', sessd));
            this.construct();
            this.writeResults();
        end
        function this = constructDagogoJack(varargin)
            this = mlglucose.MetabDirector( ...
                mlglucose.DagogoJackBuilder(varargin{:}));
        end
        function this = constructFokkerPlanck(varargin)
            this = mlglucose.MetabDirector( ...
                mlglucose.FokkerPlanckBuilder(varargin{:}));
        end
        function this = constructKPZ(varargin)
            this = mlglucose.MetabDirector( ...
                mlglucose.KPZBuilder(varargin{:}));
        end
    end

	methods
        function buildCalibration(this)
            this.calBuilder.readMeasurements();
            this.calBuilder.selectCalHierarchy();
            this.calBuilder.propagateEfficiencies();
        end
        function buildAif(this)
            this.aifBuilder.calibration = this.calBuilder.wellCounterCal;
            this.aifBuilder.readMeasurements();
            this.aifBuilder.propagateEfficiencies();
        end
        function buildRoi(this)
            this.roiBuilder.readMeasurements();
        end
        function buildTac(this)
            this.tacBuilder.calibration = this.calBuilder.scannerCal;
            this.tacBuilder.readMeasurements();
            this.tacBuilder.sampleMeasurements( ...
                this.roiBuilder);
            this.tacBuilder.propagateEfficiencies();
        end
        function buildLab(this)
            this.labBuilder.readMeasurements();
        end
        function buildModel(this)
            this.modelBuilder.sampleModel( ...
                this.aifBuilder, this.tacBuilder, this.labBuilder);
            this.modelBuilder.estimateParameters();
        end
        
 		function this = MetabDirector(varargin)
 			%% METABDIRECTOR is the director of a builder design pattern.
 			%  @param required builder is mlglucose.MetabBuilder.

 			this = this@mlkinetics.MetabDirector(varargin{:});
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

