classdef (Sealed) F18DeoxyGlucoseKit < handle & mlkinetics.TracerKit
    %% is a concrete factory ~ cyclotron
    %  
    %  Created 09-Jun-2022 12:30:26 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
    %  Developed on Matlab 9.12.0.1956245 (R2022a) Update 2 for MACI64.  Copyright 2022 John J. Lee.
    
    methods (Static)
        function this = instance(varargin)
            %% INSTANCE
            %  Args:
            %      imaging (mlfourd.IImaging): supporting NIfTI, 4dfp, ...
            %      taus (numeric): frame durations of dynamic scans.
            %      qualifier (text): 'initialize' resets singleton data.
            
            ip = inputParser;
            addRequired(ip, 'imaging', @(x) isa(x, 'mlfourd.IImaging'));
            addRequired(ip, 'taus', @isnumeric);
            addOptional(ip, 'qualifier', '', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            persistent uniqueInstance
            if (strcmp(ip.Results.qualifier, 'initialize'))
                uniqueInstance = [];
            end          
            if (isempty(uniqueInstance))
                this = mlglucose.F18DeoxyGlucoseKit(ipr.imaging, ipr.taus);
                uniqueInstance = this;
            else
                this = uniqueInstance;
            end
        end
    end 

    %% PRIVATE

    methods (Access = private)
        function this = F18DeoxyGlucoseKit(imaging, taus)
            re = regexp(imaging.fileprefix, '\S+_ses-(?<dtm>\d{14})_trc-(?<tra>\w+)_\S+', 'names');
            assert(contains(re.trc, 'fdg', 'IgnoreCase', true));
            dtm = datetime(re.dtm, 'InputFormat', 'yyyyMMddHHmmss');
            this.tracerData_ = mlkinetics.TracerData( ...
                'isotope', '18F', ...
                'tracer', 'FDG', ...
                'datetimeMeasured', dtm, ...
                'taus', taus);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
