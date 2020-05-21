classdef Huang1980Model
	%% HUANG1980MODEL provides model data and methods to the strategy design pattern comprising
    %  mlglucose.{Huang1980, Huang1980Nest, Huang1980SimulAnneal, Huang1980HMC, Huang1980 LM, Huang1980BFGS}.
    %  It operates on single voxels.

	%  $Revision$
 	%  was created 10-Apr-2020 17:55:02 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlglucose/src/+mlglucose.
 	%% It was developed on Matlab 9.7.0.1319299 (R2019b) Update 5 for MACI64.  Copyright 2020 John Joowon Lee.
 	
	properties
        artery_interpolated
        glc
        hct
        LC
 		map
        times_sampled
        v1
    end
    
    methods (Static)
        function qs       = solution(ks, v1, artery_interpolated)
            %  @param artery_interpolated is uniformly sampled with at high sampling freq. starting at time = 0.
            
            k1 = ks(1);
            k2 = ks(2);
            k3 = ks(3);
            k4 = ks(4);
            scale = 1;
            
            n = length(artery_interpolated);
            times = 0:1:n-1;
            k234 = k2 + k3 + k4;         
            bminusa = sqrt(k234^2 - 4 * k2 * k4);
            alpha = 0.5 * (k234 - bminusa);
            beta  = 0.5 * (k234 + bminusa);   
            conva = conv(exp(-alpha .* times), artery_interpolated);
            convb = conv(exp(-beta  .* times), artery_interpolated);
            conva = conva(1:n);
            convb = convb(1:n);
            conv2 = (k4 - alpha) .* conva + (beta - k4) .* convb;
            conv3 =                 conva -                convb;
            q2 = (k1 / bminusa)      * conv2;
            q3 = (k3 * k1 / bminusa) * conv3;
            qs = v1 * (artery_interpolated + scale * (q2 + q3));            
        end
        function qs       = sampled(ks, v1, artery_interpolated, times_sampled)
            %  @param artery_interpolated is uniformly sampled with at high sampling freq.
            %  @param times_sampled are samples scheduled by the time-resolved PET reconstruction
            
            import mlglucose.Huang1980Model.solution  
            qs = solution(ks, v1, artery_interpolated);
            qs = qs(round(times_sampled - times_sampled(1) + 1));
        end
        function [dqs,qs] = grad_solution(ks, v1, artery_interpolated)
            k1 = ks(1);
            k2 = ks(2);
            k3 = ks(3);
            k4 = ks(4);
            scale = 1;
            
            n = length(artery_interpolated);
            times = 0:1:n-1;
            k234 = k2 + k3 + k4;
            bminusa = sqrt(k234^2 - 4 * k2 * k4);
            alpha = 0.5 * (k234 - bminusa);
            beta  = 0.5 * (k234 + bminusa);
            conva  = conv(exp(-alpha .* times),          artery_interpolated);
            convb  = conv(exp(-beta  .* times),          artery_interpolated);
            convta = conv(exp(-alpha .* times) .* times, artery_interpolated);
            convtb = conv(exp(-beta  .* times) .* times, artery_interpolated);
            conva  = conva( 1:n);
            convb  = convb( 1:n);  
            convta = convta(1:n);
            convtb = convtb(1:n);            
            conv2  = (k4 - alpha) * conva + (beta - k4) * convb;
            conv3  =                conva -               convb;
            
            q2 = (k1 / bminusa)      * conv2;
            q3 = (k3 * k1 / bminusa) * conv3;
            qs = v1 * (artery_interpolated + scale * (q2 + q3)); 

            part_bminusa_k2 = (k234 - 2 * k4) / bminusa;
            part_a_k2 = 0.5 * (1 - part_bminusa_k2);
            part_b_k2 = 0.5 * (1 + part_bminusa_k2);

            part_bminusa_k3 =  k234 / bminusa;
            part_a_k3 = 0.5 * (1 - part_bminusa_k3);
            part_b_k3 = 0.5 * (1 + part_bminusa_k3);

            part_bminusa_k4 = (k234 - 2 * k2) / bminusa;
            part_a_k4 = 0.5 * (1 - part_bminusa_k4);
            part_b_k4 = 0.5 * (1 + part_bminusa_k4);
            
            part_q2_k1 = (1/bminusa) * conv2;
            
            part_q2_k2 = ...
                -(k1 * part_bminusa_k2 / bminusa^2) * conv2 + ...
                 (k1 / bminusa) * ( ...
                     part_a_k2 * (-conva - (k4 - alpha) * convta) + ...
                     part_b_k2 * ( convb - (beta - k4)  * convtb) ...
                 );
            
            part_q2_k3 = ...
                -(k1 * part_bminusa_k3 / bminusa^2) * conv2 + ...
                 (k1 / bminusa) * ( ...
                     part_a_k3 * (-conva - (k4 - alpha) * convta) + ...
                     part_b_k3 * ( convb - (beta - k4)  * convtb) ...
                 );
            
            part_q2_k4 = ...
                -(k1 * part_bminusa_k4 / bminusa^2) * conv2 + ...
                 (k1 / bminusa) * ( ...
                     conva - ...
                     part_a_k4 * conva -...
                     part_a_k4 * (k4 - alpha) * convta + ...
                     part_b_k4 * convb - ...
                     convb - ...
                     part_b_k4 * (beta - k4)  * convtb ...
                 );

            part_q3_k1 = (k3 / bminusa) * conv3;
            
            part_q3_k2 = ...
                -(k3 * k1 * part_bminusa_k2 / bminusa^2) * conv3 - ...
                 (k3 * k1 / bminusa) * ( ...
                     -part_a_k2 * convta + ...
                      part_b_k2 * convtb ...
                 );
            
            part_q3_k3 = ...
                (k1 / bminusa) * conv3 -  ...
                (k3 * k1 * part_bminusa_k3 / bminusa^2) * conv3 + ...
                (k3 * k1 / bminusa) * ( ...
                    -part_a_k3 * convta + ...
                     part_b_k3 * convtb ...
                );
            
            part_q3_k4 = ...
                -(k3 * k1 * part_bminusa_k4 / bminusa^2) * conv3 + ...
                 (k3 * k1 / bminusa) * ( ...
                     -part_a_k4 * convta + ...
                      part_b_k4 * convtb ...
                 );
            
            part_qs_k1 = part_q2_k1 + part_q3_k1; 
            part_qs_k2 = part_q2_k2 + part_q3_k2; 
            part_qs_k3 = part_q2_k3 + part_q3_k3; 
            part_qs_k4 = part_q2_k4 + part_q3_k4;
            dqs = v1 * scale * [part_qs_k1; part_qs_k2; part_qs_k3; part_qs_k4];
        end
        function [dqs,qs] = grad_sampled(ks, v1, artery_interpolated, times_sampled)
            import mlglucose.Huang1980Model.grad_solution  
            [dqs,qs] = grad_solution(ks, v1, artery_interpolated);
            dqs = dqs(:, ceil(times_sampled)); % kBq*s/mL
            qs  =  qs(:, ceil(times_sampled)); % kBq/mL
        end        
        function dks      = grad_ks(dqs, Z, Sigma)
            dks = dqs*Z./Sigma; % 4 x 1
        end
         
        function logp = log_likelihood(Z, Sigma)
            %% for Huang1980HMC
            
            logp = sum(-log(Sigma) - 0.5*log(2*pi) - 0.5*Z.^2); % scalar
        end
        function loss = simulanneal_objective(ks, v1, artery_interpolated, times_sampled, qs0, sigma0)
            import mlglucose.Huang1980Model.sampled          
            qs = sampled(ks, v1, artery_interpolated, times_sampled);            
            loss = 0.5 * sum((1 - qs ./ qs0).^2) / sigma0^2; % + sum(log(sigma0*qs0)); % sigma ~ sigma0 * qs0
        end   
        function m = preferredMap()
            m = containers.Map;
            m('k1')  = struct('min',  0.001, 'max', 0.5,  'init', 0.03,   'sigma', 0.02);
            m('k2')  = struct('min',  eps,   'max', 0.1,  'init', 0.003,  'sigma', 0.01);
            m('k3')  = struct('min',  eps,   'max', 0.05, 'init', 0.003,  'sigma', 0.002);
            m('k4')  = struct('min',  eps,   'max', 0.01, 'init', 0.0003, 'sigma', 0.001);
        end
        function rop  = rbc_over_plasma(t)
            %% RBCOVERPLASMA is [FDG(RBC)]/[FDG(plasma)]
            
            t   = t/60;      % sec -> min
            a0  = 0.814104;  % FINAL STATS param  a0 mean  0.814192	 std 0.004405
            a1  = 0.000680;  % FINAL STATS param  a1 mean  0.001042	 std 0.000636
            a2  = 0.103307;  % FINAL STATS param  a2 mean  0.157897	 std 0.110695
            tau = 50.052431; % FINAL STATS param tau mean  116.239401	 std 51.979195
            rop = a0 + a1*t + a2*(1 - exp(-t/tau));
        end
        function Cp   = wb2plasma(Cwb, hct, t)
            import mlglucose.Huang1980Model.rbc_over_plasma;
            Cp = Cwb./(1 + hct*(rbc_over_plasma(t) - 1));
        end
    end

	methods		  
 		function this = Huang1980Model(varargin)            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'map', this.preferredMap(), @(x) isa(x, 'containers.Map'))
            addParameter(ip, 'v1', 0.038, @isnumeric)
            addParameter(ip, 'times_sampled', [], @isnumeric)
            addParameter(ip, 'artery_interpolated', [], @isnumeric)
            addParameter(ip, 'glc', 4.7, @isnumeric)
            addParameter(ip, 'hct', 0.4, @isnumeric)
            addParameter(ip, 'LC', 0.81, @isnumeric)
            addParameter(ip, 'convert_wb2plasma', true, @islogical)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            this.map = ipr.map;
            this.v1 = ipr.v1;
            this.times_sampled = ipr.times_sampled;
            if this.times_sampled(end)+1 ~= length(ipr.artery_interpolated)
                this.artery_interpolated = ...
                    pchip(0:length(ipr.artery_interpolated)-1, ipr.artery_interpolated, 0:this.times_sampled(end));
            else
                this.artery_interpolated = ipr.artery_interpolated;
            end
            this.glc = ipr.glc;
            this.hct = ipr.hct;            
            if (this.hct > 1)
                this.hct = this.hct/100;
            end
            this.LC = ipr.LC;
            
            if ipr.convert_wb2plasma
                import mlglucose.Huang1980Model.wb2plasma
                t = 0:this.times_sampled(end);
                this.artery_interpolated = wb2plasma(this.artery_interpolated, this.hct, t);
            end
        end
        
        function fdg  = solution_simulated(this, varargin)
            ip = inputParser;
            addRequired(ip, 'ks', @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            fdg = mlglucose.Huang1980Model.sampled(ipr.ks, this.v1, this.artery_interpolated, this.times_sampled);
        end
    end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

