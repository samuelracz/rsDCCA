% rtDCCA - Real-time detrended cross-correlation analysis of multiple 
%          time series components.
%          
%           Contains a custom defined class in which 'single-click' DCCC 
%           calculations can be performed on custom analysis scales. The
%           class also contains methods to simplify data addition and 
%           removal dynamically during runtime.
%
% Usage:
%   >> obj = rtDCCA             - define a rtDCCA class object
%   >> runDCCC(obj, method)     - DCCC Calcuation via matrix or pairwise
%                                 methods = {1,2}. Default value is {1}.
%
% Inputs:
%   obj.datapoints  - input time series ([channel,datapoints] format)
%   obj.scales      - time scales (array of 2^n recommended, where 
%                     3=< n <[data length/4]; in an 1 x [number of scales] 
%                     array format, ascending order)
%
% Outputs:
%   obj.FC          - calculated DCCC/DFA metrics timestamped with
%                     the point in dataflow it was calculated at 
%
% All inputs and outputs are built into the class object.
%
% For more information, see:
%   Z. Kaposzta, A. Czoch, O. Stylianou, K. Keumbi, P. Mukli, A. Eke, and F. Racz,
%   Real-Time Algorithm for Detrended Cross-Correlation Analysis of Long-Range Coupled Processes. 
%   Frontiers in Physiology. 13. 817268. (2022), doi: 10.3389/fphys.2022.817268. 
% See example use in: 'example.m'.
% 
% Copyright (C) 2022 Frigyes S. Racz and Zalan Kaposzta 
% Semmelweis University, Hungary 
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  US.

classdef rtDCCA < matlab.mixin.SetGet
    properties
        datapoints          % data buffer
        datastream = [0 0]; % current point of data analysis
        lables              % (optional) lables
        scales              % scales
        indices             % working variable
        FC                  % DCCC/DFA results
    end
    methods
%       Add datapoints with arbitrary channel count to the data buffer
        function addtoBuffer(obj,data)
            lr = size(obj.datapoints);
            ld = size(data);
            if (ld(1) == lr(1) || lr(1) == 0)
                obj.datapoints = [obj.datapoints data];
            else
                obj.datapoints(end+1:ld(1),1:end) = NaN(ld(1)-lr(1),lr(2));
                obj.datapoints = [obj.datapoints data];
            end
        end
%       Remove an amount of datapoints from all channels
        function emptyBuffer(obj,amount)
            if nargin < 2 ||  amount > size(obj.datapoints,2)
                amount = size(obj.datapoints,2);
            end
            obj.datapoints(:,1:amount) = [];
            obj.datastream(1) = obj.datastream(1) - amount;
            obj.datastream(2) = obj.datastream(2) + amount;
            obj.indices = obj.indices-amount;
        end
%       Execute real-time DCCC Calculation on data defined via 'datapoints'
%       on scales defined via 'scales' variables.
        function runDCCC(obj,method)
            if nargin < 2
                method = 1;
            end
            % Create start and end indices
            if (isempty(obj.indices))
                obj.indices = zeros(length(obj.scales),1);
                for i=1:length(obj.indices)
                obj.indices(i,1) = obj.scales(i)*4;
                end
            end            
            if (isempty(obj.datastream))
                obj.datastream = [0 0];
            end
            % Set scales and data for current analysis
            s_min = min(obj.scales);
            if length(obj.datapoints)>=obj.datastream(1)+s_min
                for p=obj.datastream(1):s_min:size(obj.datapoints,2)
                    obj.datastream(1) = p;
                    s = 0;
                    for o=1:length(obj.scales)
                        if obj.datastream(1) >= obj.indices(o)
                            s = obj.scales(1:o);
                            obj.indices(o) = obj.indices(o)+obj.scales(o);
                        else
                            break;
                        end
                    end
                    if s~=0
                        if method == 1
                            fc = omDCCC(obj,s);
                        else
                            fc = opDCCC(obj,s);
                        end
                        obj.FC = [obj.FC, fc];
                    end
                end
            else
                disp('Reached the current end of datapoints with analysis')
            end
        end
%       Fast parallel detrended cross-correlation coefficient calculation
        function [ Fc ] = omDCCC(obj,s)
            % Initializing main variables
            data = obj.datapoints(:,...     % data currently analysed
                1+obj.datastream(1)-max(s)*4:obj.datastream(1));            
            nch = size(data,1);             % number of parallel timeseries
            N = size(data,2);               % number of datapoints
            ns = length(s);                 % number of time scales
            smin = min(s);                  % smallest time scale
            smax = max(s);                  % largest time scale
            Fc = struct('t',[],'DCCA',[]);  % output structure
            
            % Initializing temporary variables
            DCCA_tmp = zeros(nch,nch,ns,4);
            denom_tmp = zeros(nch,nch,ns);
            for si = 1:ns
                denom_tmp(:,:,si) = repmat(1/(4*smax/s(si)), [nch,nch]);
            end
            
            % Initializing helper variables
            sx = zeros(nch,ns);
            sxi = zeros(nch,ns);
            sx2 = zeros(nch,nch,ns);
            
            i = zeros(nch,ns);
            wi = 0;
            Mi = 0;
            Fi = 0;
            
            % Start point of data analysis
            for t = 1:N
                if t == 1
                    Xi = data(:,t);
                    wi = 1;
                    Mi = 1;
                else
                    Xi = Xi + data(:,t);	% cumulative summation
                end
                
                i = i + 1;                	% relative index increase
                
                % Updating helper variables
                sx = sx + repmat(Xi,[1,ns]);
                sxi = sxi + repmat(Xi,[1,ns]).*i;
                sx2 = sx2 + repmat(Xi*Xi',[1,1,ns]);
                
                if i(1) < smin              % if smallest window is not
                    continue            	% yet filled, contiune with
                end                         % next datapoint
                                             
                for w = 1:ns
                    if i(1,w) < s(w)        % only iterate through
                        break               % windows that are full
                    else
                        % One-pass detrending
                        mx = -6*(-2*sxi(:,w)+(s(w)+1)*sx(:,w))./(s(w)*...
                            (s(w)^2-1));
                        bx = sx(:,w)./s(w) - mx.*(s(w)+1)./2;
                        
                        % One-pass calculation of detrended covariance                        
                        Fw =(((s(w)+1)*(2*s(w)+1))/6).*(mx*mx') + ...
                            ((s(w)+1)/2).*((mx*bx') + (mx*bx')') + ...
                            (bx*bx') + (1/s(w)).*(squeeze(sx2(:,:,w)) - ...
                            ((mx*sxi(:,w)') + (mx*sxi(:,w)')') - ...
                            ((bx*sx(:,w)') + (bx*sx(:,w)')'));                                                
                        DCCA_tmp(:,:,w,Mi) = DCCA_tmp(:,:,w,Mi) + (Fw);
                        
                        % Resetting helper variables
                        sx(:,w) = 0;
                        sxi(:,w) = 0;
                        sx2(:,:,w) = 0;
                        i(:,w) = 0;
                    end
                end
                
                % Increase/reset window offset, increase window index
                if wi < smax/smin
                    wi = wi + 1;
                    continue % if not all windows are filled, continue
                else
                    wi = 1;
                    i = zeros(nch,ns);
                end
                
                % If 4 windows are filled, calculate values
                switch Mi
                    case 4
                        Fi = Fi + 1;
                        Fc_tmp = denom_tmp.*sum(DCCA_tmp,4);
                        for j = 1:ns
                            DCCA_s = squeeze(Fc_tmp(:,:,j));
                            DFA = sqrt(diag(DCCA_s))*sqrt(diag(DCCA_s)');
                            DFA(eye(nch)==1) = sqrt(diag(DFA));
                            Fc_tmp(:,:,j) = DCCA_s./DFA;
                        end
                        Fc(Fi).t = obj.datastream(1) + obj.datastream(2);
                        Fc(Fi).DCCA = Fc_tmp;
                        DCCA_tmp = DCCA_tmp(:,:,:,2:end);
                        DCCA_tmp(:,:,:,end+1) = zeros(nch,nch,ns,1);
                    otherwise
                        Mi = Mi + 1;
                end
            end
        end
%       Fast linear detrended cross-correlation coefficient calculation
        function [ Fc ] = opDCCC(obj,s)
            
            data = obj.datapoints(:,...
                1+obj.datastream(1)-max(s)*4:obj.datastream(1));       
            Fc = struct('t',[],'DCCA',[]);
            
            for x = 1:size(data,1)
                for y = x:size(data,1)
                    % Initializing main variables
                    ts1 = data(x,:);
                    ts2 = data(y,:);
                    ns = length(s);
                    smin = min(s);
                    smax = max(s);
                    
                    % Initializing temporary variables
                    DCCA_tmp = zeros(4,ns);
                    denom_tmp = 4*smax./s;
                    
                    % Initializing temporary variables
                    sx = zeros(1,ns);
                    sxi = zeros(1,ns);                    
                    sy = zeros(1,ns);
                    syi = zeros(1,ns);                    
                    sxy = zeros(1,ns);
                    
                    i = zeros(1,ns);
                    wi = 0;
                    Mi = 0;
                    Fi = 0;
                    
                    for t = 1:length(ts1)
                        if t == 1
                            Xi = ts1(t);
                            Yi = ts2(t);
                            wi = 1;
                            Mi = 1;
                        else
                            Xi = Xi + ts1(t); % cumulative summation
                            Yi = Yi + ts2(t); % cumulative summation
                        end
                        
                        i = i + 1;            % increase relative index
                        
                        % Updating helper variables 
                        sx = sx + Xi;
                        sxi = sxi + Xi.*i;                        
                        sy = sy + Yi;
                        syi = syi + Yi.*i;                        
                        sxy = sxy + Xi*Yi;
                        
                        if i(1) < smin        % if smallest window is not 
                            continue          % filled, contiune with next 
                        end                   % datapoint
                        
                        for w = 1:ns
                            if i(w) < s(w)    % only iterate through 
                                break         % windows that are filled
                            else
                                % One-pass detrending
                                mx = -6*(-2*sxi(w)+(s(w)+1)*sx(w))/...
                                    (s(w)*(s(w)^2-1));
                                bx = sx(w)/s(w) - mx*(s(w)+1)/2;
                                
                                my = -6*(-2*syi(w)+(s(w)+1)*sy(w))/...
                                    (s(w)*(s(w)^2-1));
                                by = sy(w)/s(w) - my*(s(w)+1)/2;
                                
                                % One-pass calculation of detrended covar.
                                Fw = (mx*my*s(w)^2)/3 + (mx*my*s(w))/2 +...
                                    (mx*my)/6 + (mx*by*s(w))/2 +...
                                    (mx*by)/2 + (my*bx*s(w))/2 +...
                                    (my*bx)/2 + bx*by + (1/s(w))*sxy(w)...
                                    - (my/s(w))*sxi(w) - (by/s(w))*sx(w)...
                                    - (mx/s(w))*syi(w) - (bx/s(w))*sy(w);                                
                                DCCA_tmp(Mi,w) = DCCA_tmp(Mi,w) + (Fw);
                                
                                % Resetting helper variables
                                sx(w) = 0;
                                sxi(w) = 0;
                                sy(w) = 0;
                                syi(w) = 0;
                                sxy(w) = 0;
                                i(w) = 0;
                            end
                        end
                        
                        % Increase/reset window offset, increase window
                        if wi < smax/smin
                            wi = wi + 1;
                            continue    % if not all windows are filled
                        else
                            wi = 1;
                            i = zeros(1,ns);
                        end
                        
                        % If 4 windows are filled, calculate values
                        if Mi == 4
                          Fi = Fi + 1;
                          Fc(Fi).t = obj.datastream(1) + obj.datastream(2);
                          Fc_temp = sqrt((1./denom_tmp).*sum(DCCA_tmp,1));
                          Fc(Fi).DCCA(x,y,1:length(s)) = Fc_temp';
                          Fc(Fi).DCCA(y,x,1:length(s)) = Fc_temp';
                          DCCA_tmp = [DCCA_tmp(2:end,:); zeros(1,ns)];
                        end
                        
                        if Mi < 4
                            Mi = Mi + 1;
                        end
                    end
                end
            end
            % Calculate DCCC if possible (non-negative DCCA values)
            for k=1:size(Fc,2)
                for c=1:size(Fc(k).DCCA,3)
                    for x=1:size(Fc(k).DCCA,2)
                        for y=1:size(Fc(k).DCCA,2)
                            if (x ~= y)
                                ts12_dcca = Fc(k).DCCA(x,y,c);
                                ts1_dfa = Fc(k).DCCA(x,x,c);
                                ts2_dfa = Fc(k).DCCA(y,y,c);
                                Fc(k).DCCA(x,y,c) = ...
                                   (ts12_dcca*ts12_dcca)/(ts1_dfa*ts2_dfa);
                            end
                        end
                    end
                end
            end
        end
%       Pairwise cross-correlation analysis between 2 time-series
        function [ Fc ] = DCCA(obj,s)
            
            data = obj.datapoints(:,...
                1+obj.datastream(1)-max(s)*4:obj.datastream(1));
            Fc = struct('t',[],'DCCA',[]);
            
            % Initializing main variables
            ts1 = data(x,:);
            ts2 = data(y,:);
            ns = length(s);
            smin = min(s);
            smax = max(s);
            
            % Initializing temporary variables
            DCCA_tmp = zeros(4,ns);
            denom_tmp = 4*smax./s;
            
            % Initializing temporary variables
            sx = zeros(1,ns);
            sxi = zeros(1,ns);
            sy = zeros(1,ns);
            syi = zeros(1,ns);
            sxy = zeros(1,ns);
            
            i = zeros(1,ns);
            wi = 0;
            Mi = 0;
            Fi = 0;
            
            for t = 1:length(ts1)
                if t == 1
                    Xi = ts1(t);
                    Yi = ts2(t);
                    wi = 1;
                    Mi = 1;
                else
                    Xi = Xi + ts1(t); % cumulative summation
                    Yi = Yi + ts2(t); % cumulative summation
                end
                
                i = i + 1;            % increase relative index
                
                % Updating helper variables
                sx = sx + Xi;
                sxi = sxi + Xi.*i;
                sy = sy + Yi;
                syi = syi + Yi.*i;
                sxy = sxy + Xi*Yi;
                
                if i(1) < smin        % if smallest window is not
                    continue          % filled, contiune with next
                end                   % datapoint
                
                for w = 1:ns
                    if i(w) < s(w)    % only iterate through
                        break         % windows that are filled
                    else
                        % One-pass detrending
                        mx = -6*(-2*sxi(w)+(s(w)+1)*sx(w))/...
                            (s(w)*(s(w)^2-1));
                        bx = sx(w)/s(w) - mx*(s(w)+1)/2;
                        
                        my = -6*(-2*syi(w)+(s(w)+1)*sy(w))/...
                            (s(w)*(s(w)^2-1));
                        by = sy(w)/s(w) - my*(s(w)+1)/2;
                        
                        % One-pass calculation of detrended covar.
                        Fw = (mx*my*s(w)^2)/3 + (mx*my*s(w))/2 +...
                            (mx*my)/6 + (mx*by*s(w))/2 +...
                            (mx*by)/2 + (my*bx*s(w))/2 +...
                            (my*bx)/2 + bx*by + (1/s(w))*sxy(w)...
                            - (my/s(w))*sxi(w) - (by/s(w))*sx(w)...
                            - (mx/s(w))*syi(w) - (bx/s(w))*sy(w);
                        DCCA_tmp(Mi,w) = DCCA_tmp(Mi,w) + (Fw);
                        
                        % Resetting helper variables
                        sx(w) = 0;
                        sxi(w) = 0;
                        sy(w) = 0;
                        syi(w) = 0;
                        sxy(w) = 0;
                        i(w) = 0;
                    end
                end
                
                % Increase/reset window offset, increase window
                if wi < smax/smin
                    wi = wi + 1;
                    continue    % if not all windows are filled
                else
                    wi = 1;
                    i = zeros(1,ns);
                end
                
                % If 4 windows are filled, calculate values
                if Mi == 4
                    Fi = Fi + 1;
                    Fc(Fi).t = obj.datastream(1) + obj.datastream(2);
                    Fc(Fi).DCCA = (1./denom_tmp).*sum(DCCA_tmp,1);
                    DCCA_tmp = [DCCA_tmp(2:end,:); zeros(1,ns)];
                end
                
                if Mi < 4
                    Mi = Mi + 1;
                end
            end
        end
    end
end







