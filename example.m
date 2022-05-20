% Example script - provides a simulated real-world use case for the code
%                  found packaged together with this example file.
% 
% For more information, see:
%   Z. Kaposzta, A. Czoch, O. Stylianou, K. Keumbi, P. Mukli, A. Eke, 
%   and F. Racz,
%   Real-Time Algorithm for Detrended Cross-Correlation Analysis of 
%   Long-Range Coupled Processes. 
%   Frontiers in Physiology. 13. 817268. (2022), 
%   doi: 10.3389/fphys.2022.817268. 
%
% See also: rtDCCA MATLAB class.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User defined variables for demonstrative purposes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chunksize = 100;    % defines how many simulated datapoints will be sent
memorylimit = 2000; % defines the limit of datapoints stored in memory
if ~exist('datachunk','var')
    datachunk = []; % creates a dummy variable which will be filled with
end                 % simulated datapoints

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example code of a theoretical real-world application
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('rt','var')
    rt = rtDCCA;                    % set up rtDCCA class
    rt.scales = [8,16,32,64,128];   % define analysis scales
end

addtoBuffer(rt,datachunk)           % add datapoints to the buffer
runDCCC(rt,1)                       % perform DCCC/DFA calculation
disp(['Analysis completed from [1 to ',int2str(rt.datastream(1)...
    +rt.datastream(2)),'].'])

if rt.datastream(1) > memorylimit
    emptyBuffer(rt,memorylimit/2)   % limit amount of memory allocated to
end                                 % the class' data buffer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [FOR DEMONSTRATIVE PURPOSES ONLY, NOT REQUIRED FOR APPLICATION]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creation and allocation of simulated datapoints
% This represents loading in data via an outside source, such as real-time
% datastream or a pre-existing recording. This code snippet is unnecessary 
% for real-world applications.
if ~exist('input_struct','var')
    addpath([cd,'\datageneration'])
    input_struct = struct('N',10000,'alpha',0.2,'delta',0.2,'gamma',1,...
                         'beta',1,'d1',0.4,'d4',0.4,'d2',0.3,'d3',0.3,...
                         'rho',0.9);
    [ data ] = func_mcARFIMA( input_struct );
    data = data';
end
prompt=['Press [Enter] to simulate loading a data chunk into the code '...
        'or press [0] to exit the simulation \n'];
x = input(prompt);
if isempty(x)
    datachunk = data(:,1:chunksize);
    data(:,1:chunksize) = [];
    if size(data,2) <= chunksize
        [ data ] = func_mcARFIMA( input_struct );
        data = data';
    end
    run('example');
else
    disp('Simulation terminated successfully.')
end