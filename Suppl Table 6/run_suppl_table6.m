function run_suppl_table6(option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script generates/replicates the Suppl. Table 6.
% %
% INPUT
% - option: integer 1 to generate item with existing results.
%           integer 2 to recreate item from original data, involving all
%           required computations.
% OUTPUT
% Suppl. Table 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 1
    error('Please provide an option (1 or 2) as an argument.');
end

if option == 1
    disp('Option 1: Generate item with existing results');
    cd './data/script/'
    create_suppl_table6;
    
elseif option == 2
    disp('Option 2: Recreate item from original data');
    cd './data_replicated/scripts/'
    
    tic;
    disp('Running: Building networks')
    create_Scere_DIP_net;

    elapsed_time = toc;
    
    disp(['Elapsed time for building networks: ', num2str(elapsed_time), ' s']);
    
    disp('Running: topological measures')
    tic;
    run_compute_statistics_original;
    elapsed_time = toc;
    
    disp(['Elapsed time for calculating topological measures: ', num2str(elapsed_time), ' s']);
    
    disp('Creating: statistics table')
    create_suppl_table6;
      
else
    error('Invalid option. Please choose 1 or 2.');
end

end