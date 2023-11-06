function run_suppl_fig7(option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script generates/replicates the Suppl. Figure 7.
% %
% INPUT
% - option: integer 1 to generate item with existing results.
%           integer 2 to recreate item from original data, involving all
%           required computations.
% OUTPUT
% Suppl. Figure 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 1
    error('Please provide an option (1 or 2) as an argument.');
end

if option == 1
    disp('Option 1: Generate item with existing results');
    cd './data/scripts/'
    plot_suppl_fig7;
    
elseif option == 2
    disp('Option 2: Recreate item from original data');
    cd './data_replicated/scripts/'
    
    tic;
    disp('Running: Building Yeast DIP network')
    create_Scere_DIP_net;
    elapsed_time = toc;
    
    disp(['Elapsed time for Yeast DIP creation: ', num2str(elapsed_time), ' s']);
    
    disp('Running: AFM-IS and AFM-piTMS scores processing')
    tic;
    system('python protein_pairs_ints_data_processing.py');
    system('python protein_pairs_pitms_data_processing.py');
    elapsed_time = toc;
    
    disp(['Elapsed time for AFM-IS and AFM-piTMS scores processing: ', num2str(elapsed_time), ' s']);
    
    disp('Plotting Suppl. Figure 7')
    plot_suppl_fig7;
    
else
    error('Invalid option. Please choose 1 or 2.');
end

end