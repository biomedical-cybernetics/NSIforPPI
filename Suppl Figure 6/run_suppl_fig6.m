function run_suppl_fig6(option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script generates/replicates the Suppl. Figure 5.
% %
% INPUT
% - option: integer 1 to generate item with existing results.
%           integer 2 to recreate item from original data, involving all
%           required computations.
% OUTPUT
% Suppl. Figure 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 1
    error('Please provide an option (1 or 2) as an argument.');
end

if option == 1
    disp('Option 1: Generate item with existing results');
    cd './data/script/'
    plot_suppl_fig6;
    
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
    
    disp('Running: Creation of list of protein pairs to predict and 1% Yeast DIP network perturbation')
    tic;
    create_list_pairs_network_perturbed_1percent_Yeast;
    elapsed_time = toc;
    
    disp(['Elapsed time for Creating list of protein pairs and 1% Yeast DIP network perturbation: ', num2str(elapsed_time), ' s']);
    
    tic;
    disp('Running: CH-based scores')
    run_CHA_linkpred1_monopartite_Yeast_net_perturbed;
    elapsed_time = toc;
    
    disp(['Elapsed time for CH link prediction: ', num2str(elapsed_time), ' s']);
    % 800.9407 s
    
    disp('Plotting Suppl. Figure 6')
    plot_suppl_fig6;
    
else
    error('Invalid option. Please choose 1 or 2.');
end

end