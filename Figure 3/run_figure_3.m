function run_figure_3(option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script generates/replicates the Figure 3.
% %
% INPUT
% - option: integer 1 to generate item with existing results.
%           integer 2 to recreate item from original data, involving all
%           required computations.
% OUTPUT
% Figure 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 1
    error('Please provide an option (1 or 2) as an argument.');
end

if option == 1
    disp('Option 1: Generate item with existing results');
    cd './data/script/'
    plot_figure_3;
    
elseif option == 2
    disp('Option 2: Recreate item from original data');
    cd './data_replicated/scripts/'
    
    disp('Running: Building Yeast DIP network')
    tic;
    create_Scere_DIP_net;
    elapsed_time = toc;
    
    disp(['Elapsed time for Yeast DIP creation: ', num2str(elapsed_time), ' s']);
    
    disp('Running: AFM-IS and AFM-piTMS scores processing')
    tic;
    system('python protein_pairs_ints_data_processing.py');
    system('python protein_pairs_pitms_data_processing.py');
    elapsed_time = toc;
    
    disp(['Elapsed time for AFM-IS and AFM-piTMS scores processing: ', num2str(elapsed_time), ' s']);
    
    disp('Running: Creation of list of protein pairs to predict and 10% Yeast DIP network perturbation')
    tic;
    create_list_GSP_GSN_and_perturbed_net_10perc;
    elapsed_time = toc;
    
    disp(['Elapsed time for Creating list of protein pairs and 10% Yeast DIP network perturbation: ', num2str(elapsed_time), ' s']);
    
    disp('Running: CH-baseed scores')
    tic;
    run_CHA_linkpred10_monopartite_Yeast_net_perturbed;
    elapsed_time = toc;
    
    disp(['Elapsed time for CH link prediction: ', num2str(elapsed_time), ' s']);
    
    disp('Plotting Figure 3')
    tic;
    plot_figure_3;
    elapsed_time = toc;
    
    disp(['Elapsed time for Figure 3: ', num2str(elapsed_time), ' s']);
    
else
    error('Invalid option. Please choose 1 or 2.');
end

end