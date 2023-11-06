function run_figure_4(option)
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
    cd './data/scripts/'
    plot_figure_4;
    
elseif option == 2
    disp('Option 2: Recreate item from original data');
    cd './data_replicated/scripts/'
    
    tic;
    disp('Running: Building networks')
    create_Scere_DIP_net;

    elapsed_time = toc;
    
    disp(['Elapsed time for building networks: ', num2str(elapsed_time), ' s']);
    
    disp('Running: Creation of list of protein pairs to predict and 1% Yeast DIP network perturbation')
    tic;
    create_list_pairs_network_perturbed_1percent_Yeast;
    elapsed_time = toc;
    
    disp(['Elapsed time for Creating list of protein pairs and 1% Yeast DIP network perturbation: ', num2str(elapsed_time), ' s']);
    
    disp('Running: CH link prediction for Yeast DIP')
    tic;
    run_CHA_Yeast_DIP_net;
    elapsed_time = toc;
    
    disp(['Elapsed time for calculating CH similarity scores for Yeast DIP: ', num2str(elapsed_time), ' s']);
    
    disp('Running: AFM-IS and AFM-piTMS scores processing')
    tic;
    system('python protein_pairs_ints_data_processing.py');
    system('python protein_pairs_pitms_data_processing.py');
    elapsed_time = toc;
    
    disp(['Elapsed time for AFM-IS and AFM-piTMS scores processing: ', num2str(elapsed_time), ' s']);
    
    disp('Running: IUPRED3')
    tic;
    system('python iupred3_AF2_Yeast.py');
    elapsed_time = toc;
    
    disp(['Elapsed time for IUPRED3: ', num2str(elapsed_time), ' s']);

    disp('Plotting: Figure 4')
    tic;
    plot_figure_4;
    elapsed_time = toc;
    
    disp(['Elapsed time for plotting Figure 4: ', num2str(elapsed_time), ' s']);
      
else
    error('Invalid option. Please choose 1 or 2.');
end

end