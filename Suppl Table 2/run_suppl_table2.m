function run_suppl_table2(option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script generates/replicates the Suppl. Table 2.
% %
% INPUT
% - option: integer 1 to generate item with existing results.
%           integer 2 to recreate item from original data, involving all
%           required computations.
% OUTPUT
% Suppl. Table 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 1
    error('Please provide an option (1 or 2) as an argument.');
end

if option == 1
    disp('Option 1: Generate item with existing results');
    cd './data/script/'
    create_suppl_table2;
    
elseif option == 2
    disp('Option 2: Recreate item from original data');
    cd './data_replicated/scripts/'
    
    tic;
    disp('Running: Building networks')
    create_Arabidopsis;
    create_Arabidopsis_BioGRID;
    create_Bioplex;
    create_Fly_BioGRID;
    create_Hein;
    create_HI_II_14;
    create_HI_III;
    create_Human_BioGRID;
    create_Interactome3D;
    create_Mouse_BioGRID;
    create_S_Pombe_BioGRID;
    create_Worm_BioGRID;
    create_Yeast;
    create_Yeast_BioGRID;
    create_Lit_BM_13;
    elapsed_time = toc;
    
    disp(['Elapsed time for building networks: ', num2str(elapsed_time), ' s']);
    
    tic;
    disp('Running: Link Removal 10%')
    run_link_removal_10;
    elapsed_time = toc;
    
    disp(['Elapsed time for link removal 10%: ', num2str(elapsed_time), ' s']);
    
    disp('Running: CH-based scores')
    tic;
    run_simulation_linkrem10_CHv2_L2_L3;
    elapsed_time = toc;
      
    disp(['Elapsed time for CH link predicition: ', num2str(elapsed_time), ' s']);
     
    disp('Plotting Suppl. Table 2')
    create_suppl_table2;
    
else
    error('Invalid option. Please choose 1 or 2.');
end

end