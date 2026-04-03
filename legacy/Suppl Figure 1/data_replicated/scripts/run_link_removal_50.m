% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2023
% Description: This MATLAB script creates sparsified matrices by removing
% at random a certain percentage of links.


networks = {'Arabidopsis', 'Arabidopsis_BioGRID', 'Bioplex', 'Fly_BioGRID', ...
    'Hein', 'HI_II_14', 'HI_III', 'Human_BioGRID', 'Interactome3D', 'Lit_BM_13', ...
    'Mouse_BioGRID', 'S_Pombe_BioGRID', 'Worm_BioGRID', 'Yeast', 'Yeast_BioGRID'};

perc = 0.5;
iters = 10;
method = 'rand';

for i = 1:length(networks)
    if ~exist(['../sparsified_matrices/' networks{i} '_linkrem50.mat'],'file')
        display(networks{i})
        load(['../matrix/' networks{i} '.mat'], 'x')
        [matrices, ~] = random_link_removal(x, perc, iters, method);
        save(['../sparsified_matrices/' networks{i} '_linkrem50.mat'], 'matrices')
    end
end