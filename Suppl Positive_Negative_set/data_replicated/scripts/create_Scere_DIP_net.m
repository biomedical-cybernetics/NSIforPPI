% MATLAB Script Documentation
% Author: Ilyes Abdelhamid, 2022
% Description: This MATLAB script reads a table of protein interactions
% from a file, extracts specific information from it, and saves the data in a desired format.

path_output_net = '../matrix/';
path_output_table = '../table/';
path_filename = '../../data//original_data/';
filename = 'Scere20170205.txt';
dt = readtable([path_filename, 'Scere20170205.txt'], 'Delimiter', '\t', 'ReadVariableNames', false);

% Column 1 of dt is ID interactor A
ID_interactor_A = 1;
% Column 2 of dt is ID interactor B
ID_interactor_B = 2;
% Column 10 of dt is taxonomy ID of interactor A 
Taxid_interactor_A = 10;
% Column 11 of dt is taxonomy ID of interactor B
Taxid_interactor_B = 11;
% Taxonomy ID of yeast (Saccharomyces cerevisiae)
Taxid_Scere = 'taxid:4932(Saccharomyces cerevisiae)';
dt = table2array(dt(:, [ID_interactor_A ID_interactor_B Taxid_interactor_A Taxid_interactor_B]));
list_uniprot_id1 = {};
list_uniprot_id2 = {};
for i = 1:size(dt,1)
    % protein A and protein B must have a Uniprot entry for better
    % reliability.
    % Researchers often turn to UniProtKB when they need reliable and
    % detailed information about individual protein sequences,
    if strfind(dt{i,1}, 'uniprotkb') & strfind(dt{i,2}, 'uniprotkb')
        % protein A and protein B must be assigned to the same correct "taxid" 
        if strcmp(Taxid_Scere, dt{i,3}) && strcmp(Taxid_Scere, dt{i,4})
            % Extract Uniprot ID for protein A
            splitStr = regexp(dt{i,1},'\|','split');
            id1 = splitStr{1};
            splitStr_uniprot = regexp(dt{i,1},'.*uniprotkb:','split');
            uniprot_id1 = splitStr_uniprot{2};
            
            % Extract Uniprot ID for protein B
            splitStr = regexp(dt{i,2},'\|','split');
            id2 = splitStr{1};
            splitStr_uniprot = regexp(dt{i,2},'.*uniprotkb:','split');
            uniprot_id2 = splitStr_uniprot{2};
            
            % Check if Protein A is different from Protein B (heterodimeric interaction)
            if ~strcmp(id1, id2)
                % Append Uniprot IDs to their respective vectors
                list_uniprot_id1 = [list_uniprot_id1(:); uniprot_id1];
                list_uniprot_id2 = [list_uniprot_id2(:); uniprot_id2];
            end
        end
        
    else
        continue
    end
    
end

% Create a matrix (x) and a map for the first component.
[x, map] = create_matrix(list_uniprot_id1(:,1), list_uniprot_id2(:,1), [], 1);
x_lcc = x{1};
map_lcc = map{1};

% save unique uniprot ids in csv file.
%unique_uniprot_ids = unique([list_uniprot_id1; list_uniprot_id2]);
fid = fopen([path_output_table, 'Scere_uniprot_ids.csv'],'w');
fprintf(fid,'%s\n', map_lcc{:});
fclose(fid);

% save network.
save([path_output_net, 'Yeast_DIP_net.mat'], 'x_lcc', 'map_lcc', '-v7.3')

