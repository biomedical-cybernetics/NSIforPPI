sets = {'Positive_set', 'Negative_set'};
filename = 'DIP_20171201.fasta';
fasta_file = fastaread(['../original_data/', filename]);
fasta_cell = struct2cell(fasta_file)';
organism = 'Yeast';
tic;
for i = 1:length(sets)
    save_files_dir = ['../FASTA_files/' organism '/' organism '_' sets{i} '/'];
    P = tdfread(['../' sets{i} '_Yeast.txt'], 'tab');
    protein_1 = cellstr(P.Protein_1);
    protein_2 = cellstr(P.Protein_2);
    num_pairs_to_test = length(protein_1);
    no_seq = 0;
    for j = 1:num_pairs_to_test
        idx_1 = find(contains(fasta_cell(:,1),['uniprot:' protein_1{j}]));
        idx_2 = find(contains(fasta_cell(:,1),['uniprot:' protein_2{j}]));
        if isempty(idx_1) || isempty(idx_2)
            no_seq = no_seq + 1;
            continue;
        end
        data(1).Sequence = fasta_cell{idx_1,2};
        data(1).Header =  protein_1{j};
        data(2).Sequence = fasta_cell{idx_2,2};
        data(2).Header = protein_2{j};
        
        if exist([save_files_dir num2str(j) '_' protein_1{j} '_' protein_2{j} '.fasta'],'file')
            continue;
        else
            fastawrite([save_files_dir num2str(j) '_' protein_1{j} '_' protein_2{j} '.fasta'], data);
        end
        
    end
    
end
toc;