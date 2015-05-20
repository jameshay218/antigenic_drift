function void = BuildTree_indiv_nexus(filename_infectionTreeData, n_seqs)
% Produce the tree file as nexus format
global epi_params;

%load(filename_simData); % Is this line needed? 20130716
load(filename_infectionTreeData);
who

folder = char(regexp(filename_infectionTreeData,'.+/','match'));
outfile_treeData = strcat(folder, 'indiv_treeData_', int2str(n_seqs));            
%outfile_tree = strcat('out/', filename_simData, '_indiv_genealogy_', int2str(n_seqs), '.tree');
outfile_tree = strcat(folder, 'indiv_genealogy_', int2str(n_seqs), '.tree');
n_tot_samples = n_seqs;

% seq_times are the death times of each individuals
[seq_times, indiv_sampled] = GetIndividualsSampled_indiv(epi_params, births, deaths, n_tot_samples);

% Should retrieve viruses binding avidities Jul 17, 2013
% Access virus_traits.mat
indiv_sampled_binding = binding(indiv_sampled);
for i = 1:n_tot_samples
    header = num2str(i);
    comment = ['[&annotation="",Netcharge.set="",Bindingscr.set="' num2str(indiv_sampled_binding(i)) '"]'];
    nexus_names{i} = strcat('sample', int2str(i), '_' ,num2str(indiv_sampled(i)), '_', num2str(seq_times(i)));
    names{i} = [header comment];
end

% The ancestral lineages of all sampled nodes
parentLineages = GetIndividualsLineages_indiv(n_seqs, indiv_sampled, parent);

b = NaN*ones((n_seqs - 1),2);
d = NaN*ones(2*n_seqs - 1, 1);
loc_b = 1;

curr_indiv = indiv_sampled;
complete_indiv = indiv_sampled;
complete_seq_times = seq_times;

while length(curr_indiv)>0
    
    % curr_indiv: sampled individuals
    % parentLineages: parental lineages of sampled individuals   
    % coal_daughters: two external tips from sampled individuals
    % coal_parent: internal node from unsampled individuals
    [coal_daughters, coal_parent, timeOfCoalescence] = FindMostRecentCoalescence_indiv(curr_indiv, parentLineages, births, deaths);
    
    if timeOfCoalescence == -Inf
        break;
    end
    
    indiv1_index = find(complete_indiv == coal_daughters(1));
    indiv2_index = find(complete_indiv == coal_daughters(2));

    % update distance matrix
    d(indiv1_index(end), 1) = complete_seq_times(indiv1_index(end)) - timeOfCoalescence;
    d(indiv2_index(end), 1) = complete_seq_times(indiv2_index(end)) - timeOfCoalescence;

    % update branch matrix
    b(loc_b, 1:2) = [indiv1_index(end) indiv2_index(end)];
    
    % create internal nodes
    %names{n_tot_samples + loc_b} = strcat('node', int2str(loc_b), '_', num2str(coal_parent) ,'_' ,num2str(timeOfCoalescence)); % add to name
    header = num2str(n_tot_samples+loc_b);
    comment = strcat('[&annotation="",Netcharge.set="",Bindingscr.set="', num2str(binding(coal_parent)), '"]');
    %nexus_names{i} = strcat('node', int2str(i), '_' ,num2str(indiv_sampled(i)), '_', num2str(seq_times(i)));
    names{n_tot_samples + loc_b} = [header comment];
    
    complete_indiv = [complete_indiv coal_parent]; % add to individual list
    complete_seq_times = [complete_seq_times timeOfCoalescence]; % add to seq time list
    
    loc_b = loc_b + 1;
    

    
    % erase coal_daughters from curr_indiv
    loc1_inCurr = find(curr_indiv == coal_daughters(1));
    curr_indiv(loc1_inCurr(1)) = [];
    loc2_inCurr = find(curr_indiv == coal_daughters(2));
    curr_indiv(loc2_inCurr(1)) = [];
    
    % add coal_parent to curr_indiv list
    curr_indiv = [curr_indiv coal_parent];
    
    % redo parentLineages
    n_seqs = length(curr_indiv)
    clear parentLineages
    parentLineages = GetIndividualsLineages_indiv(n_seqs, curr_indiv, parent);        
    
end


% the remaining lineages do not coalesce
% figure out roughly the lambda
%lambda_approx = 1/min([d(indiv1_index(end), 1), d(indiv2_index(end), 1)])
%K_approx = lambda_approx/((n_seqs + 1)*n_seqs);

% make it coalesce more rapidly!
%K_approx = K_approx*10; % make it coalesce more rapidly now!
 
% link together arbitrarily
while 1
    if n_seqs == 1
        break;
    end
    
    perm_indiv = randperm(length(curr_indiv));
    coal_daughters = curr_indiv(perm_indiv(1:2));
    
    indiv1_index = find(complete_indiv == coal_daughters(1));
    indiv2_index = find(complete_indiv == coal_daughters(2));

    %branch_length = exprnd(1/(K_approx*n_seqs*(n_seqs-1)));
    %timeOfCoalescence = min([complete_seq_times(indiv1_index) complete_seq_times(indiv2_index)]) - branch_length;
    recovery_rate = 1/5; % in days
    num_infected = 7000;
    lambda_val = n_seqs*(n_seqs-1)*recovery_rate/num_infected;
    branch_length = exprnd(1/lambda_val);
    timeOfCoalescence = min([complete_seq_times(indiv1_index) complete_seq_times(indiv2_index)]) - branch_length;
    
    who
    
    %timeOfCoalescence = -50; % just have all the remaining coalesce at time t = -50
    %timeOfCoalescence = time.start-50; % just have all the remaining coalesce at time t = -50
    %timeOfCoalescence = epi_params.tRange_stoch(1)-50; 
    
    d(indiv1_index, 1) = complete_seq_times(indiv1_index) - timeOfCoalescence;
    d(indiv2_index, 1) = complete_seq_times(indiv2_index) - timeOfCoalescence;
    
    if (find(d<0))
      disp 'error: distance is negative.';
      %disp 'change distance to zero if it is negative';
      %d(find(d<0))=0;
    end
    b(loc_b, 1:2) = [indiv1_index(end) indiv2_index(end)];
    % timeOfCoalescence < 0
    header = num2str(n_tot_samples+loc_b);
    comment = strcat('[&annotation="",Netcharge.set="",Bindingscr.set=""]');
    names{n_tot_samples + loc_b} = [header comment];
    %nexus_names{i} = strcat('node', int2str(loc_b), '_0_', num2str(timeOfCoalescence));
    loc_b = loc_b + 1;
    
    coal_parent = max(complete_indiv) + 1;
    
    complete_indiv = [complete_indiv coal_parent];
    complete_seq_times = [complete_seq_times timeOfCoalescence];
    
    % erase coal_daughters from curr_indiv
    loc1_inCurr = find(curr_indiv == coal_daughters(1));
    curr_indiv(loc1_inCurr(1)) = [];
    loc2_inCurr = find(curr_indiv == coal_daughters(2));
    curr_indiv(loc2_inCurr(1)) = [];
    
    % add coal_parent to curr_indiv list
    curr_indiv = [curr_indiv coal_parent];
    
    % redo parentLineages
    n_seqs = length(curr_indiv)
    
end
    
d(end) = 0;
% output tree array
save([folder 'tree_' num2str(n_tot_samples) '.mat'], 'b','d', 'names');

if (find(d<0))
      disp 'error: distance is negative.';
      disp 'change distance to zero if it is negative';
      d(find(d<0))=0;
end

tree = phytree(b,d, names);
view(tree)


% output tree data
save(outfile_treeData);
phytreewrite(outfile_tree, tree, 'BranchNames', false);

% output branch lengths
pos_s = regexp(outfile_tree,'\.');
treeoutfile = [outfile_tree(1:pos_s-1) '.branch' outfile_tree(pos_s:end)];

% write phytree
phytreewrite(treeoutfile, tree, 'BranchNames', true);

% write nexus tree
pos_s = regexp(outfile_tree,'\.');
nexusoutfile = [outfile_tree(1:pos_s-1) '.nx.branch' outfile_tree(pos_s:end)];

fid=fopen(nexusoutfile, 'w');
seq = ['#NEXUS\n\n' 'Begin taxa;\n' 'Dimensions ntax=' num2str(n_tot_samples) ';\n' 'Taxlabels\n'];
fprintf(fid, seq);
for i=1:length(nexus_names)
  fprintf(fid, [nexus_names{i} '\n']);
end
seq = [';\n' 'End;\n\n' 'Begin trees;\n' 'Translate\n'];
fprintf(fid, seq);
for i=1:length(nexus_names)
    if i<length(nexus_names)
      fprintf(fid, [num2str(i) ' ' nexus_names{i} ',\n']);
    else
      fprintf(fid, [num2str(i) ' ' nexus_names{i} '\n']);
    end
end
seq = [';\n' 'tree TREE1 = [&R]'];
fprintf(fid, seq);

% append phytree into nexus file 
tree_str=fileread(treeoutfile);
expression = '\';
replace = '';
new_tree_str = regexprep(tree_str,expression,replace);
fprintf(fid, new_tree_str);
fprintf(fid, ['\nEnd;\n']);
fclose all;





