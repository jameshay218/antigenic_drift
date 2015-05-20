function void = BuildTree_indiv_v2(filename_simData, filename_infectionTreeData, n_seqs)
global epi_params;

%load(filename_simData); % Is this line needed? 20130716
load(filename_infectionTreeData);
who

folder = char(regexp(filename_simData,'.+/','match'));
outfile_treeData = strcat(folder, 'indiv_treeData_', int2str(n_seqs));            
%outfile_tree = strcat('out/', filename_simData, '_indiv_genealogy_', int2str(n_seqs), '.tree');
outfile_tree = strcat(folder, 'indiv_genealogy_', int2str(n_seqs), '.tree');
n_tot_samples = n_seqs;

[seq_times, indiv_sampled] = GetIndividualsSampled_indiv_v2(epi_params, births, deaths, n_tot_samples);

for i = 1:n_tot_samples
    names{i} = strcat('sample', int2str(i), '_', num2str(seq_times(i)));
end

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
    % coal_daughters: two external tips
    % coal_parent: internal node
    [coal_daughters, coal_parent, timeOfCoalescence] = FindMostRecentCoalescence_indiv(curr_indiv, parentLineages, births, deaths);
    
    if timeOfCoalescence == -Inf
        break;
    end
    
    indiv1_index = find(complete_indiv == coal_daughters(1));
    indiv2_index = find(complete_indiv == coal_daughters(2));

    d(indiv1_index(end), 1) = complete_seq_times(indiv1_index(end)) - timeOfCoalescence;
    d(indiv2_index(end), 1) = complete_seq_times(indiv2_index(end)) - timeOfCoalescence;

    b(loc_b, 1:2) = [indiv1_index(end) indiv2_index(end)];
    
    % create internal nodes
    names{n_tot_samples + loc_b} = strcat('node', int2str(loc_b), '_', num2str(timeOfCoalescence)); % add to name
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
lambda_approx = 1/min([d(indiv1_index(end), 1), d(indiv2_index(end), 1)])
K_approx = lambda_approx/((n_seqs + 1)*n_seqs);

% make it coalesce more rapidly!
K_approx = K_approx*10; % make it coalesce more rapidly now!
 
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
     
    who
    
    %timeOfCoalescence = -50; % just have all the remaining coalesce at time t = -50
    %timeOfCoalescence = time.start-50; % just have all the remaining coalesce at time t = -50
    timeOfCoalescence = epi_params.tRange_stoch(1)-50; 
    
    d(indiv1_index, 1) = complete_seq_times(indiv1_index) - timeOfCoalescence;
    d(indiv2_index, 1) = complete_seq_times(indiv2_index) - timeOfCoalescence;
    
    if (find(d<0))
      disp 'error: distance is negative';
    end
    b(loc_b, 1:2) = [indiv1_index(end) indiv2_index(end)];
    names{n_tot_samples + loc_b} = strcat('node', int2str(loc_b), '_', num2str(timeOfCoalescence));
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


tree = phytree(b,d, names);
view(tree)

save(outfile_treeData);
phytreewrite(outfile_tree, tree, 'BranchNames', false);

pos_s = regexp(outfile_tree,'\.');
phytreewrite([outfile_tree(1:pos_s-1) '.branch' outfile_tree(pos_s:end)], tree, 'BranchNames', true);

