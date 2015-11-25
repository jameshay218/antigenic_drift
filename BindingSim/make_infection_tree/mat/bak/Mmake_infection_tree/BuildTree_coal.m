function void = BuildTree_coal(filename_simData, n_seqs)

load(filename_simData);

outfile_treeData = strcat(filename_simData, '_coal_treeData_', int2str(n_seqs));           
outfile_tree = strcat(filename_simData, '_coal_genealogy_', int2str(n_seqs), '.tree');

n_tot_samples = n_seqs;

seq_times = (GetIndividualsSampled_coal(timeSeries, n_tot_samples))';
indiv_sampled = (1:n_seqs)';

for i = 1:n_tot_samples
    names{i} = strcat('sample', int2str(i), '_', num2str(seq_times(i)));
end

locs_max = find(seq_times == max(seq_times));
n_curr_lin = length(locs_max);
last_seq = n_seqs;
b = NaN*ones((n_seqs - 1),2); loc_b = 1;
d = NaN*ones(2*n_seqs - 1, 1); d(end) = 0;

t_curr = seq_times(locs_max(1));
listNodes = locs_max;
dateNodes = seq_times(locs_max);
seq_times(locs_max) = [];

locs_begin = find(timeSeries.t_vals > t_curr);
curr_time_loc = locs_begin(1);

epi_params.beta = epi_params.mean_R0/epi_params.mean_durationInfection;

while curr_time_loc > 1

    if isempty(find(isnan(d)))
        break;
    end
    
    lambda_curr = n_curr_lin*(n_curr_lin-1)*(epi_params.beta*timeSeries.S_vals(curr_time_loc)/epi_params.N)*(1 + 1/epi_params.kbeta)/timeSeries.I_vals(curr_time_loc); 
    coal_time = exprnd(1/lambda_curr);
    
    coal_time_real = t_curr - coal_time;
    
    %timeSeries.t_vals(curr_time_loc - 1)
    
    %max(seq_times)

    event = find([coal_time_real timeSeries.t_vals(curr_time_loc - 1) max(seq_times)] == max([coal_time_real timeSeries.t_vals(curr_time_loc - 1) max(seq_times)]));
    
    switch event(1)
        case 1
            % coalescent event occurred- figure out which lineages coalesced and update b and d appropriately!
            locs = randperm(n_curr_lin);
            last_seq = last_seq + 1;
            t_curr = coal_time_real;
        
            listNodes = [listNodes last_seq];
            dateNodes = [dateNodes t_curr];
        
            b(loc_b, :) = listNodes(locs(1:2));
            d(listNodes(locs(1))) = dateNodes(locs(1)) - t_curr;
            d(listNodes(locs(2))) = dateNodes(locs(2)) - t_curr;
        
            listNodes(locs(1:2)) = [];
            dateNodes(locs(1:2)) = [];
            n_curr_lin = n_curr_lin - 1
            loc_b = loc_b + 1;
            
            t_curr = coal_time_real;
            
        case 2
            curr_time_loc = curr_time_loc - 1;
            t_curr = timeSeries.t_vals(curr_time_loc);
            
        case 3
            
            locs_max = find(seq_times == max(seq_times));
            n_curr_lin = n_curr_lin + length(locs_max);
            
%             listNodes
%             dateNodes
%             seq_times(locs_max)
            
            t_curr = seq_times(locs_max(1));
            listNodes = [listNodes locs_max'];
            dateNodes = [dateNodes seq_times(locs_max)'];
            seq_times(locs_max) = [];
           
        otherwise 
            error('event needs to be 1, 2, or 3')
    end
           
end


% % % the remaining lineages do not coalesce
% % % figure out roughly the lambda
% % lambda_approx = 1/min([d(indiv1_index(end), 1), d(indiv2_index(end), 1)])
% % K_approx = lambda_approx/((n_seqs + 1)*n_seqs);
% % 
% % % make it coalesce more rapidly!
% % K_approx = K_approx*10; % make it coalesce more rapidly now!


% link together arbitrarily
while 1
    
    if isempty(find(isnan(d)))
        break;
    end
    
%     lambda_curr = n_curr_lin*(n_curr_lin-1)*(epi_params.beta*timeSeries.S_vals(1)/epi_params.N)*(1 + 1/epi_params.kbeta)/timeSeries.I_vals(1); 
%     coal_time = exprnd(1/lambda_curr);
%     
%     t_curr
%     
%     coal_time_real = t_curr - coal_time
    %coal_time_real = -50; % just have all the remaining coalesce at time t = -50
    coal_time_real = time.start-50; % just have all the remaining coalesce at time t = -50
        
%     timeSeries.t_vals(1)
%     
%     max(seq_times)
% 
%     event = find([coal_time_real timeSeries.t_vals(curr_time_loc - 1) max(seq_times)] == max([coal_time_real timeSeries.t_vals(curr_time_loc - 1) max(seq_times)]))
%     
%     
    % coalescent event occurred- figure out which lineages coalesced and update b and d appropriately!
    locs = randperm(n_curr_lin);
    last_seq = last_seq + 1;
    t_curr = coal_time_real;
        
    listNodes = [listNodes last_seq];
    dateNodes = [dateNodes t_curr];
        
    b(loc_b, :) = listNodes(locs(1:2));
    d(listNodes(locs(1))) = dateNodes(locs(1)) - t_curr;
    d(listNodes(locs(2))) = dateNodes(locs(2)) - t_curr;
        
    listNodes(locs(1:2)) = [];
    dateNodes(locs(1:2)) = [];
    n_curr_lin = n_curr_lin - 1
    loc_b = loc_b + 1;
            
    t_curr = coal_time_real;
            
    
end
    
d(end) = 0;
b
d

tree = phytree(b,d, names);
view(tree)

save(outfile_treeData);
phytreewrite(outfile_tree, tree, 'BranchNames', false)
