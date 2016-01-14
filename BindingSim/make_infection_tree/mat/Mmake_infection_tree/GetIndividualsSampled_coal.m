function seq_times = GetIndividualsSampled_coal(timeSeries, n_tot_samples)

if 0
    % set sampling time points, and pick random individuals who are
    % infected at those sampling time points.
    % this breaks down sometimes-- when there is a coalescence between a
    % sampled individual and its parent and its parent is also sampled, but at a time that's 
    % earlier than the coalescence time 
    
    % times of isolation are chosen according to incidence level
    seq_times = GetTimesOfIsolation(timeSeries, n_tot_samples);
    
    locsNotPrevUsed = 1:length(births);
    for i = 1:n_tot_samples
        locsBornOK = union(find(births < seq_times(i)), find(isnan(births)));
        locsDeathOK = find(deaths > seq_times(i));
        locsOK = intersect(intersect(locsBornOK, locsDeathOK), locsNotPrevUsed);
        randpermlist = randperm(length(locsOK));
        indiv_sampled(i) = locsOK(randpermlist(1));
        locsNotPrevUsed = setdiff(locsNotPrevUsed, indiv_sampled(i));
    end
else
    % set sampling time points by time of recovery
    
    % times of isolation depend on recoveries
    n_time_intervals = length(timeSeries.R_vals);
    %prop_recoveries = timeSeries.R_vals/timeSeries.R_vals(end);
    rand_nums = rand(n_tot_samples, 1);

    % SAMPLE ALL SEQUENCE ON OR ABOUT time point = 225!!!
    %locs = find(timeSeries.t_vals >= 225);

    %locsNotPrevUsed = 1:length(births);
    for i = 1:n_tot_samples
        [i n_tot_samples]
        %locs = find(prop_incidence > rand_nums(i)); % pick times
        %proportional to recoveries
   
        % WE'RE PICKING UNIFORMLY FROM THE EPIDEMIC RIGHT NOW!!!
        locs = randi(n_time_intervals);  % pick uniformly from time during epidemic
   
        %date_of_isolation = timeSeries.t_vals(locs(1));
        
        seq_times(i) = timeSeries.t_vals(locs(1));
        
        % this is the ideal situation: pick from very last bit of the epidemic!!
        %locs = round(0.75*n_time_intervals + (n_time_intervals - (0.75*n_time_intervals))*rand);
%    
%         locsBornOK = union(find(births < date_of_isolation), find(isnan(births)));
%         locsDeathOK = find(deaths > date_of_isolation);
%         locsOK = intersect(intersect(locsBornOK, locsDeathOK), locsNotPrevUsed);
%         
%         % choose the individual who has died earliest (so closest to the
%         % sampling time
%         loc_indiv = find(deaths(locsOK) == min(deaths(locsOK)));
%         
%         randperm_array = randperm(length(loc_indiv));
%         
%         indiv_sampled(i) = locsOK(loc_indiv(randperm_array(1)));
%         seq_times(i) = deaths(indiv_sampled(i));
%         
%         locsNotPrevUsed = setdiff(locsNotPrevUsed, indiv_sampled(i));
%         
    end
end

[seq_times, indices] = sort(seq_times);
%indiv_sampled = indiv_sampled(indices);
    