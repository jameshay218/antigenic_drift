function [seq_times, indiv_sampled] = GetIndividualsSampled_indiv(epi_params, births, deaths, n_tot_samples)
indiv_sampled = [];
if 1
    % set sampling time points, and pick random individuals who are
    % infected at those sampling time points.
    % this breaks down sometimes-- when there is a coalescence between a
    % sampled individual and its parent and its parent is also sampled, but at a time that's 
    % earlier than the coalescence time 
    
    % times of isolation are chosen according to incidence level
    seq_times = epi_params.tRange_stoch(1) + (epi_params.tRange_stoch(2) - epi_params.tRange_stoch(1))*rand(1, n_tot_samples);
    %seq_times = GetTimesOfIsolation(timeSeries, n_tot_samples);
    
    locsNotPrevUsed = 1:length(births);
    %random stream used for randperm
    s = RandStream('mt19937ar','Seed',sum(100*clock));
    for i = 1:n_tot_samples
        %i
        locsBornOK = union(find(births < seq_times(i)), find(isnan(births)));
        locsDeathOK = find(deaths > seq_times(i));
        locsOK = intersect(intersect(locsBornOK, locsDeathOK), locsNotPrevUsed);
        randpermlist = randperm(s,length(locsOK));
        indiv_sampled(i) = locsOK(randpermlist(1));
        locsNotPrevUsed = setdiff(locsNotPrevUsed, indiv_sampled(i));
    end
end

[seq_times, indices] = sort(seq_times);
indiv_sampled = indiv_sampled(indices);
    