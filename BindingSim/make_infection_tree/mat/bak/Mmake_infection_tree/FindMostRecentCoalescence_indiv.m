function [coal_samples, coal_parent, mostRecentTimeOfCoalescence] = FindMostRecentCoalescence_indiv(curr_indiv, parentLineages, births, deaths)

coal_samples = [NaN NaN];
coal_parent = NaN;
mostRecentTimeOfCoalescence = -Inf;

n_individuals_sampled = length(curr_indiv);

for i = 1:n_individuals_sampled
    for j = (i+1):n_individuals_sampled
        parents_in_common = setdiff(intersect(parentLineages(i).lin, parentLineages(j).lin), NaN);
        if ~isempty(parents_in_common)
            this_parent = max(parents_in_common);
            loc1 = find(parentLineages(i).lin == this_parent);
            loc2 = find(parentLineages(j).lin == this_parent);

            if curr_indiv(i) == this_parent
                coalTime1 = deaths(curr_indiv(i));
            else
                coalTime1 = births(parentLineages(i).lin(loc1 - 1));
            end
            
            if curr_indiv(j) == this_parent
                coalTime2 = deaths(curr_indiv(j));
            else
                coalTime2 = births(parentLineages(j).lin(loc2 - 1));
            end
            
            thisTimeOfCoalesence = min(coalTime1, coalTime2);
            
            if thisTimeOfCoalesence > mostRecentTimeOfCoalescence
                coal_samples = [curr_indiv(i) curr_indiv(j)];
                coal_parent = this_parent;
                mostRecentTimeOfCoalescence = thisTimeOfCoalesence;
            end
        end
    end
end
