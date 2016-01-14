function [sortedBirths, sortedDeaths, sortedParents] = GetInfectionTree(filename)

load(filename);

C_vals_floor = floor(timeSeries.C_vals);
new_infections = diff(C_vals_floor);

R_vals_floor = floor(timeSeries.R_vals);
new_recoveries = diff(R_vals_floor);

mean_beta = epi_params.mean_R0/epi_params.mean_durationInfection;
scale_parameter = mean_beta/epi_params.kbeta;

births(1:round(epi_params.I_init), 1) = NaN;
deaths(1:round(epi_params.I_init), 1) = NaN;
parent(1:round(epi_params.I_init), 1) = NaN;
beta_val(1:round(epi_params.I_init), 1) = gamrnd(epi_params.kbeta, scale_parameter, epi_params.I_init, 1);

index(1:round(epi_params.I_init), 1) = (1:round(epi_params.I_init))';
curr_cum_infected = round(epi_params.I_init);

curr_cum_recovered = 0;

birthsEntire = [];
deathsEntire = [];
parentEntire = [];
indexEntire = [];

n_steps = length(new_infections);

for cntr = 1:n_steps

    [cntr n_steps]
    
    if new_infections(cntr) > 0
        new_births(1:new_infections(cntr), 1) = timeSeries.t_vals(cntr);
        new_deaths(1:new_infections(cntr), 1) = NaN;
        % AGH!! a random person was the parent!!
        %new_parent(1:new_infections(cntr), 1) = index(randi(length(births),new_infections(cntr),1));
        % NOW THE PARENT IS CHOSEN ACCORDING TO THEIR TRANSMISSION RATE!
        betas_norm = beta_val/sum(beta_val);
 %       new_infections(cntr)
        
        n_newInfections = (mnrnd(new_infections(cntr), betas_norm))';
        locs_pos = find(n_newInfections > 0)'; % should always be a row vector!!
        infection_cntr = 1;
        for loc = locs_pos
            newParentIndices(infection_cntr:(infection_cntr + n_newInfections(loc) - 1),1) = loc;
            infection_cntr = infection_cntr + n_newInfections(loc);
        end
%         newParentIndices
%         beta_val(newParentIndices)
%         mean(beta_val)
%         pause
%        new_infections(cntr)
        new_parent(1:new_infections(cntr), 1) = index(newParentIndices,1);
        newParentIndices = [];
        % end of addition
        new_index(1:new_infections(cntr), 1) = curr_cum_infected + (1:new_infections(cntr))';
        % giving new individuals their own transmission rates.
        new_beta_val(1:new_infections(cntr), 1) = gamrnd(epi_params.kbeta, scale_parameter, new_infections(cntr), 1);
        curr_cum_infected = curr_cum_infected + new_infections(cntr);
    else
        new_births = [];
        new_deaths = [];
        new_parent = [];
        new_index = [];
    end
    
    if new_recoveries(cntr) > 0
        vals = randperm(length(births));
        recov_locs = vals(1:new_recoveries(cntr));
        
        birthsEntire = [birthsEntire; births(recov_locs, 1)]; births(recov_locs, :) = [];
        deaths(recov_locs, 1) = timeSeries.t_vals(cntr);
        deathsEntire = [deathsEntire; deaths(recov_locs, 1)]; deaths(recov_locs, :) = [];
        parentEntire = [parentEntire; parent(recov_locs, 1)]; parent(recov_locs, :) = [];
        indexEntire = [indexEntire; index(recov_locs, 1)]; index(recov_locs, :) = [];
        beta_val(recov_locs, :) = [];

        curr_cum_recovered = curr_cum_recovered + new_recoveries(cntr);
        
    end
        
    births = [births; new_births]; new_births = [];
    deaths = [deaths; new_deaths]; new_deaths = [];
    parent = [parent; new_parent]; new_parent = [];
    beta_val = [beta_val; new_beta_val]; new_beta_val = [];
    index = [index; new_index]; new_index = [];

end

birthsEntire = [birthsEntire; births]; 
deaths(:, 1) = max(timeSeries.t_vals);
deathsEntire = [deathsEntire; deaths]; 
parentEntire = [parentEntire; parent]; 
indexEntire = [indexEntire; index];

[sortedIndexVals, IX] = sort(indexEntire);
sortedBirths = birthsEntire(IX);
sortedDeaths = deathsEntire(IX);
sortedParents = parentEntire(IX);
