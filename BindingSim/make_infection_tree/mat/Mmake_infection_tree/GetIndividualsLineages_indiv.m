function parentLineages = GetIndividualsLineages_indiv(n_seqs, indiv_sampled, parent)
    parentLineages = [];
for i = 1:n_seqs
    parentLineages(i).lin = GetSampleLineage(indiv_sampled(i), parent);
    parentLineages(i).tip = indiv_sampled(i);
end