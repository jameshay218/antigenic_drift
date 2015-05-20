function parentLineages = GetIndividualsLineages_indiv(n_seqs, indiv_sampled, parent)

for i = 1:n_seqs
    parentLineages(i).lin = GetSampleLineage(indiv_sampled(i), parent);
end