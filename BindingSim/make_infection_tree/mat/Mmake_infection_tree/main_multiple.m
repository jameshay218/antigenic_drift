function [] = main_multiple(tipno)
addpath(genpath('../plot_trait_change'));
figure; hold on;
%proj = '20130812_med/';

%proj1 = '20130816_med_n8_c08/';
%proj1 = '20130816_med_n4_c09/';
%proj2 = '20130814_med_fixb/';

%proj1 = '20130818_lar_n4_c07/';
%proj2 = '20130819_lar_fixb28/';

proj1 = '20130827_lar_n4_c08_N5/samp400/';
proj2 = '20130827_lar_fixb_c05/';

%proj1 = '20130828_lar_n4_c08_kc01/';
%proj2 = '20130828_lar_fixb_c05_kc01/';

tips = [400 399 398 397 396];
plot_lineages(proj1, 'b', tipno, tips);
tips = [400 399 398 397 396];
plot_lineages(proj2, 'r', tipno, tips);

yrs_e = 45;
set(gca,'XTick',0:365*10:365*yrs_e);
set(gca,'XTickLabel',{'0','10','20','30','40','45'});

end



function [] = plot_lineages(proj, color, tipno, tips)
[tree taxa] = annotate_tree(['dat/' proj 'indiv_genealogy_' num2str(tipno) '.branch.tree'], tipno);
[t, n, b, nm] = read_newick_seq(['dat/' proj 'indiv_genealogy_' num2str(tipno) '.tree']);
v = zeros(length(b),1)
for i=1:length(taxa)
    if i<=tipno
        id = str2num(taxa(1,i).name);
    else
        id = i;
    end
    str = taxa(1,i).annotation;
    str1 = regexp(str, 'Bindingscr.set=\"\d.\d+\"', 'match')
    if ~isempty(char(str1))
        v(id,1)=str2num(cell2mat(regexp(char(str1), '\d.\d+', 'match')));
    else
        v(id,1)=0;
    end
end
pairs = string2pairs(['dat/' proj 'indiv_genealogy_' num2str(tipno) '.topo.tree'],tipno);
disp pairs

dist_root = [];
n_total = length(pairs);
for i=1:length(pairs)
    if i == tipno
      disp test
    end
    [lin]= getParentList(pairs, i);
    distance = 0;
    distance = sum(b(lin(1:end-1)));
    %for j=1:length(lin)
    %    if lin(j) == 0
    %        break;
    %    end
    %    distance = distance + b(lin(j)); 
    %end
    dist_root(i,:) = distance;
end

save(['dat/' proj 'indiv_genealogy_' num2str(tipno) '_taxa.mat'], 't', 'n', 'b', 'nm', 'taxa', 'pairs', 'dist_root', 'v');
disp distance
for i1=1:1
v(tipno*2-i1)=0.4;
end
%v(598)=0.4;
%v(599)=0.4;

%for tip=tipno:-1:tipno-4
for i=1:length(tips)
    tip = tips(i);
    [lin]= getParentList(pairs, tip);
    plot(dist_root(lin(1:end-1)),v(lin(1:end-1)),[color '.-']);
end
%plot(dist_root(1:tipno),v(1:tipno),[color 'x']);
end