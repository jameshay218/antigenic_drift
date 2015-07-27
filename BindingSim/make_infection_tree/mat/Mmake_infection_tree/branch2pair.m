function [ pair ] = branch2pair( treefilename, traitfilename, ancnode)
%UNTITLED Summary of branch2pair goes here
% The function convert branch structure into pair.
%   Detailed explanation goes here
% ex: branch2pair('20130718/indiv_treeData_300.mat','20130718/virus_traits.mat')
% treefilename contains branches. traitfilename contains binding avidity.

folder = char(regexp(treefilename,'(.*)/','match'));
dat = load(['dat/' treefilename]);
b = dat.b;
names = dat.names;
clear dat;

ID_num = [];
for i=1:length(names)
    tokens = regexp(char(names(i)),'_','split');
    %if str2num(tokens{1,2})>0
        ID_num(i,1) = str2num(tokens{1,2});
    %end
end
%ID_num = regexp(names,'_','split');


dat2 = load(['dat/' traitfilename]);
traits = dat2.dat_VirusesArray(:,[1 5 8]); %col1 vid, col8 final virus netcharge
%TF = zeros(length(traits),1);
%TF(ID_num) = 1;

charge_txt = [];
%traits(logical(TF))
for i=1:length(ID_num)
    i
    loc_v = find(traits(:,1)==ID_num(i));
    if ~isempty(loc_v)
        charge_txt(i,1) = traits(loc_v,3);
        charge_txt(i,2) = traits(loc_v,2);
    else
        charge_txt(i,1) = NaN;
        charge_txt(i,2) = NaN;
    end
end


n = length(b) + 1; % number of tips
pair = zeros(2.*n-1,1);
%b_list = reshape(b',numel(b),1);
for i=1:length(b)
  index_tip1 = b(i,1);
  index_tip2 = b(i,2);
  parent = n+i;
  pair(index_tip1) = parent;
  pair(index_tip2) = parent;
end

% Update index to be 0 when the it is not in the lineage
new_pair = pair;
%for i=1:length(pair)
%    if i==599
%      disp test;
%    end
% [anc]= getParentList(pair, i);
% define an internal node that is the root of a new lineage
% if ancnode <= 300
%    anc_stop = pair(ancnode);
% else
%    anc_stop = ancnode; 
% end
% if isempty(find(anc==ancnode))
%   new_pair(i) = 0;
% end
%end

disp 'Remove all the nodes not within the lineage';
save(['dat/' folder 'pair.mat'], 'new_pair', 'pair', 'charge_txt', 'ancnode');



