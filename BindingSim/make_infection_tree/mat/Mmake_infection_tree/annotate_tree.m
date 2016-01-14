function [tr taxa] = annotate_tree(file, tipno)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



%TF = strcmp('ABG88729',headercell)


% function [tree,total,branches,names] = read_newick(tstr)
% Reads tree from standart file in Newick format 
% 
% Description of this format can be found at
% http://evolution.genetics.washington.edu/phylip/newicktree.html
%
% n -- the number of OTUs
% tree -- a string with the tree
% branches -- branch length
% names -- OTU names
% If no branch length provided, thay are set to ones.
%
%
% PHYLLAB toolbox v1.1
% ex. [tree n b nm, taxa] = read_newick_seq('hm_h3n2_mcc.nx.trees')

% opening file

ff = fopen(file,'r');
if ff == -1 
   n=-1;	tree = [];	
   disp('File open error in read_newick.');
   return;
end;

% The tree can be recorded in multiple lines. Removing LF, CR, etc...
[tr, k] = fscanf(ff,'%c'); 
%If total node# is 300
if tipno==300
tr1 = regexprep(tr,'300[&','xxx[&');
tr2 = regexprep(tr1,'[345]\d\d[&','[&');
tr = regexprep(tr2,'xxx','300');
end
%If total node# is 400
if tipno==400
tr1 = regexprep(tr,'400[&','xxx[&');
tr2 = regexprep(tr1,'[4567]\d\d[&','[&');
tr = regexprep(tr2,'xxx','400');
end

k = length(tr); %update the length
n=0;
i=1;
while i<k 
    if tr(i) <= ' '
       tr(i) =[]; 
       k=k-1;
    else
       i=i+1;
    end
    if tr(i) == ')'  
        n=n+1; 
    end
end
n=n+1; %  TEMP!!!!!!!!!!!!!!!!!!!!!!!!!

% Now replacing all symbolic names to numbers!
% Even if the name is already number it will be replaced
posO = 1; % position of last '(' - to determin the start of the nest 
flagO=0; % the status of previous symbol: 1 - if punctuation, 2 - other
flagA=0; % the status of previous symbol: 1 - if annotation, 0 - normal
nm = char(n,n);
taxa = struct;
cur_tip=1;
cur_inter_node=n+1;

[l,k] = size(tr);
i=1;
tr
%  position of '[' and ']'.
pos_annotation_s = 0;
pos_annotation_e = 0;
annotation = '';
while i<k
    if i>length(tr)
      break;
    end
    %skip the processes if tr(i) is in the bracket [...] 
    if (tr(i) == '[' & flagA ==0)
       pos_annotation_s = i;
       flagA = 1;
       i=i+1;
       continue;
    end
    if i == 4177
     disp test;
    end
    if (tr(i) == ']' & flagA ==1)
       pos_annotation_e = i;
       flagA = 0;
       i=i+1;
       continue;
    end
    %if (pos_annotation_s > 1 && flagA == 0)
      %annotation = tr(pos_annotation_s:pos_annotation_e);
    %end
    if (flagA == 0)
    %parse internal nodes
    if (tr(i) == ':' | tr(i) == ';' ) & flagO == 1
        taxa(cur_inter_node).id = cur_inter_node;
        taxa(cur_inter_node).name = '';
        if(pos_annotation_s == 0 | pos_annotation_e == 0)
            taxa(cur_inter_node).annotation = '';
        else
            taxa(cur_inter_node).annotation = tr(pos_annotation_s:pos_annotation_e);
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Need to find the netcharge given internal ID. Jul 2, 2013

        cur_inter_str = '';
        cur_inter_str = [num2str(cur_inter_node)];
        tr = sprintf('%s%s%s',tr(1:posO),cur_inter_str,tr(i:end)); 
        %tr = sprintf('%s%d%s',tr(1:posO),[],tr(i:end))
        [j,i]=size(sprintf('%s%s',tr(1:posO),cur_inter_str));
        %[j,i]=size(sprintf('%s%d',tr(1:posO),cur_inter_node));
        i=i+1; %shift the char index
        cur_inter_node=cur_inter_node+1;
    end
    
	%if (tr(i) == ')' | tr(i) == '(' | tr(i) == ',' | tr(i) == ':' ) & flagO == 2
    if (tr(i) == ')' | tr(i) == '(' | tr(i) == ',' | tr(i) == ':') & flagO == 2
      taxa(cur_tip).id = cur_tip;
      if pos_annotation_s~=0
      taxa(cur_tip).name = tr(posO+1:pos_annotation_s-1);
      taxa(cur_tip).annotation = tr(pos_annotation_s:pos_annotation_e);
      else
      taxa(cur_tip).name = tr(posO+1:i-1);
      end
      if cur_tip > 1 
          nm=strvcat(nm,tr(posO+1:i-1));
          else nm=tr(posO+1:i-1);
      end;
  
% Replace with the original tip ID in nexus file. The original tip ID can be found in nm.
% Search the netcharge for cur_tip.
      cur_tip %disp       
      pos_e = regexp(nm(end,1:end), '\[')-1; %skip the previous annotation[...]
      if isempty(pos_e) 
        pos_e = length(taxa(cur_tip).name); %keep same tip name
      end
      cur_tip_str = '';
      % add new annotation
      cur_tip_str = [num2str(nm(end,1:pos_e))];
      tr = sprintf('%s%s%s',tr(1:posO),cur_tip_str,tr(i:end));
      [j,i]=size(sprintf('%s%s',tr(1:posO),cur_tip_str));
      %tr = sprintf('%s%d%s',tr(1:posO),cur_tip,tr(i:end))
      %[j,i]=size(sprintf('%s%d',tr(1:posO),cur_tip));
      i=i+1;
      [j,k]=size(tr);
      cur_tip=cur_tip+1;
      flagO = 3;%%% added by HY
   end;
   if (tr(i) == '(' | tr(i) == ')' | tr(i) == ',')
        flagO = 1;posO = i;
   end;
   if tr(i) == ':' 
      flagO = 3;
   end;
   if (tr(i) ~= ')' & tr(i) ~= '(' & tr(i) ~= ',' & tr(i) ~= ';' & tr(i) ~= ':' & flagO ~= 3)
      flagO = 2;
   end;
   end
	i=i+1;  
end;


%***************
  % disp(tr);
%***************

% Save trees file with viral trait
outfile = [file(1:regexp(file, '\.')-1) '.traits' file(regexp(file, '\.'):end)];
fileID = fopen(outfile,'w');
fprintf(fileID,'%s',tr);
fclose(fileID);

end