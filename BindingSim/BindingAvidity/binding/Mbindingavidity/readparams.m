function  [pa] = readparams(filename)
%UNTITLED Summary of this function goes here
%Read the parameter files

 if exist(filename)
    %pa_array = textread(filename,'%s')
    fid = fopen(filename)
    p = textscan(fid,'%s');
    pa_array = p{1};
 elseif exist(['../' filename])
    pa_array = textread(['../' filename],'%s');
 end
    pa = [];
 for i=1:length(pa_array)
   if strcmp(pa_array(i),'steps')
       pa.steps = str2num(char(pa_array(i+1)));
   end
      if strcmp(pa_array(i),'N')
       pa.N = str2num(char(pa_array(i+1)));
   end
   if strcmp(pa_array(i),'N_Strains')
       pa.N_Strains = str2num(char(pa_array(i+1)));
   end
   if strcmp(pa_array(i),'mutation')
       pa.mutation = str2num(char(pa_array(i+1)));
   end
   if strcmp(pa_array(i),'b_mutation')
       pa.b_mutation = str2num(char(pa_array(i+1)));
   end
   if strcmp(pa_array(i),'N_Infect')
       pa.N_Infect = str2num(char(pa_array(i+1)));
   end
   if strcmp(pa_array(i),'N_Binding')
       pa.N_Binding = str2num(char(pa_array(i+1)));
   end
   if strcmp(pa_array(i),'immunity')
       pa.immunity = str2num(char(pa_array(i+1)));
   end
   if strcmp(pa_array(i),'Vg')
       pa.Vg = str2num(char(pa_array(i+1)));
   end
   if strcmp(pa_array(i),'wan')
       pa.wan = str2num(char(pa_array(i+1)));
   end
   if strcmp(pa_array(i),'mu')
       pa.mu = str2num(char(pa_array(i+1)));
   end
   if strcmp(pa_array(i),'gamma')
       pa.gamma = str2num(char(pa_array(i+1)));
   end
   if strcmp(pa_array(i),'p')
       pa.p = str2num(char(pa_array(i+1)));
   end
   if strcmp(pa_array(i),'r')
       pa.r = str2num(char(pa_array(i+1)));
   end
   if strcmp(pa_array(i),'V0')
       pa.v0 = str2num(char(pa_array(i+1)));
   end
   if strcmp(pa_array(i),'a')
       pa.a = str2num(char(pa_array(i+1)));
   end
   if strcmp(pa_array(i),'b')
       pa.b = str2num(char(pa_array(i+1)));
   end
   if strcmp(pa_array(i),'c')
       pa.c = str2num(char(pa_array(i+1)));
   end
   if strcmp(pa_array(i),'kc')
       pa.kc = str2num(char(pa_array(i+1)));
   end
   if strcmp(pa_array(i),'ks')
       pa.ks = str2num(char(pa_array(i+1)));
   end
   if strcmp(pa_array(i),'nv')
       pa.nv = str2num(char(pa_array(i+1)));
   end 
   if strcmp(pa_array(i),'proj')
       pa.proj = char(pa_array(i+1));
   end
   if strcmp(pa_array(i),'AbB')
       pa.AbB = str2num(char(pa_array(i+1)));
   end
end
end

