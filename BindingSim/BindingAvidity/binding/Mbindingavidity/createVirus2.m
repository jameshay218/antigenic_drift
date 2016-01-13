%Create a virus particle with birth time, death time and source virus.
%Birth time: the time when individual is infected
%Death time
%   Detailed explanation goes here
%Mar 21, 2013

function [ Virus ] = createVirus2( vid, birth, death, sourceVirus, infectionK, beta, initialV, currentV )

field1 = 'vid';
value1 = num2cell(vid);
field2 = 'birth';
value2 = birth;
field3 = 'death';
value3 = death;
field4 = 'parent';
value4 = sourceVirus;
field5 = 'infectionK';
value5 = num2cell(infectionK);
field6 = 'beta';
value6 = beta;
field7 = 'initialV';
value7 = initialV;
field8 = 'currentV';
value8 = currentV;

Virus = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7,field8,value8);
end

