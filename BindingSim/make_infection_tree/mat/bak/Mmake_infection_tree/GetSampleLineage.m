function this_lineage = GetSampleLineage(focal_sample, parent)

this_lineage = focal_sample;
cntr = 1;

while 1 
%     if this_lineage(cntr) == 0
%         return;
%     end
    lineage_parent = parent(this_lineage(cntr));
    this_lineage(cntr+1) = lineage_parent;
    if (isnan(lineage_parent)) % || (lineage_parent == 0)
        return;
    end
    cntr = cntr + 1;
end