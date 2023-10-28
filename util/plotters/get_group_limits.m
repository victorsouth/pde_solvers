%%
function [ymin, ymax] = get_group_limits(group_names, lims)

result = lims(group_names{1});

for i = 2:length(group_names)
    ylimits = lims(group_names{i});
    
    if ylimits(1) < result(1)
        result(1) = ylimits(1);
    end
    if ylimits(2) > result(2)
        result(2) = ylimits(2);
    end
end

ymin = result(1);
ymax = result(2);

end