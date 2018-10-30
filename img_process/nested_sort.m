function sort_data = nested_sort(sort_mat, data)


s1 = sort_mat(:,1);

[~, ix] = sort(s1);

sort_data = data(ix);

if size(sort_mat,2)==2
    s2 = sort_mat(:,2);
    sorted_s2 = sort(s2);
    for n=1:numel(sorted_s2)
        sort_data(sorted_s2(n)==s2)= nested_sort(sort_mat(:,2:end),...
            sort_data(sorted_s2(n)==s2));
    end
    
end

end