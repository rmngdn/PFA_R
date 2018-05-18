function E_list = compute_err( w, sample_num, Y, L_list, xList, sList, k )

%init
for i = 1:k
    E_list{i} = [];
end

for i=1:sample_num
    for j = 1:k
        x_j = xList{j};
        E_list{j} = [E_list{j}; (1/sList{j})*w((j-1)*sample_num+i,1)*(Y(:,i)-L_list{j}*x_j(:,i))'*(Y(:,i)-L_list{j}*x_j(:,i))];
    end
end

end

