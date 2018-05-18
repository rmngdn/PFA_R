function w_1= Algorithm_2( sample_num, w, E_list, lam_1, k)
%% fix Y, solve W based Algorithm 2
E_or = [];
for i =1:k
    E_or = [E_or;E_list{i}];
end
E = (E_or)/sum(E_or);
[E_add,index]=sort(E);

p=inf;

for i=k*sample_num:-1:1
    o=(2*lam_1+sum(E_add(1:i,1)))/i - E_add(i);
    
    if o>=0
        p=i;
        break
    end
end

o=(2*lam_1+sum(E_add(1:p,1)))/p;
w(1:p,1)=(o-E_add(1:p,1))/(2*lam_1);
w((p+1):(k*sample_num))=0;

w_1 = -100*zeros(size(w,1),1);

w_1(index)=w;

end

