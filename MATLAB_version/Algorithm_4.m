 function [Y,w, L_list] = Algorithm_4(xList, sample_num, iter_num, lam_1, d_num, k)
% x_1, x_2 and x_3 is the local sample-spectrum
% Y is the global sample-spectrum

w=ones(k*sample_num,1)./(k*sample_num);

Y_final = [];
w_final = [];

final_err= inf;

%% doing Algorithm 3
for iter=1:iter_num

wList = cell(1,k);
sList = cell(1,k);
M = zeros(sample_num);
for i = 1:k
    wList{i} = diag(sqrt(w(((i-1)*sample_num+1):i*sample_num,1)));
    sList{i} = sum(diag(((eye(sample_num)-(xList{i})\(xList{i})))*((eye(sample_num)-(xList{i})\(xList{i})))'));
    M = M + (1/sList{i})*(wList{i}*(eye(sample_num)-(xList{i}*wList{i})\(xList{i}*wList{i})))*(wList{i}*(eye(sample_num)-(xList{i}*wList{i})\(xList{i}*wList{i})))';
end


% compute the global sample-spectrum (Y) based on the eigenvalue decomposition
[Y,Eigen_Value_all]=Find_K_Min_Eigen(M,d_num+1);

Y=Y(:,2:(d_num+1))';

L_list = cell(1,k);
for i = 1:k
    L_list{i} = (Y*wList{i})/(xList{i}*wList{i}); 
end

err = sum( diag(Y*M*Y') );
if err-final_err<0
    final_err = err;
    Y_final = Y;
    w_final = w;
else
    break
        
end
%% obtain err
E_list = compute_err( w, sample_num, Y, L_list, xList, sList, k );

%% doing Algorithm 2, solve w

w_tmp= Algorithm_2( sample_num, w, E_list, lam_1, k);
w=w_tmp;
 


end

Y=Y_final;
w=w_final;


 end


