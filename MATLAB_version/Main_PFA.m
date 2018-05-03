%% set input data and output filename
gene_fileName = 'GLIO_Gene_Expression.txt';
methy_fileName = 'GLIO_Methy_Expression.txt';
mirna_fileName = 'GLIO_Mirna_Expression.txt';
test = 'data_mirna_expression.csv'
res_fileName = 'global_sample_spectrum.csv';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read input data for program
data_test = csvread(test,1,1)
data_gene = readtable(gene_fileName, 'ReadRowNames', 1, 'HeaderLines', 0); 
data_gene = table2array(data_gene(:,1:215));%last column is an artefact of importation
data_Methy = readtable(methy_fileName, 'ReadRowNames', 1, 'HeaderLines', 0); 
data_Methy = table2array(data_Methy(:,1:215));
data_Mirna = readtable(mirna_fileName, 'ReadRowNames', 1, 'HeaderLines', 0); 
data_Mirna = table2array(data_Mirna(:,1:215));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% capture the local sample-spectrum for each biological data type by Algorithm_1
data_1 = data_gene;
data_2 = data_Methy;
data_3 = data_Mirna;

sample_num = size(data_1,2); 

[u_1,eig_v_1]=Algorithm_1(data_1,sample_num);
[u_2,eig_v_2]=Algorithm_1(data_2,sample_num);
[u_3,eig_v_3]=Algorithm_1(data_3,sample_num);

d_1 =1;
for d_1=1:sample_num
    if sum(eig_v_1(1:d_1))/sum(eig_v_1)>0.8
        break;
    end
end

d_2 =1;
for d_2=1:sample_num
    if sum(eig_v_2(1:d_2))/sum(eig_v_2)>0.8
        break;
    end
end

d_3 =1;
for d_3=1:sample_num
    if sum(eig_v_3(1:d_3))/sum(eig_v_3)>0.8
        break;
    end
end
  
[u_1,eig_v_1]=Algorithm_1(data_1,d_1);
[u_2,eig_v_2]=Algorithm_1(data_2,d_2);
[u_3,eig_v_3]=Algorithm_1(data_3,d_3);

x_1 = u_1'*data_1; %% the local sample-spectrum
x_2 = u_2'*data_2; %% the local sample-spectrum
x_3 = u_3'*data_3; %% the local sample-spectrum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% capture the global sample-spectrum according to Algorithm_4
d_num = min([d_1,d_2,d_3]);

iter_num =1000;
lam_1 = 1;

[Y,w,L_1,L_2,L_3] = Algorithm_4(x_1, x_2, x_3, sample_num, iter_num,  lam_1, d_num);%% Y is the global sample-spectrum

csvwrite(res_fileName,Y);