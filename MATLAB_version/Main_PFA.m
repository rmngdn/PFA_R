% I made some modifications on the initial version available at: 
% http://sysbio.sibcb.ac.cn/cb/chenlab/images/PFApackage_0.1.rar



%% LOAD YOUR DATA:

% %% TEST
% % File names:
% 
% gene_fileName = 'TestData/GLIO_Gene_Expression.txt';
% methy_fileName = 'TestData/GLIO_Methy_Expression.txt';
% mirna_fileName = 'TestData/GLIO_Mirna_Expression.txt';
% 
% 
% 
% 
% Output file name
res_fileName = 'Y_asthmadata.csv';
% 
% 
% % Have to adapt the importation method to your data format
% 
% data_gene = readtable(gene_fileName,'ReadRowNames', 1, 'HeaderLines', 0);
% data_gene = table2array(data_gene(:,1:215)); %last column is an artefact of importation
% 
% 
% data_Methy = readtable(methy_fileName, 'ReadRowNames', 1, 'HeaderLines', 0); 
% data_Methy = table2array(data_Methy(:,1:215)); %last column is an artefact of importation
% 
% data_Mirna = readtable(mirna_fileName, 'ReadRowNames', 1, 'HeaderLines', 0); 
% data_Mirna = table2array(data_Mirna(:,1:215)); %last column is an artefact of importation

%% ASTHMA
proteo = table2array(readtable('TestData/proteomics.csv','ReadRowNames', 0, 'HeaderLines', 0));
transc = table2array(readtable('TestData/transcriptomics.csv','ReadRowNames', 0, 'HeaderLines', 0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create dataList

% characteristics of the data

sample_num = 91;
k = 2;

dataList = cell(1,k);
dataList{1}=proteo;
dataList{2}=transc;




% Centralization of the data if needed
for i = 1:k
    dataList{i} = dataList{i} - (1/sample_num)*dataList{i}*ones(sample_num);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% capture the local sample-spectrum for each biological data type by Algorithm_1


%Threshold for accuracy
threshold=0.8; 

% init dmin
d_num = inf;

xList = cell(1,k);

for i = 1:k
    [u_i, eig_v_i] = Algorithm_1(dataList{i},sample_num);
    
    d_i = 1;
    for d_i=1:sample_num
        if sum(eig_v_i(1:d_i))/sum(eig_v_i)>threshold
            break;
        end
    end
    % re-compute Algo 1 with the number d_i:
    [u_i, eig_v_i] = Algorithm_1(dataList{i},d_i);
    
    % update dmin:
    if d_i < d_num
        d_num = d_i;
    end
    % compute x_i 
    xList{i} = u_i'*dataList{i}; %% the local sample-spectrum
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% capture the global sample-spectrum according to Algorithm_4


iter_num =1000;
lam_1 = 1;

[Y,w,L_list] = Algorithm_4(xList, sample_num, iter_num,  lam_1, d_num, k);%% Y is the global sample-spectrum


csvwrite(res_fileName,Y);