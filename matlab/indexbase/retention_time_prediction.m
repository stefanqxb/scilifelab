close all
clear all

tic;
% calculating features
[aa_index,aa] =xlsread('rt_index.xls');
[retention_time,peptide,orginal] = xlsread('retention_time_peptide.xlsx');
% [data,peptide_2,orginal_2] = xlsread('hydrophobicity.xlsx');
% length_of_peptide = data(:,1);
% hydrophobicity = data(:,2);

% split data reduction

% method 1 using  3 words as feature


data_set = extract_words(peptide);
% changable 
% [~,q] = size(data_set);
% data_set = data_set(1:round(q/2));
%
score_mat1 =make_score_mat(aa_index); % order of aa is different
score_mat2 = eye(20);
score_mat = score_mat1+ score_mat2;
counting_mat = extract_aal(peptide,data_set,score_mat,aa_index); % A R N D C Q E G H I L K M F P S T W Y V

% demension reduction

feature_mat = cell2mat(counting_mat);

row = size(feature_mat,1);

for i= 1:row
    feature_mat(i,:) = feature_mat(i,:)/sum(feature_mat(i,:)); 
end



% set up train test set using bow as features
ratio = 0.9;
train_row = round(ratio*row);

train_set = feature_mat(1:train_row,:);
train_targ = retention_time(1:train_row,:);
test_set = feature_mat(train_row+1:end,:);
test_targ = retention_time(train_row+1:row,:);


% SVR
model = svmtrain(train_targ,train_set,'-s 3 -t 2 -g 40 -p 0.1 -c 100 -h 0 '); 
[predicted_label, accuracy, decision_values] = svmpredict(test_targ, test_set, model);

figure;
plot(test_targ,predicted_label,'r*');
xlabel('observed retention time');
ylabel('predicted retention time');
title('predictor of retention time');

gamma = corrcoef(test_targ,predicted_label);
% 
% diff_mat = abs(predicted_label - test_targ);
% max_t = max(test_targ);
% min_t = min(test_targ);
% step=100;
% 
% figure(2);
% hist(diff_mat,step);
% my_hist = hist(diff_mat,step);
% time_interval = time_95_diff(my_hist,max_t,min_t,step);

toc;

















