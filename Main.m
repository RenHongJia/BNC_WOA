%load data
data = xlsread('data.xlsx'); 
data =data';
classify_node_num = size(data,1);
%param set
nPop =20; 
MaxIt=100; 
param_set.max_fan_in = 3;

[zol,correct_rate,precision,recall,F1] = learn_struct_and_classified_bnc_swarm(data,classify_node_num,nPop,MaxIt,param_set);

fprintf("zol:%.4f\n correct_rate:%.4f\n precision:%.4f\n recall:%.4f\n F1:%.4f\n",zol,correct_rate,precision,recall,F1)