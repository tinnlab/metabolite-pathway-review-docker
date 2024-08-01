%% iMSEA - a strategy to decipher drug interaction

clear
clc
tic
load('/code/Data.mat')                                          % Users can replace this with their own data
load('/code/MRAN.mat')
load('/code/Pathway.mat')
load('/code/Metabolite.mat')

AllData = sample;
num_meta = size(AllData, 2);
num_samples = size(AllData, 1);

grp_A = grp_OXA; % A - OXA treatment group
grp_B = grp_VC; % B - VC treatment group
grp_AB = grp_OXA_VC; % AB - OXA+VC combination treatment group

num_A = size(grp_OXA, 1);
num_B = size(grp_VC, 1);
num_AB = size(grp_OXA_VC, 1);
num_model = size(grp_model, 1);

% prepare the data for PLS-DA model construction
X_A = [grp_A ; grp_model];
X_B = [grp_B ; grp_model];
X_AB = [grp_AB ; grp_model];
yA_true = [ones(num_A, 1); zeros(num_model, 1)];
yB_true = [ones(num_B, 1); zeros(num_model, 1)];
yAB_true = [ones(num_AB, 1); zeros(num_model, 1)];

VIP_A = getVIP(X_A, yA_true, 4);
VIP_B = getVIP(X_B, yB_true, 2);
VIP_AB = getVIP(X_AB, yAB_true, 4);

VIP_A = VIP_A .* VIP_A / num_meta;
VIP_B = VIP_B .* VIP_B / num_meta;
VIP_AB = VIP_AB .* VIP_AB / num_meta;

% Calculte Transition Matrix
Net = matData.fullNetwork;
Net = Net';
AllNode = matData.allNode;
index_reaction = pathway.index_reaction;
num_pathway = size(index_reaction, 1);
% get the index of the detected metabolites in ALLNODES
seed_nodes = zeros(num_meta, 1);
for i = 1:num_meta
    for j = 1:size(AllNode)
        if(strcmp(meta_KEGGID{i}, AllNode{j}))
            seed_nodes(i) = j;
        end
    end
end
clear i j
index = find(index_seed~=0);
startNodes = index_seed(index);
name_metabolites = name_metabolites(index);
meta_KEGGID = meta_KEGGID(index);

num_allnode = size(AllNode, 1);
num_seeds = size(startNodes, 1);
num_pathways = size(index_reaction, 1);
matTransition = zeros(num_seeds, num_pathways);
afterwalk = zeros(num_seeds, num_allnode);
for i = 1:num_seeds
    afterwalk(i, :) = RWR_new(Net, startNodes(i), 1, 0.3);                % The authors called RWR_new without supported files or calling an external library properly
    for j = 1:num_pathways
        index_pathway =  index_reaction{j};
        matTransition(i, j) = pathActivity_new(index_pathway, afterwalk(i, :));
    end
end
clear i j

VIP_A = VIP_A(index);
VIP_B = VIP_B(index);
VIP_AB = VIP_AB(index);

pathScore_A = matTransition' * VIP_A;
pathScore_B  = matTransition' * VIP_B;
pathScore_AB = matTransition' * VIP_AB;

% Calculate CI based on PA

EA = pathScore_A;
EB = pathScore_B;
EAB = pathScore_AB;

CI = zeros(num_pathways, 1);
for i = 1:num_pathways
    ma = max(EA(i), EB(i));
    CI(i) = ma / EAB(i);
end

% Calculte the perturbed PA of all groups
% Perturb the label 1000 times and Get the VIP of different drug treatment group
num_pert = 1000;
VIP_A_perturb = zeros(num_pert, num_meta);
VIP_B_perturb = zeros(num_pert, num_meta);
VIP_AB_perturb = zeros(num_pert, num_meta);
for i = 1:num_pert
    % perturb the group information  
    index1_perturb = randperm(size(yA_true,1))';
    index2_perturb = randperm(size(yB_true,1))';
    index3_perturb = randperm(size(yAB_true,1))';
    yA_perturb = yA_true(index1_perturb);
    yB_perturb = yB_true(index2_perturb);
    yAB_perturb = yAB_true(index3_perturb);
    VIP_A_perturb(i,:) = getVIP(X_A, yA_perturb, 3);
    VIP_B_perturb(i,:) = getVIP(X_B, yB_perturb, 3);
    VIP_AB_perturb(i,:) = getVIP(X_AB, yAB_perturb, 3);
end
VIP_OXA_perturb = VIP_A_perturb;
VIP_VC_perturb = VIP_B_perturb;
VIP_NMN_SS_perturb = VIP_AB_perturb;

VIP_A_perturb = VIP_A_perturb .* VIP_A_perturb / num_meta;
VIP_B_perturb = VIP_B_perturb .* VIP_B_perturb / num_meta;
VIP_AB_perturb = VIP_AB_perturb .* VIP_AB_perturb / num_meta;
VIP_A_perturb = VIP_A_perturb(:, index);
VIP_B_perturb = VIP_B_perturb(:, index);
VIP_AB_perturb = VIP_AB_perturb(:, index);

pathScore_A_perturb = VIP_A_perturb * matTransition;
pathScore_B_perturb  = VIP_B_perturb * matTransition;
pathScore_AB_perturb = VIP_AB_perturb * matTransition;

% Calculate the p value of each pathway
num_larger1 = zeros(num_pathway, 1);
num_larger2 = zeros(num_pathway, 1);
num_larger3 = zeros(num_pathway, 1);
for i = 1: num_pathway
    for j = 1:num_pert
        if(pathScore_A_perturb(j,i) >= pathScore_A(i))
            num_larger1(i) = num_larger1(i) + 1;
        end
        if(pathScore_B_perturb(j,i) >= pathScore_B(i))
            num_larger2(i) = num_larger2(i) + 1;
        end
        if(pathScore_AB_perturb(j,i) >= pathScore_AB(i))
            num_larger3(i) = num_larger3(i) + 1;
        end
    end
end
P_A = num_larger1/num_pert;
P_B = num_larger2/num_pert;
P_AB = num_larger3/num_pert;
toc

% save('.\Result\PA.mat','CI', 'P_A', 'P_B', 'P_AB', 'pathScore_A','pathScore_B',  'pathScore_AB')




