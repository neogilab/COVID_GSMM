addpath(genpath('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/Test/Human-GEM/'))
addpath(genpath('/home/anoop/Downloads/RAVEN-2.4.0'))

initCobraToolbox(false) 
changeCobraSolver ('gurobi', 'all');

load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_1_HC.mat')


COVID19_1_HC_GEM = addReaction(COVID19_1_HC_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_1_HC_GEM, [], true, false, false, essentialTasks);


COVID19_1_HC_GEM = simplifyModel(COVID19_1_HC_GEM);
COVID19_1_HC_GEM.b = COVID19_1_HC_GEM.b(:,1);
COVID19_1_HC_GEM = setParam(COVID19_1_HC_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_1_HC.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_1_HC_GEM = setParam(COVID19_1_HC_GEM,'lb',cRxns,cLB);
COVID19_1_HC_GEM = setParam(COVID19_1_HC_GEM,'ub',cRxns,cUB);

sol_COVID19_1_HC_GEM = solveLP(COVID19_1_HC_GEM);
writematrix(sol_COVID19_1_HC_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_1_HC_GEM.txt')
writecell(COVID19_1_HC_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_1_HC_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_2_HC.mat')


COVID19_2_HC_GEM = addReaction(COVID19_2_HC_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_2_HC_GEM, [], true, false, false, essentialTasks);


COVID19_2_HC_GEM = simplifyModel(COVID19_2_HC_GEM);
COVID19_2_HC_GEM.b = COVID19_2_HC_GEM.b(:,1);
COVID19_2_HC_GEM = setParam(COVID19_2_HC_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_2_HC.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_2_HC_GEM = setParam(COVID19_2_HC_GEM,'lb',cRxns,cLB);
COVID19_2_HC_GEM = setParam(COVID19_2_HC_GEM,'ub',cRxns,cUB);

sol_COVID19_2_HC_GEM = solveLP(COVID19_2_HC_GEM);
writematrix(sol_COVID19_2_HC_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_2_HC_GEM.txt')
writecell(COVID19_2_HC_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_2_HC_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_4_HC.mat')


COVID19_4_HC_GEM = addReaction(COVID19_4_HC_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_4_HC_GEM, [], true, false, false, essentialTasks);


COVID19_4_HC_GEM = simplifyModel(COVID19_4_HC_GEM);
COVID19_4_HC_GEM.b = COVID19_4_HC_GEM.b(:,1);
COVID19_4_HC_GEM = setParam(COVID19_4_HC_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_4_HC.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_4_HC_GEM = setParam(COVID19_4_HC_GEM,'lb',cRxns,cLB);
COVID19_4_HC_GEM = setParam(COVID19_4_HC_GEM,'ub',cRxns,cUB);

sol_COVID19_4_HC_GEM = solveLP(COVID19_4_HC_GEM);
writematrix(sol_COVID19_4_HC_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_4_HC_GEM.txt')
writecell(COVID19_4_HC_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_4_HC_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_5_HC.mat')


COVID19_5_HC_GEM = addReaction(COVID19_5_HC_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_5_HC_GEM, [], true, false, false, essentialTasks);


COVID19_5_HC_GEM = simplifyModel(COVID19_5_HC_GEM);
COVID19_5_HC_GEM.b = COVID19_5_HC_GEM.b(:,1);
COVID19_5_HC_GEM = setParam(COVID19_5_HC_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_5_HC.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_5_HC_GEM = setParam(COVID19_5_HC_GEM,'lb',cRxns,cLB);
COVID19_5_HC_GEM = setParam(COVID19_5_HC_GEM,'ub',cRxns,cUB);

sol_COVID19_5_HC_GEM = solveLP(COVID19_5_HC_GEM);
writematrix(sol_COVID19_5_HC_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_5_HC_GEM.txt')
writecell(COVID19_5_HC_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_5_HC_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_6_HC.mat')


COVID19_6_HC_GEM = addReaction(COVID19_6_HC_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_6_HC_GEM, [], true, false, false, essentialTasks);


COVID19_6_HC_GEM = simplifyModel(COVID19_6_HC_GEM);
COVID19_6_HC_GEM.b = COVID19_6_HC_GEM.b(:,1);
COVID19_6_HC_GEM = setParam(COVID19_6_HC_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_6_HC.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_6_HC_GEM = setParam(COVID19_6_HC_GEM,'lb',cRxns,cLB);
COVID19_6_HC_GEM = setParam(COVID19_6_HC_GEM,'ub',cRxns,cUB);

sol_COVID19_6_HC_GEM = solveLP(COVID19_6_HC_GEM);
writematrix(sol_COVID19_6_HC_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_6_HC_GEM.txt')
writecell(COVID19_6_HC_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_6_HC_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_7_HC.mat')


COVID19_7_HC_GEM = addReaction(COVID19_7_HC_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_7_HC_GEM, [], true, false, false, essentialTasks);


COVID19_7_HC_GEM = simplifyModel(COVID19_7_HC_GEM);
COVID19_7_HC_GEM.b = COVID19_7_HC_GEM.b(:,1);
COVID19_7_HC_GEM = setParam(COVID19_7_HC_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_7_HC.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_7_HC_GEM = setParam(COVID19_7_HC_GEM,'lb',cRxns,cLB);
COVID19_7_HC_GEM = setParam(COVID19_7_HC_GEM,'ub',cRxns,cUB);

sol_COVID19_7_HC_GEM = solveLP(COVID19_7_HC_GEM);
writematrix(sol_COVID19_7_HC_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_7_HC_GEM.txt')
writecell(COVID19_7_HC_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_7_HC_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_8_HC.mat')


COVID19_8_HC_GEM = addReaction(COVID19_8_HC_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_8_HC_GEM, [], true, false, false, essentialTasks);


COVID19_8_HC_GEM = simplifyModel(COVID19_8_HC_GEM);
COVID19_8_HC_GEM.b = COVID19_8_HC_GEM.b(:,1);
COVID19_8_HC_GEM = setParam(COVID19_8_HC_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_8_HC.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_8_HC_GEM = setParam(COVID19_8_HC_GEM,'lb',cRxns,cLB);
COVID19_8_HC_GEM = setParam(COVID19_8_HC_GEM,'ub',cRxns,cUB);

sol_COVID19_8_HC_GEM = solveLP(COVID19_8_HC_GEM);
writematrix(sol_COVID19_8_HC_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_8_HC_GEM.txt')
writecell(COVID19_8_HC_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_8_HC_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_9_HC.mat')


COVID19_9_HC_GEM = addReaction(COVID19_9_HC_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_9_HC_GEM, [], true, false, false, essentialTasks);


COVID19_9_HC_GEM = simplifyModel(COVID19_9_HC_GEM);
COVID19_9_HC_GEM.b = COVID19_9_HC_GEM.b(:,1);
COVID19_9_HC_GEM = setParam(COVID19_9_HC_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_9_HC.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_9_HC_GEM = setParam(COVID19_9_HC_GEM,'lb',cRxns,cLB);
COVID19_9_HC_GEM = setParam(COVID19_9_HC_GEM,'ub',cRxns,cUB);

sol_COVID19_9_HC_GEM = solveLP(COVID19_9_HC_GEM);
writematrix(sol_COVID19_9_HC_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_9_HC_GEM.txt')
writecell(COVID19_9_HC_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_9_HC_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_11_HC.mat')


COVID19_11_HC_GEM = addReaction(COVID19_11_HC_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_11_HC_GEM, [], true, false, false, essentialTasks);


COVID19_11_HC_GEM = simplifyModel(COVID19_11_HC_GEM);
COVID19_11_HC_GEM.b = COVID19_11_HC_GEM.b(:,1);
COVID19_11_HC_GEM = setParam(COVID19_11_HC_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_11_HC.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_11_HC_GEM = setParam(COVID19_11_HC_GEM,'lb',cRxns,cLB);
COVID19_11_HC_GEM = setParam(COVID19_11_HC_GEM,'ub',cRxns,cUB);

sol_COVID19_11_HC_GEM = solveLP(COVID19_11_HC_GEM);
writematrix(sol_COVID19_11_HC_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_11_HC_GEM.txt')
writecell(COVID19_11_HC_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_11_HC_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_16_HC.mat')


COVID19_16_HC_GEM = addReaction(COVID19_16_HC_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_16_HC_GEM, [], true, false, false, essentialTasks);


COVID19_16_HC_GEM = simplifyModel(COVID19_16_HC_GEM);
COVID19_16_HC_GEM.b = COVID19_16_HC_GEM.b(:,1);
COVID19_16_HC_GEM = setParam(COVID19_16_HC_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_16_HC.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_16_HC_GEM = setParam(COVID19_16_HC_GEM,'lb',cRxns,cLB);
COVID19_16_HC_GEM = setParam(COVID19_16_HC_GEM,'ub',cRxns,cUB);

sol_COVID19_16_HC_GEM = solveLP(COVID19_16_HC_GEM);
writematrix(sol_COVID19_16_HC_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_16_HC_GEM.txt')
writecell(COVID19_16_HC_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_16_HC_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_17_HC.mat')


COVID19_17_HC_GEM = addReaction(COVID19_17_HC_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_17_HC_GEM, [], true, false, false, essentialTasks);


COVID19_17_HC_GEM = simplifyModel(COVID19_17_HC_GEM);
COVID19_17_HC_GEM.b = COVID19_17_HC_GEM.b(:,1);
COVID19_17_HC_GEM = setParam(COVID19_17_HC_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_17_HC.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_17_HC_GEM = setParam(COVID19_17_HC_GEM,'lb',cRxns,cLB);
COVID19_17_HC_GEM = setParam(COVID19_17_HC_GEM,'ub',cRxns,cUB);

sol_COVID19_17_HC_GEM = solveLP(COVID19_17_HC_GEM);
writematrix(sol_COVID19_17_HC_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_17_HC_GEM.txt')
writecell(COVID19_17_HC_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_17_HC_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_19_HC.mat')


COVID19_19_HC_GEM = addReaction(COVID19_19_HC_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_19_HC_GEM, [], true, false, false, essentialTasks);


COVID19_19_HC_GEM = simplifyModel(COVID19_19_HC_GEM);
COVID19_19_HC_GEM.b = COVID19_19_HC_GEM.b(:,1);
COVID19_19_HC_GEM = setParam(COVID19_19_HC_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_19_HC.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_19_HC_GEM = setParam(COVID19_19_HC_GEM,'lb',cRxns,cLB);
COVID19_19_HC_GEM = setParam(COVID19_19_HC_GEM,'ub',cRxns,cUB);

sol_COVID19_19_HC_GEM = solveLP(COVID19_19_HC_GEM);
writematrix(sol_COVID19_19_HC_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_19_HC_GEM.txt')
writecell(COVID19_19_HC_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_19_HC_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_20_HC.mat')


COVID19_20_HC_GEM = addReaction(COVID19_20_HC_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_20_HC_GEM, [], true, false, false, essentialTasks);


COVID19_20_HC_GEM = simplifyModel(COVID19_20_HC_GEM);
COVID19_20_HC_GEM.b = COVID19_20_HC_GEM.b(:,1);
COVID19_20_HC_GEM = setParam(COVID19_20_HC_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_20_HC.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_20_HC_GEM = setParam(COVID19_20_HC_GEM,'lb',cRxns,cLB);
COVID19_20_HC_GEM = setParam(COVID19_20_HC_GEM,'ub',cRxns,cUB);

sol_COVID19_20_HC_GEM = solveLP(COVID19_20_HC_GEM);
writematrix(sol_COVID19_20_HC_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_20_HC_GEM.txt')
writecell(COVID19_20_HC_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_20_HC_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_21_HC.mat')


COVID19_21_HC_GEM = addReaction(COVID19_21_HC_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_21_HC_GEM, [], true, false, false, essentialTasks);


COVID19_21_HC_GEM = simplifyModel(COVID19_21_HC_GEM);
COVID19_21_HC_GEM.b = COVID19_21_HC_GEM.b(:,1);
COVID19_21_HC_GEM = setParam(COVID19_21_HC_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_21_HC.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_21_HC_GEM = setParam(COVID19_21_HC_GEM,'lb',cRxns,cLB);
COVID19_21_HC_GEM = setParam(COVID19_21_HC_GEM,'ub',cRxns,cUB);

sol_COVID19_21_HC_GEM = solveLP(COVID19_21_HC_GEM);
writematrix(sol_COVID19_21_HC_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_21_HC_GEM.txt')
writecell(COVID19_21_HC_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_21_HC_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_22_HC.mat')


COVID19_22_HC_GEM = addReaction(COVID19_22_HC_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_22_HC_GEM, [], true, false, false, essentialTasks);


COVID19_22_HC_GEM = simplifyModel(COVID19_22_HC_GEM);
COVID19_22_HC_GEM.b = COVID19_22_HC_GEM.b(:,1);
COVID19_22_HC_GEM = setParam(COVID19_22_HC_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_22_HC.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_22_HC_GEM = setParam(COVID19_22_HC_GEM,'lb',cRxns,cLB);
COVID19_22_HC_GEM = setParam(COVID19_22_HC_GEM,'ub',cRxns,cUB);

sol_COVID19_22_HC_GEM = solveLP(COVID19_22_HC_GEM);
writematrix(sol_COVID19_22_HC_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_22_HC_GEM.txt')
writecell(COVID19_22_HC_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_22_HC_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_25_HC.mat')


COVID19_25_HC_GEM = addReaction(COVID19_25_HC_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_25_HC_GEM, [], true, false, false, essentialTasks);


COVID19_25_HC_GEM = simplifyModel(COVID19_25_HC_GEM);
COVID19_25_HC_GEM.b = COVID19_25_HC_GEM.b(:,1);
COVID19_25_HC_GEM = setParam(COVID19_25_HC_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_25_HC.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_25_HC_GEM = setParam(COVID19_25_HC_GEM,'lb',cRxns,cLB);
COVID19_25_HC_GEM = setParam(COVID19_25_HC_GEM,'ub',cRxns,cUB);

sol_COVID19_25_HC_GEM = solveLP(COVID19_25_HC_GEM);
writematrix(sol_COVID19_25_HC_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_25_HC_GEM.txt')
writecell(COVID19_25_HC_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_25_HC_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_27_HC.mat')


COVID19_27_HC_GEM = addReaction(COVID19_27_HC_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_27_HC_GEM, [], true, false, false, essentialTasks);


COVID19_27_HC_GEM = simplifyModel(COVID19_27_HC_GEM);
COVID19_27_HC_GEM.b = COVID19_27_HC_GEM.b(:,1);
COVID19_27_HC_GEM = setParam(COVID19_27_HC_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_27_HC.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_27_HC_GEM = setParam(COVID19_27_HC_GEM,'lb',cRxns,cLB);
COVID19_27_HC_GEM = setParam(COVID19_27_HC_GEM,'ub',cRxns,cUB);

sol_COVID19_27_HC_GEM = solveLP(COVID19_27_HC_GEM);
writematrix(sol_COVID19_27_HC_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_27_HC_GEM.txt')
writecell(COVID19_27_HC_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_27_HC_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_28_HC.mat')


COVID19_28_HC_GEM = addReaction(COVID19_28_HC_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_28_HC_GEM, [], true, false, false, essentialTasks);


COVID19_28_HC_GEM = simplifyModel(COVID19_28_HC_GEM);
COVID19_28_HC_GEM.b = COVID19_28_HC_GEM.b(:,1);
COVID19_28_HC_GEM = setParam(COVID19_28_HC_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_28_HC.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_28_HC_GEM = setParam(COVID19_28_HC_GEM,'lb',cRxns,cLB);
COVID19_28_HC_GEM = setParam(COVID19_28_HC_GEM,'ub',cRxns,cUB);

sol_COVID19_28_HC_GEM = solveLP(COVID19_28_HC_GEM);
writematrix(sol_COVID19_28_HC_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_28_HC_GEM.txt')
writecell(COVID19_28_HC_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_28_HC_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_29_HC.mat')


COVID19_29_HC_GEM = addReaction(COVID19_29_HC_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_29_HC_GEM, [], true, false, false, essentialTasks);


COVID19_29_HC_GEM = simplifyModel(COVID19_29_HC_GEM);
COVID19_29_HC_GEM.b = COVID19_29_HC_GEM.b(:,1);
COVID19_29_HC_GEM = setParam(COVID19_29_HC_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_29_HC.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_29_HC_GEM = setParam(COVID19_29_HC_GEM,'lb',cRxns,cLB);
COVID19_29_HC_GEM = setParam(COVID19_29_HC_GEM,'ub',cRxns,cUB);

sol_COVID19_29_HC_GEM = solveLP(COVID19_29_HC_GEM);
writematrix(sol_COVID19_29_HC_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_29_HC_GEM.txt')
writecell(COVID19_29_HC_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_29_HC_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_24_HC.mat')


COVID19_24_HC_GEM = addReaction(COVID19_24_HC_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_24_HC_GEM, [], true, false, false, essentialTasks);


COVID19_24_HC_GEM = simplifyModel(COVID19_24_HC_GEM);
COVID19_24_HC_GEM.b = COVID19_24_HC_GEM.b(:,1);
COVID19_24_HC_GEM = setParam(COVID19_24_HC_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_24_HC.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_24_HC_GEM = setParam(COVID19_24_HC_GEM,'lb',cRxns,cLB);
COVID19_24_HC_GEM = setParam(COVID19_24_HC_GEM,'ub',cRxns,cUB);

sol_COVID19_24_HC_GEM = solveLP(COVID19_24_HC_GEM);
writematrix(sol_COVID19_24_HC_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_24_HC_GEM.txt')
writecell(COVID19_24_HC_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_24_HC_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_31_HC.mat')


COVID19_31_HC_GEM = addReaction(COVID19_31_HC_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_31_HC_GEM, [], true, false, false, essentialTasks);


COVID19_31_HC_GEM = simplifyModel(COVID19_31_HC_GEM);
COVID19_31_HC_GEM.b = COVID19_31_HC_GEM.b(:,1);
COVID19_31_HC_GEM = setParam(COVID19_31_HC_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_31_HC.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_31_HC_GEM = setParam(COVID19_31_HC_GEM,'lb',cRxns,cLB);
COVID19_31_HC_GEM = setParam(COVID19_31_HC_GEM,'ub',cRxns,cUB);

sol_COVID19_31_HC_GEM = solveLP(COVID19_31_HC_GEM);
writematrix(sol_COVID19_31_HC_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_31_HC_GEM.txt')
writecell(COVID19_31_HC_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_31_HC_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_3_Conv.mat')


COVID19_3_Conv_GEM = addReaction(COVID19_3_Conv_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_3_Conv_GEM, [], true, false, false, essentialTasks);


COVID19_3_Conv_GEM = simplifyModel(COVID19_3_Conv_GEM);
COVID19_3_Conv_GEM.b = COVID19_3_Conv_GEM.b(:,1);
COVID19_3_Conv_GEM = setParam(COVID19_3_Conv_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_3_Conv.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_3_Conv_GEM = setParam(COVID19_3_Conv_GEM,'lb',cRxns,cLB);
COVID19_3_Conv_GEM = setParam(COVID19_3_Conv_GEM,'ub',cRxns,cUB);

sol_COVID19_3_Conv_GEM = solveLP(COVID19_3_Conv_GEM);
writematrix(sol_COVID19_3_Conv_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_3_Conv_GEM.txt')
writecell(COVID19_3_Conv_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_3_Conv_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_10_Conv.mat')


COVID19_10_Conv_GEM = addReaction(COVID19_10_Conv_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_10_Conv_GEM, [], true, false, false, essentialTasks);


COVID19_10_Conv_GEM = simplifyModel(COVID19_10_Conv_GEM);
COVID19_10_Conv_GEM.b = COVID19_10_Conv_GEM.b(:,1);
COVID19_10_Conv_GEM = setParam(COVID19_10_Conv_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_10_Conv.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_10_Conv_GEM = setParam(COVID19_10_Conv_GEM,'lb',cRxns,cLB);
COVID19_10_Conv_GEM = setParam(COVID19_10_Conv_GEM,'ub',cRxns,cUB);

sol_COVID19_10_Conv_GEM = solveLP(COVID19_10_Conv_GEM);
writematrix(sol_COVID19_10_Conv_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_10_Conv_GEM.txt')
writecell(COVID19_10_Conv_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_10_Conv_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_12_Conv.mat')


COVID19_12_Conv_GEM = addReaction(COVID19_12_Conv_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_12_Conv_GEM, [], true, false, false, essentialTasks);


COVID19_12_Conv_GEM = simplifyModel(COVID19_12_Conv_GEM);
COVID19_12_Conv_GEM.b = COVID19_12_Conv_GEM.b(:,1);
COVID19_12_Conv_GEM = setParam(COVID19_12_Conv_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_12_Conv.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_12_Conv_GEM = setParam(COVID19_12_Conv_GEM,'lb',cRxns,cLB);
COVID19_12_Conv_GEM = setParam(COVID19_12_Conv_GEM,'ub',cRxns,cUB);

sol_COVID19_12_Conv_GEM = solveLP(COVID19_12_Conv_GEM);
writematrix(sol_COVID19_12_Conv_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_12_Conv_GEM.txt')
writecell(COVID19_12_Conv_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_12_Conv_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_13_Conv.mat')


COVID19_13_Conv_GEM = addReaction(COVID19_13_Conv_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_13_Conv_GEM, [], true, false, false, essentialTasks);


COVID19_13_Conv_GEM = simplifyModel(COVID19_13_Conv_GEM);
COVID19_13_Conv_GEM.b = COVID19_13_Conv_GEM.b(:,1);
COVID19_13_Conv_GEM = setParam(COVID19_13_Conv_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_13_Conv.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_13_Conv_GEM = setParam(COVID19_13_Conv_GEM,'lb',cRxns,cLB);
COVID19_13_Conv_GEM = setParam(COVID19_13_Conv_GEM,'ub',cRxns,cUB);

sol_COVID19_13_Conv_GEM = solveLP(COVID19_13_Conv_GEM);
writematrix(sol_COVID19_13_Conv_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_13_Conv_GEM.txt')
writecell(COVID19_13_Conv_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_13_Conv_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_14_Conv.mat')


COVID19_14_Conv_GEM = addReaction(COVID19_14_Conv_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_14_Conv_GEM, [], true, false, false, essentialTasks);


COVID19_14_Conv_GEM = simplifyModel(COVID19_14_Conv_GEM);
COVID19_14_Conv_GEM.b = COVID19_14_Conv_GEM.b(:,1);
COVID19_14_Conv_GEM = setParam(COVID19_14_Conv_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_14_Conv.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_14_Conv_GEM = setParam(COVID19_14_Conv_GEM,'lb',cRxns,cLB);
COVID19_14_Conv_GEM = setParam(COVID19_14_Conv_GEM,'ub',cRxns,cUB);

sol_COVID19_14_Conv_GEM = solveLP(COVID19_14_Conv_GEM);
writematrix(sol_COVID19_14_Conv_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_14_Conv_GEM.txt')
writecell(COVID19_14_Conv_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_14_Conv_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_15_Conv.mat')


COVID19_15_Conv_GEM = addReaction(COVID19_15_Conv_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_15_Conv_GEM, [], true, false, false, essentialTasks);


COVID19_15_Conv_GEM = simplifyModel(COVID19_15_Conv_GEM);
COVID19_15_Conv_GEM.b = COVID19_15_Conv_GEM.b(:,1);
COVID19_15_Conv_GEM = setParam(COVID19_15_Conv_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_15_Conv.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_15_Conv_GEM = setParam(COVID19_15_Conv_GEM,'lb',cRxns,cLB);
COVID19_15_Conv_GEM = setParam(COVID19_15_Conv_GEM,'ub',cRxns,cUB);

sol_COVID19_15_Conv_GEM = solveLP(COVID19_15_Conv_GEM);
writematrix(sol_COVID19_15_Conv_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_15_Conv_GEM.txt')
writecell(COVID19_15_Conv_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_15_Conv_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_18_Conv.mat')


COVID19_18_Conv_GEM = addReaction(COVID19_18_Conv_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_18_Conv_GEM, [], true, false, false, essentialTasks);


COVID19_18_Conv_GEM = simplifyModel(COVID19_18_Conv_GEM);
COVID19_18_Conv_GEM.b = COVID19_18_Conv_GEM.b(:,1);
COVID19_18_Conv_GEM = setParam(COVID19_18_Conv_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_18_Conv.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_18_Conv_GEM = setParam(COVID19_18_Conv_GEM,'lb',cRxns,cLB);
COVID19_18_Conv_GEM = setParam(COVID19_18_Conv_GEM,'ub',cRxns,cUB);

sol_COVID19_18_Conv_GEM = solveLP(COVID19_18_Conv_GEM);
writematrix(sol_COVID19_18_Conv_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_18_Conv_GEM.txt')
writecell(COVID19_18_Conv_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_18_Conv_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_23_Conv.mat')


COVID19_23_Conv_GEM = addReaction(COVID19_23_Conv_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_23_Conv_GEM, [], true, false, false, essentialTasks);


COVID19_23_Conv_GEM = simplifyModel(COVID19_23_Conv_GEM);
COVID19_23_Conv_GEM.b = COVID19_23_Conv_GEM.b(:,1);
COVID19_23_Conv_GEM = setParam(COVID19_23_Conv_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_23_Conv.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_23_Conv_GEM = setParam(COVID19_23_Conv_GEM,'lb',cRxns,cLB);
COVID19_23_Conv_GEM = setParam(COVID19_23_Conv_GEM,'ub',cRxns,cUB);

sol_COVID19_23_Conv_GEM = solveLP(COVID19_23_Conv_GEM);
writematrix(sol_COVID19_23_Conv_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_23_Conv_GEM.txt')
writecell(COVID19_23_Conv_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_23_Conv_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_26_Conv.mat')


COVID19_26_Conv_GEM = addReaction(COVID19_26_Conv_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_26_Conv_GEM, [], true, false, false, essentialTasks);


COVID19_26_Conv_GEM = simplifyModel(COVID19_26_Conv_GEM);
COVID19_26_Conv_GEM.b = COVID19_26_Conv_GEM.b(:,1);
COVID19_26_Conv_GEM = setParam(COVID19_26_Conv_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_26_Conv.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_26_Conv_GEM = setParam(COVID19_26_Conv_GEM,'lb',cRxns,cLB);
COVID19_26_Conv_GEM = setParam(COVID19_26_Conv_GEM,'ub',cRxns,cUB);

sol_COVID19_26_Conv_GEM = solveLP(COVID19_26_Conv_GEM);
writematrix(sol_COVID19_26_Conv_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_26_Conv_GEM.txt')
writecell(COVID19_26_Conv_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_26_Conv_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_32_Conv.mat')


COVID19_32_Conv_GEM = addReaction(COVID19_32_Conv_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_32_Conv_GEM, [], true, false, false, essentialTasks);


COVID19_32_Conv_GEM = simplifyModel(COVID19_32_Conv_GEM);
COVID19_32_Conv_GEM.b = COVID19_32_Conv_GEM.b(:,1);
COVID19_32_Conv_GEM = setParam(COVID19_32_Conv_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_32_Conv.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_32_Conv_GEM = setParam(COVID19_32_Conv_GEM,'lb',cRxns,cLB);
COVID19_32_Conv_GEM = setParam(COVID19_32_Conv_GEM,'ub',cRxns,cUB);

sol_COVID19_32_Conv_GEM = solveLP(COVID19_32_Conv_GEM);
writematrix(sol_COVID19_32_Conv_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_32_Conv_GEM.txt')
writecell(COVID19_32_Conv_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_32_Conv_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_002_Mild.mat')


COVID19_002_Mild_GEM = addReaction(COVID19_002_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_002_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_002_Mild_GEM = simplifyModel(COVID19_002_Mild_GEM);
COVID19_002_Mild_GEM.b = COVID19_002_Mild_GEM.b(:,1);
COVID19_002_Mild_GEM = setParam(COVID19_002_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_002_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_002_Mild_GEM = setParam(COVID19_002_Mild_GEM,'lb',cRxns,cLB);
COVID19_002_Mild_GEM = setParam(COVID19_002_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_002_Mild_GEM = solveLP(COVID19_002_Mild_GEM);
writematrix(sol_COVID19_002_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_002_Mild_GEM.txt')
writecell(COVID19_002_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_002_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_006_Mild.mat')


COVID19_006_Mild_GEM = addReaction(COVID19_006_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_006_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_006_Mild_GEM = simplifyModel(COVID19_006_Mild_GEM);
COVID19_006_Mild_GEM.b = COVID19_006_Mild_GEM.b(:,1);
COVID19_006_Mild_GEM = setParam(COVID19_006_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_006_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_006_Mild_GEM = setParam(COVID19_006_Mild_GEM,'lb',cRxns,cLB);
COVID19_006_Mild_GEM = setParam(COVID19_006_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_006_Mild_GEM = solveLP(COVID19_006_Mild_GEM);
writematrix(sol_COVID19_006_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_006_Mild_GEM.txt')
writecell(COVID19_006_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_006_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_008_Mild.mat')


COVID19_008_Mild_GEM = addReaction(COVID19_008_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_008_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_008_Mild_GEM = simplifyModel(COVID19_008_Mild_GEM);
COVID19_008_Mild_GEM.b = COVID19_008_Mild_GEM.b(:,1);
COVID19_008_Mild_GEM = setParam(COVID19_008_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_008_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_008_Mild_GEM = setParam(COVID19_008_Mild_GEM,'lb',cRxns,cLB);
COVID19_008_Mild_GEM = setParam(COVID19_008_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_008_Mild_GEM = solveLP(COVID19_008_Mild_GEM);
writematrix(sol_COVID19_008_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_008_Mild_GEM.txt')
writecell(COVID19_008_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_008_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_010_Mild.mat')


COVID19_010_Mild_GEM = addReaction(COVID19_010_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_010_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_010_Mild_GEM = simplifyModel(COVID19_010_Mild_GEM);
COVID19_010_Mild_GEM.b = COVID19_010_Mild_GEM.b(:,1);
COVID19_010_Mild_GEM = setParam(COVID19_010_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_010_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_010_Mild_GEM = setParam(COVID19_010_Mild_GEM,'lb',cRxns,cLB);
COVID19_010_Mild_GEM = setParam(COVID19_010_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_010_Mild_GEM = solveLP(COVID19_010_Mild_GEM);
writematrix(sol_COVID19_010_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_010_Mild_GEM.txt')
writecell(COVID19_010_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_010_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_011_Mild.mat')


COVID19_011_Mild_GEM = addReaction(COVID19_011_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_011_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_011_Mild_GEM = simplifyModel(COVID19_011_Mild_GEM);
COVID19_011_Mild_GEM.b = COVID19_011_Mild_GEM.b(:,1);
COVID19_011_Mild_GEM = setParam(COVID19_011_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_011_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_011_Mild_GEM = setParam(COVID19_011_Mild_GEM,'lb',cRxns,cLB);
COVID19_011_Mild_GEM = setParam(COVID19_011_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_011_Mild_GEM = solveLP(COVID19_011_Mild_GEM);
writematrix(sol_COVID19_011_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_011_Mild_GEM.txt')
writecell(COVID19_011_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_011_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_012_Mild.mat')


COVID19_012_Mild_GEM = addReaction(COVID19_012_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_012_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_012_Mild_GEM = simplifyModel(COVID19_012_Mild_GEM);
COVID19_012_Mild_GEM.b = COVID19_012_Mild_GEM.b(:,1);
COVID19_012_Mild_GEM = setParam(COVID19_012_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_012_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_012_Mild_GEM = setParam(COVID19_012_Mild_GEM,'lb',cRxns,cLB);
COVID19_012_Mild_GEM = setParam(COVID19_012_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_012_Mild_GEM = solveLP(COVID19_012_Mild_GEM);
writematrix(sol_COVID19_012_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_012_Mild_GEM.txt')
writecell(COVID19_012_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_012_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_013_Mild.mat')


COVID19_013_Mild_GEM = addReaction(COVID19_013_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_013_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_013_Mild_GEM = simplifyModel(COVID19_013_Mild_GEM);
COVID19_013_Mild_GEM.b = COVID19_013_Mild_GEM.b(:,1);
COVID19_013_Mild_GEM = setParam(COVID19_013_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_013_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_013_Mild_GEM = setParam(COVID19_013_Mild_GEM,'lb',cRxns,cLB);
COVID19_013_Mild_GEM = setParam(COVID19_013_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_013_Mild_GEM = solveLP(COVID19_013_Mild_GEM);
writematrix(sol_COVID19_013_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_013_Mild_GEM.txt')
writecell(COVID19_013_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_013_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_014_Mild.mat')


COVID19_014_Mild_GEM = addReaction(COVID19_014_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_014_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_014_Mild_GEM = simplifyModel(COVID19_014_Mild_GEM);
COVID19_014_Mild_GEM.b = COVID19_014_Mild_GEM.b(:,1);
COVID19_014_Mild_GEM = setParam(COVID19_014_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_014_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_014_Mild_GEM = setParam(COVID19_014_Mild_GEM,'lb',cRxns,cLB);
COVID19_014_Mild_GEM = setParam(COVID19_014_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_014_Mild_GEM = solveLP(COVID19_014_Mild_GEM);
writematrix(sol_COVID19_014_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_014_Mild_GEM.txt')
writecell(COVID19_014_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_014_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_015_Mild.mat')


COVID19_015_Mild_GEM = addReaction(COVID19_015_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_015_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_015_Mild_GEM = simplifyModel(COVID19_015_Mild_GEM);
COVID19_015_Mild_GEM.b = COVID19_015_Mild_GEM.b(:,1);
COVID19_015_Mild_GEM = setParam(COVID19_015_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_015_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_015_Mild_GEM = setParam(COVID19_015_Mild_GEM,'lb',cRxns,cLB);
COVID19_015_Mild_GEM = setParam(COVID19_015_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_015_Mild_GEM = solveLP(COVID19_015_Mild_GEM);
writematrix(sol_COVID19_015_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_015_Mild_GEM.txt')
writecell(COVID19_015_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_015_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_016_Mild.mat')


COVID19_016_Mild_GEM = addReaction(COVID19_016_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_016_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_016_Mild_GEM = simplifyModel(COVID19_016_Mild_GEM);
COVID19_016_Mild_GEM.b = COVID19_016_Mild_GEM.b(:,1);
COVID19_016_Mild_GEM = setParam(COVID19_016_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_016_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_016_Mild_GEM = setParam(COVID19_016_Mild_GEM,'lb',cRxns,cLB);
COVID19_016_Mild_GEM = setParam(COVID19_016_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_016_Mild_GEM = solveLP(COVID19_016_Mild_GEM);
writematrix(sol_COVID19_016_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_016_Mild_GEM.txt')
writecell(COVID19_016_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_016_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_017_Mild.mat')


COVID19_017_Mild_GEM = addReaction(COVID19_017_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_017_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_017_Mild_GEM = simplifyModel(COVID19_017_Mild_GEM);
COVID19_017_Mild_GEM.b = COVID19_017_Mild_GEM.b(:,1);
COVID19_017_Mild_GEM = setParam(COVID19_017_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_017_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_017_Mild_GEM = setParam(COVID19_017_Mild_GEM,'lb',cRxns,cLB);
COVID19_017_Mild_GEM = setParam(COVID19_017_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_017_Mild_GEM = solveLP(COVID19_017_Mild_GEM);
writematrix(sol_COVID19_017_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_017_Mild_GEM.txt')
writecell(COVID19_017_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_017_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_018_Mild.mat')


COVID19_018_Mild_GEM = addReaction(COVID19_018_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_018_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_018_Mild_GEM = simplifyModel(COVID19_018_Mild_GEM);
COVID19_018_Mild_GEM.b = COVID19_018_Mild_GEM.b(:,1);
COVID19_018_Mild_GEM = setParam(COVID19_018_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_018_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_018_Mild_GEM = setParam(COVID19_018_Mild_GEM,'lb',cRxns,cLB);
COVID19_018_Mild_GEM = setParam(COVID19_018_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_018_Mild_GEM = solveLP(COVID19_018_Mild_GEM);
writematrix(sol_COVID19_018_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_018_Mild_GEM.txt')
writecell(COVID19_018_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_018_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_019_Mild.mat')


COVID19_019_Mild_GEM = addReaction(COVID19_019_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_019_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_019_Mild_GEM = simplifyModel(COVID19_019_Mild_GEM);
COVID19_019_Mild_GEM.b = COVID19_019_Mild_GEM.b(:,1);
COVID19_019_Mild_GEM = setParam(COVID19_019_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_019_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_019_Mild_GEM = setParam(COVID19_019_Mild_GEM,'lb',cRxns,cLB);
COVID19_019_Mild_GEM = setParam(COVID19_019_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_019_Mild_GEM = solveLP(COVID19_019_Mild_GEM);
writematrix(sol_COVID19_019_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_019_Mild_GEM.txt')
writecell(COVID19_019_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_019_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_023_Mild.mat')


COVID19_023_Mild_GEM = addReaction(COVID19_023_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_023_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_023_Mild_GEM = simplifyModel(COVID19_023_Mild_GEM);
COVID19_023_Mild_GEM.b = COVID19_023_Mild_GEM.b(:,1);
COVID19_023_Mild_GEM = setParam(COVID19_023_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_023_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_023_Mild_GEM = setParam(COVID19_023_Mild_GEM,'lb',cRxns,cLB);
COVID19_023_Mild_GEM = setParam(COVID19_023_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_023_Mild_GEM = solveLP(COVID19_023_Mild_GEM);
writematrix(sol_COVID19_023_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_023_Mild_GEM.txt')
writecell(COVID19_023_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_023_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_029_Mild.mat')


COVID19_029_Mild_GEM = addReaction(COVID19_029_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_029_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_029_Mild_GEM = simplifyModel(COVID19_029_Mild_GEM);
COVID19_029_Mild_GEM.b = COVID19_029_Mild_GEM.b(:,1);
COVID19_029_Mild_GEM = setParam(COVID19_029_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_029_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_029_Mild_GEM = setParam(COVID19_029_Mild_GEM,'lb',cRxns,cLB);
COVID19_029_Mild_GEM = setParam(COVID19_029_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_029_Mild_GEM = solveLP(COVID19_029_Mild_GEM);
writematrix(sol_COVID19_029_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_029_Mild_GEM.txt')
writecell(COVID19_029_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_029_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_030_Mild.mat')


COVID19_030_Mild_GEM = addReaction(COVID19_030_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_030_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_030_Mild_GEM = simplifyModel(COVID19_030_Mild_GEM);
COVID19_030_Mild_GEM.b = COVID19_030_Mild_GEM.b(:,1);
COVID19_030_Mild_GEM = setParam(COVID19_030_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_030_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_030_Mild_GEM = setParam(COVID19_030_Mild_GEM,'lb',cRxns,cLB);
COVID19_030_Mild_GEM = setParam(COVID19_030_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_030_Mild_GEM = solveLP(COVID19_030_Mild_GEM);
writematrix(sol_COVID19_030_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_030_Mild_GEM.txt')
writecell(COVID19_030_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_030_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_032_Mild.mat')


COVID19_032_Mild_GEM = addReaction(COVID19_032_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_032_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_032_Mild_GEM = simplifyModel(COVID19_032_Mild_GEM);
COVID19_032_Mild_GEM.b = COVID19_032_Mild_GEM.b(:,1);
COVID19_032_Mild_GEM = setParam(COVID19_032_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_032_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_032_Mild_GEM = setParam(COVID19_032_Mild_GEM,'lb',cRxns,cLB);
COVID19_032_Mild_GEM = setParam(COVID19_032_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_032_Mild_GEM = solveLP(COVID19_032_Mild_GEM);
writematrix(sol_COVID19_032_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_032_Mild_GEM.txt')
writecell(COVID19_032_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_032_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_036_Mild.mat')


COVID19_036_Mild_GEM = addReaction(COVID19_036_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_036_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_036_Mild_GEM = simplifyModel(COVID19_036_Mild_GEM);
COVID19_036_Mild_GEM.b = COVID19_036_Mild_GEM.b(:,1);
COVID19_036_Mild_GEM = setParam(COVID19_036_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_036_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_036_Mild_GEM = setParam(COVID19_036_Mild_GEM,'lb',cRxns,cLB);
COVID19_036_Mild_GEM = setParam(COVID19_036_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_036_Mild_GEM = solveLP(COVID19_036_Mild_GEM);
writematrix(sol_COVID19_036_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_036_Mild_GEM.txt')
writecell(COVID19_036_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_036_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_039_Mild.mat')


COVID19_039_Mild_GEM = addReaction(COVID19_039_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_039_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_039_Mild_GEM = simplifyModel(COVID19_039_Mild_GEM);
COVID19_039_Mild_GEM.b = COVID19_039_Mild_GEM.b(:,1);
COVID19_039_Mild_GEM = setParam(COVID19_039_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_039_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_039_Mild_GEM = setParam(COVID19_039_Mild_GEM,'lb',cRxns,cLB);
COVID19_039_Mild_GEM = setParam(COVID19_039_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_039_Mild_GEM = solveLP(COVID19_039_Mild_GEM);
writematrix(sol_COVID19_039_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_039_Mild_GEM.txt')
writecell(COVID19_039_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_039_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_041_Mild.mat')


COVID19_041_Mild_GEM = addReaction(COVID19_041_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_041_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_041_Mild_GEM = simplifyModel(COVID19_041_Mild_GEM);
COVID19_041_Mild_GEM.b = COVID19_041_Mild_GEM.b(:,1);
COVID19_041_Mild_GEM = setParam(COVID19_041_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_041_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_041_Mild_GEM = setParam(COVID19_041_Mild_GEM,'lb',cRxns,cLB);
COVID19_041_Mild_GEM = setParam(COVID19_041_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_041_Mild_GEM = solveLP(COVID19_041_Mild_GEM);
writematrix(sol_COVID19_041_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_041_Mild_GEM.txt')
writecell(COVID19_041_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_041_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_001_Mild.mat')


COVID19_001_Mild_GEM = addReaction(COVID19_001_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_001_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_001_Mild_GEM = simplifyModel(COVID19_001_Mild_GEM);
COVID19_001_Mild_GEM.b = COVID19_001_Mild_GEM.b(:,1);
COVID19_001_Mild_GEM = setParam(COVID19_001_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_001_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_001_Mild_GEM = setParam(COVID19_001_Mild_GEM,'lb',cRxns,cLB);
COVID19_001_Mild_GEM = setParam(COVID19_001_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_001_Mild_GEM = solveLP(COVID19_001_Mild_GEM);
writematrix(sol_COVID19_001_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_001_Mild_GEM.txt')
writecell(COVID19_001_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_001_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_004_Mild.mat')


COVID19_004_Mild_GEM = addReaction(COVID19_004_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_004_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_004_Mild_GEM = simplifyModel(COVID19_004_Mild_GEM);
COVID19_004_Mild_GEM.b = COVID19_004_Mild_GEM.b(:,1);
COVID19_004_Mild_GEM = setParam(COVID19_004_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_004_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_004_Mild_GEM = setParam(COVID19_004_Mild_GEM,'lb',cRxns,cLB);
COVID19_004_Mild_GEM = setParam(COVID19_004_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_004_Mild_GEM = solveLP(COVID19_004_Mild_GEM);
writematrix(sol_COVID19_004_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_004_Mild_GEM.txt')
writecell(COVID19_004_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_004_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_005_Mild.mat')


COVID19_005_Mild_GEM = addReaction(COVID19_005_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_005_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_005_Mild_GEM = simplifyModel(COVID19_005_Mild_GEM);
COVID19_005_Mild_GEM.b = COVID19_005_Mild_GEM.b(:,1);
COVID19_005_Mild_GEM = setParam(COVID19_005_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_005_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_005_Mild_GEM = setParam(COVID19_005_Mild_GEM,'lb',cRxns,cLB);
COVID19_005_Mild_GEM = setParam(COVID19_005_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_005_Mild_GEM = solveLP(COVID19_005_Mild_GEM);
writematrix(sol_COVID19_005_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_005_Mild_GEM.txt')
writecell(COVID19_005_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_005_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_022_Mild.mat')


COVID19_022_Mild_GEM = addReaction(COVID19_022_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_022_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_022_Mild_GEM = simplifyModel(COVID19_022_Mild_GEM);
COVID19_022_Mild_GEM.b = COVID19_022_Mild_GEM.b(:,1);
COVID19_022_Mild_GEM = setParam(COVID19_022_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_022_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_022_Mild_GEM = setParam(COVID19_022_Mild_GEM,'lb',cRxns,cLB);
COVID19_022_Mild_GEM = setParam(COVID19_022_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_022_Mild_GEM = solveLP(COVID19_022_Mild_GEM);
writematrix(sol_COVID19_022_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_022_Mild_GEM.txt')
writecell(COVID19_022_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_022_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_024_Mild.mat')


COVID19_024_Mild_GEM = addReaction(COVID19_024_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_024_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_024_Mild_GEM = simplifyModel(COVID19_024_Mild_GEM);
COVID19_024_Mild_GEM.b = COVID19_024_Mild_GEM.b(:,1);
COVID19_024_Mild_GEM = setParam(COVID19_024_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_024_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_024_Mild_GEM = setParam(COVID19_024_Mild_GEM,'lb',cRxns,cLB);
COVID19_024_Mild_GEM = setParam(COVID19_024_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_024_Mild_GEM = solveLP(COVID19_024_Mild_GEM);
writematrix(sol_COVID19_024_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_024_Mild_GEM.txt')
writecell(COVID19_024_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_024_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_027_Mild.mat')


COVID19_027_Mild_GEM = addReaction(COVID19_027_Mild_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_027_Mild_GEM, [], true, false, false, essentialTasks);


COVID19_027_Mild_GEM = simplifyModel(COVID19_027_Mild_GEM);
COVID19_027_Mild_GEM.b = COVID19_027_Mild_GEM.b(:,1);
COVID19_027_Mild_GEM = setParam(COVID19_027_Mild_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_027_Mild.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_027_Mild_GEM = setParam(COVID19_027_Mild_GEM,'lb',cRxns,cLB);
COVID19_027_Mild_GEM = setParam(COVID19_027_Mild_GEM,'ub',cRxns,cUB);

sol_COVID19_027_Mild_GEM = solveLP(COVID19_027_Mild_GEM);
writematrix(sol_COVID19_027_Mild_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_027_Mild_GEM.txt')
writecell(COVID19_027_Mild_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_027_Mild_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_031_Severe.mat')


COVID19_031_Severe_GEM = addReaction(COVID19_031_Severe_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_031_Severe_GEM, [], true, false, false, essentialTasks);


COVID19_031_Severe_GEM = simplifyModel(COVID19_031_Severe_GEM);
COVID19_031_Severe_GEM.b = COVID19_031_Severe_GEM.b(:,1);
COVID19_031_Severe_GEM = setParam(COVID19_031_Severe_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_031_Severe.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_031_Severe_GEM = setParam(COVID19_031_Severe_GEM,'lb',cRxns,cLB);
COVID19_031_Severe_GEM = setParam(COVID19_031_Severe_GEM,'ub',cRxns,cUB);

sol_COVID19_031_Severe_GEM = solveLP(COVID19_031_Severe_GEM);
writematrix(sol_COVID19_031_Severe_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_031_Severe_GEM.txt')
writecell(COVID19_031_Severe_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_031_Severe_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_009_Severe.mat')


COVID19_009_Severe_GEM = addReaction(COVID19_009_Severe_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_009_Severe_GEM, [], true, false, false, essentialTasks);


COVID19_009_Severe_GEM = simplifyModel(COVID19_009_Severe_GEM);
COVID19_009_Severe_GEM.b = COVID19_009_Severe_GEM.b(:,1);
COVID19_009_Severe_GEM = setParam(COVID19_009_Severe_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_009_Severe.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_009_Severe_GEM = setParam(COVID19_009_Severe_GEM,'lb',cRxns,cLB);
COVID19_009_Severe_GEM = setParam(COVID19_009_Severe_GEM,'ub',cRxns,cUB);

sol_COVID19_009_Severe_GEM = solveLP(COVID19_009_Severe_GEM);
writematrix(sol_COVID19_009_Severe_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_009_Severe_GEM.txt')
writecell(COVID19_009_Severe_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_009_Severe_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_028_Severe.mat')


COVID19_028_Severe_GEM = addReaction(COVID19_028_Severe_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_028_Severe_GEM, [], true, false, false, essentialTasks);


COVID19_028_Severe_GEM = simplifyModel(COVID19_028_Severe_GEM);
COVID19_028_Severe_GEM.b = COVID19_028_Severe_GEM.b(:,1);
COVID19_028_Severe_GEM = setParam(COVID19_028_Severe_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_028_Severe.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_028_Severe_GEM = setParam(COVID19_028_Severe_GEM,'lb',cRxns,cLB);
COVID19_028_Severe_GEM = setParam(COVID19_028_Severe_GEM,'ub',cRxns,cUB);

sol_COVID19_028_Severe_GEM = solveLP(COVID19_028_Severe_GEM);
writematrix(sol_COVID19_028_Severe_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_028_Severe_GEM.txt')
writecell(COVID19_028_Severe_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_028_Severe_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_033_Severe.mat')


COVID19_033_Severe_GEM = addReaction(COVID19_033_Severe_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_033_Severe_GEM, [], true, false, false, essentialTasks);


COVID19_033_Severe_GEM = simplifyModel(COVID19_033_Severe_GEM);
COVID19_033_Severe_GEM.b = COVID19_033_Severe_GEM.b(:,1);
COVID19_033_Severe_GEM = setParam(COVID19_033_Severe_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_033_Severe.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_033_Severe_GEM = setParam(COVID19_033_Severe_GEM,'lb',cRxns,cLB);
COVID19_033_Severe_GEM = setParam(COVID19_033_Severe_GEM,'ub',cRxns,cUB);

sol_COVID19_033_Severe_GEM = solveLP(COVID19_033_Severe_GEM);
writematrix(sol_COVID19_033_Severe_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_033_Severe_GEM.txt')
writecell(COVID19_033_Severe_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_033_Severe_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_034_Severe.mat')


COVID19_034_Severe_GEM = addReaction(COVID19_034_Severe_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_034_Severe_GEM, [], true, false, false, essentialTasks);


COVID19_034_Severe_GEM = simplifyModel(COVID19_034_Severe_GEM);
COVID19_034_Severe_GEM.b = COVID19_034_Severe_GEM.b(:,1);
COVID19_034_Severe_GEM = setParam(COVID19_034_Severe_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_034_Severe.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_034_Severe_GEM = setParam(COVID19_034_Severe_GEM,'lb',cRxns,cLB);
COVID19_034_Severe_GEM = setParam(COVID19_034_Severe_GEM,'ub',cRxns,cUB);

sol_COVID19_034_Severe_GEM = solveLP(COVID19_034_Severe_GEM);
writematrix(sol_COVID19_034_Severe_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_034_Severe_GEM.txt')
writecell(COVID19_034_Severe_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_034_Severe_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_035_Severe.mat')


COVID19_035_Severe_GEM = addReaction(COVID19_035_Severe_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_035_Severe_GEM, [], true, false, false, essentialTasks);


COVID19_035_Severe_GEM = simplifyModel(COVID19_035_Severe_GEM);
COVID19_035_Severe_GEM.b = COVID19_035_Severe_GEM.b(:,1);
COVID19_035_Severe_GEM = setParam(COVID19_035_Severe_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_035_Severe.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_035_Severe_GEM = setParam(COVID19_035_Severe_GEM,'lb',cRxns,cLB);
COVID19_035_Severe_GEM = setParam(COVID19_035_Severe_GEM,'ub',cRxns,cUB);

sol_COVID19_035_Severe_GEM = solveLP(COVID19_035_Severe_GEM);
writematrix(sol_COVID19_035_Severe_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_035_Severe_GEM.txt')
writecell(COVID19_035_Severe_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_035_Severe_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_026_Severe.mat')


COVID19_026_Severe_GEM = addReaction(COVID19_026_Severe_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_026_Severe_GEM, [], true, false, false, essentialTasks);


COVID19_026_Severe_GEM = simplifyModel(COVID19_026_Severe_GEM);
COVID19_026_Severe_GEM.b = COVID19_026_Severe_GEM.b(:,1);
COVID19_026_Severe_GEM = setParam(COVID19_026_Severe_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_026_Severe.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_026_Severe_GEM = setParam(COVID19_026_Severe_GEM,'lb',cRxns,cLB);
COVID19_026_Severe_GEM = setParam(COVID19_026_Severe_GEM,'ub',cRxns,cUB);

sol_COVID19_026_Severe_GEM = solveLP(COVID19_026_Severe_GEM);
writematrix(sol_COVID19_026_Severe_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_026_Severe_GEM.txt')
writecell(COVID19_026_Severe_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_026_Severe_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_037_Severe.mat')


COVID19_037_Severe_GEM = addReaction(COVID19_037_Severe_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_037_Severe_GEM, [], true, false, false, essentialTasks);


COVID19_037_Severe_GEM = simplifyModel(COVID19_037_Severe_GEM);
COVID19_037_Severe_GEM.b = COVID19_037_Severe_GEM.b(:,1);
COVID19_037_Severe_GEM = setParam(COVID19_037_Severe_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_037_Severe.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_037_Severe_GEM = setParam(COVID19_037_Severe_GEM,'lb',cRxns,cLB);
COVID19_037_Severe_GEM = setParam(COVID19_037_Severe_GEM,'ub',cRxns,cUB);

sol_COVID19_037_Severe_GEM = solveLP(COVID19_037_Severe_GEM);
writematrix(sol_COVID19_037_Severe_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_037_Severe_GEM.txt')
writecell(COVID19_037_Severe_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_037_Severe_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_007_Severe.mat')


COVID19_007_Severe_GEM = addReaction(COVID19_007_Severe_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_007_Severe_GEM, [], true, false, false, essentialTasks);


COVID19_007_Severe_GEM = simplifyModel(COVID19_007_Severe_GEM);
COVID19_007_Severe_GEM.b = COVID19_007_Severe_GEM.b(:,1);
COVID19_007_Severe_GEM = setParam(COVID19_007_Severe_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_007_Severe.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_007_Severe_GEM = setParam(COVID19_007_Severe_GEM,'lb',cRxns,cLB);
COVID19_007_Severe_GEM = setParam(COVID19_007_Severe_GEM,'ub',cRxns,cUB);

sol_COVID19_007_Severe_GEM = solveLP(COVID19_007_Severe_GEM);
writematrix(sol_COVID19_007_Severe_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_007_Severe_GEM.txt')
writecell(COVID19_007_Severe_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_007_Severe_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_020_Severe.mat')


COVID19_020_Severe_GEM = addReaction(COVID19_020_Severe_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_020_Severe_GEM, [], true, false, false, essentialTasks);


COVID19_020_Severe_GEM = simplifyModel(COVID19_020_Severe_GEM);
COVID19_020_Severe_GEM.b = COVID19_020_Severe_GEM.b(:,1);
COVID19_020_Severe_GEM = setParam(COVID19_020_Severe_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_020_Severe.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_020_Severe_GEM = setParam(COVID19_020_Severe_GEM,'lb',cRxns,cLB);
COVID19_020_Severe_GEM = setParam(COVID19_020_Severe_GEM,'ub',cRxns,cUB);

sol_COVID19_020_Severe_GEM = solveLP(COVID19_020_Severe_GEM);
writematrix(sol_COVID19_020_Severe_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_020_Severe_GEM.txt')
writecell(COVID19_020_Severe_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_020_Severe_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
load('/home/anoop/Desktop/COVID_Omics/Personal_GEM/COVID19_021_Severe.mat')


COVID19_021_Severe_GEM = addReaction(COVID19_021_Severe_GEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',-1000,'upperBound',0);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(COVID19_021_Severe_GEM, [], true, false, false, essentialTasks);


COVID19_021_Severe_GEM = simplifyModel(COVID19_021_Severe_GEM);
COVID19_021_Severe_GEM.b = COVID19_021_Severe_GEM.b(:,1);
COVID19_021_Severe_GEM = setParam(COVID19_021_Severe_GEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/COVID_Omics/Personal_GEM/Bounds/FC/Converted/COVID19_021_Severe.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

COVID19_021_Severe_GEM = setParam(COVID19_021_Severe_GEM,'lb',cRxns,cLB);
COVID19_021_Severe_GEM = setParam(COVID19_021_Severe_GEM,'ub',cRxns,cUB);

sol_COVID19_021_Severe_GEM = solveLP(COVID19_021_Severe_GEM);
writematrix(sol_COVID19_021_Severe_GEM.x,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Fba_COVID19_021_Severe_GEM.txt')
writecell(COVID19_021_Severe_GEM.rxns,'/home/anoop/Desktop/COVID_Omics/Personal_GEM/Without_VBOF/FBA/Reaction_COVID19_021_Severe_GEM.txt')

clear
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%