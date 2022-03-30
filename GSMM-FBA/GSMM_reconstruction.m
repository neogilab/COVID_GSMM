addpath(genpath('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/Test/Human-GEM/'))
addpath(genpath('/home/anoop/Downloads/RAVEN-2.4.0'))

initCobraToolbox(false) 
changeCobraSolver ('gurobi', 'all');

COVID19_7_HC_gex= readtable('TPM_SampleID_GeneName.txt');
COVID19_7_HC.tissues={'pbmc'};
COVID19_7_HC.genes=COVID19_7_HC_gex.ID;
COVID19_7_HC.levels=table2array(COVID19_7_HC_gex(:,'COVID19_7_HC'));
COVID19_7_HC.threshold=1;


Params.timeLimit = 4500.0;

load('Human-GEM.mat');
ihuman = addBoundaryMets(ihuman);
essentialTasks = parseTaskList('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/Test/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.xlsx');
checkTasks(ihuman, [], true, false, false, essentialTasks);

refModel = ihuman;  % the reference model from which the GEM will be extracted
tissue = 'pbmc';  % must match the tissue name in data_struct.tissues
celltype = [];  % used if tissues are subdivided into cell type, which is not the case here
hpaData = [];  % data structure containing protein abundance information (not used here)
arrayData = COVID19_7_HC;  % data structure with gene (RNA) abundance information
metabolomicsData = [];  % list of metabolite names if metabolomic data is available
removeGenes = true;  % (default) remove lowly/non-expressed genes from the extracted GEM
taskFile = [];  % we already loaded the task file, so this input is not required
useScoresForTasks = true;  % (default) use expression data to decide which reactions to keep
printReport = true;  % (default) print status/completion report to screen
taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
params =Params;  % additional optimization parameters for the INIT algorithm
paramsFT = [];
COVID19_7_HC_GEM = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);

COVID19_7_HC_GEM.id='pbmc_COVID19_7_HC'
checkTasks(COVID19_7_HC_GEM, [], true, false, false, essentialTasks)
save('COVID19_7_HC.mat','COVID19_7_HC_GEM')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COVID19_8_HC_gex= readtable('TPM_SampleID_GeneName.txt');
COVID19_8_HC.tissues={'pbmc'};
COVID19_8_HC.genes=COVID19_8_HC_gex.ID;
COVID19_8_HC.levels=table2array(COVID19_8_HC_gex(:,'COVID19_8_HC'));
COVID19_8_HC.threshold=1;

Params.timeLimit = 4500.0;

load('Human-GEM.mat');
ihuman = addBoundaryMets(ihuman);
essentialTasks = parseTaskList('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/Test/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.xlsx');
checkTasks(ihuman, [], true, false, false, essentialTasks);

refModel = ihuman;  % the reference model from which the GEM will be extracted
tissue = 'pbmc';  % must match the tissue name in data_struct.tissues
celltype = [];  % used if tissues are subdivided into cell type, which is not the case here
hpaData = [];  % data structure containing protein abundance information (not used here)
arrayData = COVID19_8_HC;  % data structure with gene (RNA) abundance information
metabolomicsData = [];  % list of metabolite names if metabolomic data is available
removeGenes = true;  % (default) remove lowly/non-expressed genes from the extracted GEM
taskFile = [];  % we already loaded the task file, so this input is not required
useScoresForTasks = true;  % (default) use expression data to decide which reactions to keep
printReport = true;  % (default) print status/completion report to screen
taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
params =Params;  % additional optimization parameters for the INIT algorithm
paramsFT = [];
COVID19_8_HC_GEM = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);

COVID19_8_HC_GEM.id='pbmc_COVID19_8_HC'
checkTasks(COVID19_8_HC_GEM, [], true, false, false, essentialTasks)
save('COVID19_8_HC.mat','COVID19_8_HC_GEM')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COVID19_9_HC_gex= readtable('TPM_SampleID_GeneName.txt');
COVID19_9_HC.tissues={'pbmc'};
COVID19_9_HC.genes=COVID19_9_HC_gex.ID;
COVID19_9_HC.levels=table2array(COVID19_9_HC_gex(:,'COVID19_9_HC'));
COVID19_9_HC.threshold=1;

Params.timeLimit = 4500.0;

load('Human-GEM.mat');
ihuman = addBoundaryMets(ihuman);
essentialTasks = parseTaskList('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/Test/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.xlsx');
checkTasks(ihuman, [], true, false, false, essentialTasks);

refModel = ihuman;  % the reference model from which the GEM will be extracted
tissue = 'pbmc';  % must match the tissue name in data_struct.tissues
celltype = [];  % used if tissues are subdivided into cell type, which is not the case here
hpaData = [];  % data structure containing protein abundance information (not used here)
arrayData = COVID19_9_HC;  % data structure with gene (RNA) abundance information
metabolomicsData = [];  % list of metabolite names if metabolomic data is available
removeGenes = true;  % (default) remove lowly/non-expressed genes from the extracted GEM
taskFile = [];  % we already loaded the task file, so this input is not required
useScoresForTasks = true;  % (default) use expression data to decide which reactions to keep
printReport = true;  % (default) print status/completion report to screen
taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
params =Params;  % additional optimization parameters for the INIT algorithm
paramsFT = [];
COVID19_9_HC_GEM = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);

COVID19_9_HC_GEM.id='pbmc_COVID19_9_HC'
checkTasks(COVID19_9_HC_GEM, [], true, false, false, essentialTasks)
save('COVID19_9_HC.mat','COVID19_9_HC_GEM')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COVID19_11_HC_gex= readtable('TPM_SampleID_GeneName.txt');
COVID19_11_HC.tissues={'pbmc'};
COVID19_11_HC.genes=COVID19_11_HC_gex.ID;
COVID19_11_HC.levels=table2array(COVID19_11_HC_gex(:,'COVID19_11_HC'));
COVID19_11_HC.threshold=1;

Params.timeLimit = 4500.0;

load('Human-GEM.mat');
ihuman = addBoundaryMets(ihuman);
essentialTasks = parseTaskList('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/Test/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.xlsx');
checkTasks(ihuman, [], true, false, false, essentialTasks);

refModel = ihuman;  % the reference model from which the GEM will be extracted
tissue = 'pbmc';  % must match the tissue name in data_struct.tissues
celltype = [];  % used if tissues are subdivided into cell type, which is not the case here
hpaData = [];  % data structure containing protein abundance information (not used here)
arrayData = COVID19_11_HC;  % data structure with gene (RNA) abundance information
metabolomicsData = [];  % list of metabolite names if metabolomic data is available
removeGenes = true;  % (default) remove lowly/non-expressed genes from the extracted GEM
taskFile = [];  % we already loaded the task file, so this input is not required
useScoresForTasks = true;  % (default) use expression data to decide which reactions to keep
printReport = true;  % (default) print status/completion report to screen
taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
params =Params;  % additional optimization parameters for the INIT algorithm
paramsFT = [];
COVID19_11_HC_GEM = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);

COVID19_11_HC_GEM.id='pbmc_COVID19_11_HC'
checkTasks(COVID19_11_HC_GEM, [], true, false, false, essentialTasks)
save('COVID19_11_HC.mat','COVID19_11_HC_GEM')

clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COVID19_16_HC_gex= readtable('TPM_SampleID_GeneName.txt');
COVID19_16_HC.tissues={'pbmc'};
COVID19_16_HC.genes=COVID19_16_HC_gex.ID;
COVID19_16_HC.levels=table2array(COVID19_16_HC_gex(:,'COVID19_16_HC'));
COVID19_16_HC.threshold=1;

Params.timeLimit = 4500.0;

load('Human-GEM.mat');
ihuman = addBoundaryMets(ihuman);
essentialTasks = parseTaskList('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/Test/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.xlsx');
checkTasks(ihuman, [], true, false, false, essentialTasks);

refModel = ihuman;  % the reference model from which the GEM will be extracted
tissue = 'pbmc';  % must match the tissue name in data_struct.tissues
celltype = [];  % used if tissues are subdivided into cell type, which is not the case here
hpaData = [];  % data structure containing protein abundance information (not used here)
arrayData = COVID19_16_HC;  % data structure with gene (RNA) abundance information
metabolomicsData = [];  % list of metabolite names if metabolomic data is available
removeGenes = true;  % (default) remove lowly/non-expressed genes from the extracted GEM
taskFile = [];  % we already loaded the task file, so this input is not required
useScoresForTasks = true;  % (default) use expression data to decide which reactions to keep
printReport = true;  % (default) print status/completion report to screen
taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
params =Params;  % additional optimization parameters for the INIT algorithm
paramsFT = [];
COVID19_16_HC_GEM = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);

COVID19_16_HC_GEM.id='pbmc_COVID19_16_HC'
checkTasks(COVID19_16_HC_GEM, [], true, false, false, essentialTasks)
save('COVID19_16_HC.mat','COVID19_16_HC_GEM')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


COVID19_17_HC_gex= readtable('TPM_SampleID_GeneName.txt');
COVID19_17_HC.tissues={'pbmc'};
COVID19_17_HC.genes=COVID19_17_HC_gex.ID;
COVID19_17_HC.levels=table2array(COVID19_17_HC_gex(:,'COVID19_17_HC'));
COVID19_17_HC.threshold=1;

Params.timeLimit = 4500.0;

load('Human-GEM.mat');
ihuman = addBoundaryMets(ihuman);
essentialTasks = parseTaskList('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/Test/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.xlsx');
checkTasks(ihuman, [], true, false, false, essentialTasks);

refModel = ihuman;  % the reference model from which the GEM will be extracted
tissue = 'pbmc';  % must match the tissue name in data_struct.tissues
celltype = [];  % used if tissues are subdivided into cell type, which is not the case here
hpaData = [];  % data structure containing protein abundance information (not used here)
arrayData = COVID19_17_HC;  % data structure with gene (RNA) abundance information
metabolomicsData = [];  % list of metabolite names if metabolomic data is available
removeGenes = true;  % (default) remove lowly/non-expressed genes from the extracted GEM
taskFile = [];  % we already loaded the task file, so this input is not required
useScoresForTasks = true;  % (default) use expression data to decide which reactions to keep
printReport = true;  % (default) print status/completion report to screen
taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
params =Params;  % additional optimization parameters for the INIT algorithm
paramsFT = [];
COVID19_17_HC_GEM = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);

COVID19_17_HC_GEM.id='pbmc_COVID19_17_HC'
checkTasks(COVID19_17_HC_GEM, [], true, false, false, essentialTasks)
save('COVID19_17_HC.mat','COVID19_17_HC_GEM')

clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COVID19_19_HC_gex= readtable('TPM_SampleID_GeneName.txt');
COVID19_19_HC.tissues={'pbmc'};
COVID19_19_HC.genes=COVID19_19_HC_gex.ID;
COVID19_19_HC.levels=table2array(COVID19_19_HC_gex(:,'COVID19_19_HC'));
COVID19_19_HC.threshold=1;

Params.timeLimit = 4500.0;

load('Human-GEM.mat');
ihuman = addBoundaryMets(ihuman);
essentialTasks = parseTaskList('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/Test/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.xlsx');
checkTasks(ihuman, [], true, false, false, essentialTasks);

refModel = ihuman;  % the reference model from which the GEM will be extracted
tissue = 'pbmc';  % must match the tissue name in data_struct.tissues
celltype = [];  % used if tissues are subdivided into cell type, which is not the case here
hpaData = [];  % data structure containing protein abundance information (not used here)
arrayData = COVID19_19_HC;  % data structure with gene (RNA) abundance information
metabolomicsData = [];  % list of metabolite names if metabolomic data is available
removeGenes = true;  % (default) remove lowly/non-expressed genes from the extracted GEM
taskFile = [];  % we already loaded the task file, so this input is not required
useScoresForTasks = true;  % (default) use expression data to decide which reactions to keep
printReport = true;  % (default) print status/completion report to screen
taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
params =Params;  % additional optimization parameters for the INIT algorithm
paramsFT = [];
COVID19_19_HC_GEM = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);

COVID19_19_HC_GEM.id='pbmc_COVID19_19_HC'
checkTasks(COVID19_19_HC_GEM, [], true, false, false, essentialTasks)
save('COVID19_19_HC.mat','COVID19_19_HC_GEM')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COVID19_20_HC_gex= readtable('TPM_SampleID_GeneName.txt');
COVID19_20_HC.tissues={'pbmc'};
COVID19_20_HC.genes=COVID19_20_HC_gex.ID;
COVID19_20_HC.levels=table2array(COVID19_20_HC_gex(:,'COVID19_20_HC'));
COVID19_20_HC.threshold=1;

Params.timeLimit = 4500.0;

load('Human-GEM.mat');
ihuman = addBoundaryMets(ihuman);
essentialTasks = parseTaskList('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/Test/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.xlsx');
checkTasks(ihuman, [], true, false, false, essentialTasks);

refModel = ihuman;  % the reference model from which the GEM will be extracted
tissue = 'pbmc';  % must match the tissue name in data_struct.tissues
celltype = [];  % used if tissues are subdivided into cell type, which is not the case here
hpaData = [];  % data structure containing protein abundance information (not used here)
arrayData = COVID19_20_HC;  % data structure with gene (RNA) abundance information
metabolomicsData = [];  % list of metabolite names if metabolomic data is available
removeGenes = true;  % (default) remove lowly/non-expressed genes from the extracted GEM
taskFile = [];  % we already loaded the task file, so this input is not required
useScoresForTasks = true;  % (default) use expression data to decide which reactions to keep
printReport = true;  % (default) print status/completion report to screen
taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
params =Params;  % additional optimization parameters for the INIT algorithm
paramsFT = [];
COVID19_20_HC_GEM = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);

COVID19_20_HC_GEM.id='pbmc_COVID19_20_HC'
checkTasks(COVID19_20_HC_GEM, [], true, false, false, essentialTasks)
save('COVID19_20_HC.mat','COVID19_20_HC_GEM')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COVID19_21_HC_gex= readtable('TPM_SampleID_GeneName.txt');
COVID19_21_HC.tissues={'pbmc'};
COVID19_21_HC.genes=COVID19_21_HC_gex.ID;
COVID19_21_HC.levels=table2array(COVID19_21_HC_gex(:,'COVID19_21_HC'));
COVID19_21_HC.threshold=1;

Params.timeLimit = 4500.0;

load('Human-GEM.mat');
ihuman = addBoundaryMets(ihuman);
essentialTasks = parseTaskList('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/Test/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.xlsx');
checkTasks(ihuman, [], true, false, false, essentialTasks);

refModel = ihuman;  % the reference model from which the GEM will be extracted
tissue = 'pbmc';  % must match the tissue name in data_struct.tissues
celltype = [];  % used if tissues are subdivided into cell type, which is not the case here
hpaData = [];  % data structure containing protein abundance information (not used here)
arrayData = COVID19_21_HC;  % data structure with gene (RNA) abundance information
metabolomicsData = [];  % list of metabolite names if metabolomic data is available
removeGenes = true;  % (default) remove lowly/non-expressed genes from the extracted GEM
taskFile = [];  % we already loaded the task file, so this input is not required
useScoresForTasks = true;  % (default) use expression data to decide which reactions to keep
printReport = true;  % (default) print status/completion report to screen
taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
params =Params;  % additional optimization parameters for the INIT algorithm
paramsFT = [];
COVID19_21_HC_GEM = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);

COVID19_21_HC_GEM.id='pbmc_COVID19_21_HC'
checkTasks(COVID19_21_HC_GEM, [], true, false, false, essentialTasks)
save('COVID19_21_HC.mat','COVID19_21_HC_GEM')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COVID19_22_HC_gex= readtable('TPM_SampleID_GeneName.txt');
COVID19_22_HC.tissues={'pbmc'};
COVID19_22_HC.genes=COVID19_22_HC_gex.ID;
COVID19_22_HC.levels=table2array(COVID19_22_HC_gex(:,'COVID19_22_HC'));
COVID19_22_HC.threshold=1;

Params.timeLimit = 4500.0;

load('Human-GEM.mat');
ihuman = addBoundaryMets(ihuman);
essentialTasks = parseTaskList('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/Test/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.xlsx');
checkTasks(ihuman, [], true, false, false, essentialTasks);

refModel = ihuman;  % the reference model from which the GEM will be extracted
tissue = 'pbmc';  % must match the tissue name in data_struct.tissues
celltype = [];  % used if tissues are subdivided into cell type, which is not the case here
hpaData = [];  % data structure containing protein abundance information (not used here)
arrayData = COVID19_22_HC;  % data structure with gene (RNA) abundance information
metabolomicsData = [];  % list of metabolite names if metabolomic data is available
removeGenes = true;  % (default) remove lowly/non-expressed genes from the extracted GEM
taskFile = [];  % we already loaded the task file, so this input is not required
useScoresForTasks = true;  % (default) use expression data to decide which reactions to keep
printReport = true;  % (default) print status/completion report to screen
taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
params =Params;  % additional optimization parameters for the INIT algorithm
paramsFT = [];
COVID19_22_HC_GEM = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);

COVID19_22_HC_GEM.id='pbmc_COVID19_22_HC'
checkTasks(COVID19_22_HC_GEM, [], true, false, false, essentialTasks)
save('COVID19_22_HC.mat','COVID19_22_HC_GEM')

%%%%%%%%%%%%%%%%%%%%%%%
COVID19_25_HC_gex= readtable('TPM_SampleID_GeneName.txt');
COVID19_25_HC.tissues={'pbmc'};
COVID19_25_HC.genes=COVID19_25_HC_gex.ID;
COVID19_25_HC.levels=table2array(COVID19_25_HC_gex(:,'COVID19_25_HC'));
COVID19_25_HC.threshold=1;

Params.timeLimit = 4500.0;

load('Human-GEM.mat');
ihuman = addBoundaryMets(ihuman);
essentialTasks = parseTaskList('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/Test/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.xlsx');
checkTasks(ihuman, [], true, false, false, essentialTasks);

refModel = ihuman;  % the reference model from which the GEM will be extracted
tissue = 'pbmc';  % must match the tissue name in data_struct.tissues
celltype = [];  % used if tissues are subdivided into cell type, which is not the case here
hpaData = [];  % data structure containing protein abundance information (not used here)
arrayData = COVID19_25_HC;  % data structure with gene (RNA) abundance information
metabolomicsData = [];  % list of metabolite names if metabolomic data is available
removeGenes = true;  % (default) remove lowly/non-expressed genes from the extracted GEM
taskFile = [];  % we already loaded the task file, so this input is not required
useScoresForTasks = true;  % (default) use expression data to decide which reactions to keep
printReport = true;  % (default) print status/completion report to screen
taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
params =Params;  % additional optimization parameters for the INIT algorithm
paramsFT = [];
COVID19_25_HC_GEM = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);

COVID19_25_HC_GEM.id='pbmc_COVID19_25_HC'
checkTasks(COVID19_25_HC_GEM, [], true, false, false, essentialTasks)
save('COVID19_25_HC.mat','COVID19_25_HC_GEM')

clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COVID19_27_HC_gex= readtable('TPM_SampleID_GeneName.txt');
COVID19_27_HC.tissues={'pbmc'};
COVID19_27_HC.genes=COVID19_27_HC_gex.ID;
COVID19_27_HC.levels=table2array(COVID19_27_HC_gex(:,'COVID19_27_HC'));
COVID19_27_HC.threshold=1;

Params.timeLimit = 4500.0;

load('Human-GEM.mat');
ihuman = addBoundaryMets(ihuman);
essentialTasks = parseTaskList('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/Test/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.xlsx');
checkTasks(ihuman, [], true, false, false, essentialTasks);

refModel = ihuman;  % the reference model from which the GEM will be extracted
tissue = 'pbmc';  % must match the tissue name in data_struct.tissues
celltype = [];  % used if tissues are subdivided into cell type, which is not the case here
hpaData = [];  % data structure containing protein abundance information (not used here)
arrayData = COVID19_27_HC;  % data structure with gene (RNA) abundance information
metabolomicsData = [];  % list of metabolite names if metabolomic data is available
removeGenes = true;  % (default) remove lowly/non-expressed genes from the extracted GEM
taskFile = [];  % we already loaded the task file, so this input is not required
useScoresForTasks = true;  % (default) use expression data to decide which reactions to keep
printReport = true;  % (default) print status/completion report to screen
taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
params =Params;  % additional optimization parameters for the INIT algorithm
paramsFT = [];
COVID19_27_HC_GEM = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);

COVID19_27_HC_GEM.id='pbmc_COVID19_27_HC'
checkTasks(COVID19_27_HC_GEM, [], true, false, false, essentialTasks)
save('COVID19_27_HC.mat','COVID19_27_HC_GEM')

%%%%%%%%%%%%%%%%%%%%%

COVID19_28_HC_gex= readtable('TPM_SampleID_GeneName.txt');
COVID19_28_HC.tissues={'pbmc'};
COVID19_28_HC.genes=COVID19_28_HC_gex.ID;
COVID19_28_HC.levels=table2array(COVID19_28_HC_gex(:,'COVID19_28_HC'));
COVID19_28_HC.threshold=1;

Params.timeLimit = 4500.0;

load('Human-GEM.mat');
ihuman = addBoundaryMets(ihuman);
essentialTasks = parseTaskList('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/Test/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.xlsx');
checkTasks(ihuman, [], true, false, false, essentialTasks);

refModel = ihuman;  % the reference model from which the GEM will be extracted
tissue = 'pbmc';  % must match the tissue name in data_struct.tissues
celltype = [];  % used if tissues are subdivided into cell type, which is not the case here
hpaData = [];  % data structure containing protein abundance information (not used here)
arrayData = COVID19_28_HC;  % data structure with gene (RNA) abundance information
metabolomicsData = [];  % list of metabolite names if metabolomic data is available
removeGenes = true;  % (default) remove lowly/non-expressed genes from the extracted GEM
taskFile = [];  % we already loaded the task file, so this input is not required
useScoresForTasks = true;  % (default) use expression data to decide which reactions to keep
printReport = true;  % (default) print status/completion report to screen
taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
params =Params;  % additional optimization parameters for the INIT algorithm
paramsFT = [];
COVID19_28_HC_GEM = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);

COVID19_28_HC_GEM.id='pbmc_COVID19_28_HC'
checkTasks(COVID19_28_HC_GEM, [], true, false, false, essentialTasks)
save('COVID19_28_HC.mat','COVID19_28_HC_GEM')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COVID19_29_HC_gex= readtable('TPM_SampleID_GeneName.txt');
COVID19_29_HC.tissues={'pbmc'};
COVID19_29_HC.genes=COVID19_29_HC_gex.ID;
COVID19_29_HC.levels=table2array(COVID19_29_HC_gex(:,'COVID19_29_HC'));
COVID19_29_HC.threshold=1;

Params.timeLimit = 4500.0;

load('Human-GEM.mat');
ihuman = addBoundaryMets(ihuman);
essentialTasks = parseTaskList('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/Test/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.xlsx');
checkTasks(ihuman, [], true, false, false, essentialTasks);

refModel = ihuman;  % the reference model from which the GEM will be extracted
tissue = 'pbmc';  % must match the tissue name in data_struct.tissues
celltype = [];  % used if tissues are subdivided into cell type, which is not the case here
hpaData = [];  % data structure containing protein abundance information (not used here)
arrayData = COVID19_29_HC;  % data structure with gene (RNA) abundance information
metabolomicsData = [];  % list of metabolite names if metabolomic data is available
removeGenes = true;  % (default) remove lowly/non-expressed genes from the extracted GEM
taskFile = [];  % we already loaded the task file, so this input is not required
useScoresForTasks = true;  % (default) use expression data to decide which reactions to keep
printReport = true;  % (default) print status/completion report to screen
taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
params =Params;  % additional optimization parameters for the INIT algorithm
paramsFT = [];
COVID19_29_HC_GEM = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);

COVID19_29_HC_GEM.id='pbmc_COVID19_29_HC'
checkTasks(COVID19_29_HC_GEM, [], true, false, false, essentialTasks)
save('COVID19_29_HC.mat','COVID19_29_HC_GEM')

%%%%%%%%%%%%%%%%%%%%%%%

COVID19_24_HC_gex= readtable('TPM_SampleID_GeneName.txt');
COVID19_24_HC.tissues={'pbmc'};
COVID19_24_HC.genes=COVID19_24_HC_gex.ID;
COVID19_24_HC.levels=table2array(COVID19_24_HC_gex(:,'COVID19_24_HC'));
COVID19_24_HC.threshold=1;

Params.timeLimit = 4500.0;

load('Human-GEM.mat');
ihuman = addBoundaryMets(ihuman);
essentialTasks = parseTaskList('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/Test/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.xlsx');
checkTasks(ihuman, [], true, false, false, essentialTasks);

refModel = ihuman;  % the reference model from which the GEM will be extracted
tissue = 'pbmc';  % must match the tissue name in data_struct.tissues
celltype = [];  % used if tissues are subdivided into cell type, which is not the case here
hpaData = [];  % data structure containing protein abundance information (not used here)
arrayData = COVID19_24_HC;  % data structure with gene (RNA) abundance information
metabolomicsData = [];  % list of metabolite names if metabolomic data is available
removeGenes = true;  % (default) remove lowly/non-expressed genes from the extracted GEM
taskFile = [];  % we already loaded the task file, so this input is not required
useScoresForTasks = true;  % (default) use expression data to decide which reactions to keep
printReport = true;  % (default) print status/completion report to screen
taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
params =Params;  % additional optimization parameters for the INIT algorithm
paramsFT = [];
COVID19_24_HC_GEM = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);

COVID19_24_HC_GEM.id='pbmc_COVID19_24_HC'
checkTasks(COVID19_24_HC_GEM, [], true, false, false, essentialTasks)
save('COVID19_24_HC.mat','COVID19_24_HC_GEM')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COVID19_31_HC_gex= readtable('TPM_SampleID_GeneName.txt');
COVID19_31_HC.tissues={'pbmc'};
COVID19_31_HC.genes=COVID19_31_HC_gex.ID;
COVID19_31_HC.levels=table2array(COVID19_31_HC_gex(:,'COVID19_31_HC'));
COVID19_31_HC.threshold=1;

Params.timeLimit = 4500.0;

load('Human-GEM.mat');
ihuman = addBoundaryMets(ihuman);
essentialTasks = parseTaskList('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/Test/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.xlsx');
checkTasks(ihuman, [], true, false, false, essentialTasks);

refModel = ihuman;  % the reference model from which the GEM will be extracted
tissue = 'pbmc';  % must match the tissue name in data_struct.tissues
celltype = [];  % used if tissues are subdivided into cell type, which is not the case here
hpaData = [];  % data structure containing protein abundance information (not used here)
arrayData = COVID19_31_HC;  % data structure with gene (RNA) abundance information
metabolomicsData = [];  % list of metabolite names if metabolomic data is available
removeGenes = true;  % (default) remove lowly/non-expressed genes from the extracted GEM
taskFile = [];  % we already loaded the task file, so this input is not required
useScoresForTasks = true;  % (default) use expression data to decide which reactions to keep
printReport = true;  % (default) print status/completion report to screen
taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
params =Params;  % additional optimization parameters for the INIT algorithm
paramsFT = [];
COVID19_31_HC_GEM = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);

COVID19_31_HC_GEM.id='pbmc_CCOVID19_31_HC'
checkTasks(COVID19_31_HC_GEM, [], true, false, false, essentialTasks)
save('COVID19_31_HC.mat','COVID19_31_HC_GEM')

clear