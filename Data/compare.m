%% demo_mixtures
clear;

%% Initialize variables.
filename = './Data/mix_data_50.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
Al = dataArray{:, 1};
As = dataArray{:, 2};
Ba = dataArray{:, 3};
bc = dataArray{:, 4};
Br = dataArray{:, 5};
Ca = dataArray{:, 6};
Cl = dataArray{:, 7};
Cr = dataArray{:, 8};
Cu = dataArray{:, 9};
Fe = dataArray{:, 10};
K = dataArray{:, 11};
Mn = dataArray{:, 12};
Ni = dataArray{:, 13};
Pb = dataArray{:, 14};
S = dataArray{:, 15};
Se = dataArray{:, 16};
Si = dataArray{:, 17};
Ti = dataArray{:, 18};
V = dataArray{:, 19};
Zn = dataArray{:, 20};

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%% BEGIN PCP
%X = [pm25 pm1 Al As Ba bc Br Ca Cl Cr Cu Fe K Mn Ni Pb S Se Si Ti V Zn];
X = [Al As Ba bc Br Ca Cl Cr Cu Fe K  Mn  Ni  Pb  S  Se  Si Ti  V Zn];
n = [ 1  2  3  4  5  6  7  8  9 10 11 12  13  14  15 16  17 18 19 20];
% [] makes a vector, take columns from data set, put into matrix as columns

numMissingPerRow = sum( isnan(X), 2 ); 
%get rid of rows with NANs
goodRows = find( numMissingPerRow == 0 ); 
% good rows without missing data

X = X(goodRows,:); 
%semicolon means it doesnt output the results

[m,n] = size(X);

lambda = 1/sqrt(m); 
%weight parameter for pcp, 1 / sqrt(number of rows)

%% Vector LOD
% delta = vector of 50% quantiles for 20 variables
delta50 = [0.04195, 3e-04, 0.0066, 0.5733649, 2e-04, 0.02645, 0.0024, ...
    4e-04, 0.003, 0.0537, 0.0326, 0, 0.0015, 0.005, 0.72385, 0, 0.0574, ...
    0.0028, 0.0019, 0.0081];

[Lv, Sv, lossv] = pcp_lod(X, lambda, 10, delta50); 

save('./Data/lowrank_50v.mat', 'Lv')

save('./Data/sparse_50v.mat', 'Sv')

%% Scalar LOD
[Ls,Ss, losss] = pcp_lod(X, lambda, 10, 0.01); 

save('./Data/lowrank_50s.mat', 'Ls')
save('./Data/sparse_50s.mat', 'Ss')

%% Matrix LOD
% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 20);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Al", "As", "Ba", "bc", "Br", "Ca", "Cl", "Cr", "Cu", "Fe", "K", "Mn", "Ni", "Pb", "S", "Se", "Si", "Ti", "V", "Zn"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
mixdatalod50 = readtable("./Data/matrix50.csv", opts);

%% Convert to output type
matrix50 = table2array(mixdatalod50);

%% Clear temporary variables
clear opts

[Lm,Sm, lossm] = pcp_lod(X, lambda, 10, matrix50); 

save('./Data/lowrank_50m.mat', 'Lm')
save('./Data/sparse_50m.mat', 'Sm')
