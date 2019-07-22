%Code for Bonferroni Corrected Pearson's Correlation Coefficient network
analysis
%This section clears all existing entities in the workspace, command window
lines and open menus and windows to prevent potential overlapping of data.
clc
close all
clear all
%This section loads in all the data files: DatasetName =
importdata('filename.csv')
%This section generates Pearson's Correlation Coefficient symmetry
correlation matrices with corresponding p-values.
tic
MasterData = importdata ('22_MasterData.csv');
[R1,P1]= corrcoef(MasterData,'Rows','pairwise'); %R1 is the correlation
matrix for MasterData set, P1 is the probability matrix for MasterData set.
Surv6M = importdata ('22_6MData.csv');
[R2,P2]= corrcoef(Surv6M,'Rows','pairwise'); %R2 is the correlation matrix
for Surv6M set, P2 is the probability matrix for Surv6M.
NonSurv6M = importdata ('22_6MNData.csv');
[R3,P3]= corrcoef(NonSurv6M,'Rows','pairwise'); %R3 is the correlation
matrix for NonSurv6M set, P3 is the probability matrix for NonSurv6M set.
Surv12M = importdata ('22_12MData.csv');
[R4,P4]= corrcoef(Surv12M,'Rows','pairwise'); %R4 is the correlation matrix
for Surv12M set, P4 is the probability matrix for Surv12M set.
NonSurv12M = importdata ('22_12MNData.csv');
[R5,P5]= corrcoef(NonSurv12M,'Rows','pairwise'); %R5 is the correlation
matrix for NonSurv12M set, P5 is the probability matrix for NonSurv12M set.
Surv18M = importdata ('22_18MData.csv');
[R6,P6]= corrcoef(Surv18M,'Rows','pairwise'); %R6 is the correlation matrix
for Surv18M set, P6 is the probability matrix for Surv18M set.
NonSurv18M = importdata ('22_18MNData.csv');
[R7,P7]= corrcoef(NonSurv18M,'Rows','pairwise'); %R7 is the correlation
matrix for NonSurv18M set, P7 is the probability matrix for NonSurv18M set.
toc
%This section filters and removes non-significant correlations from the
symmetry matrix and generates the network graph of the survivor and nonsurvivor group for a specific time period. The code shown is generating the
network graphs for the time period of 6 months.
K=231 %This is the Bonferroni correction K value which is applied to the pvalue during filtering.
P2B=NaN(22); %Generates new matrix of 22x22 NaN values P2B is a matrix of
31x31 of NaN values to prevent overriding of the original matrix .
for j=1:22
 for i=1:22
 if P2(i,j)<=0.05/K; %Condition: if value = x, then replace with y,
else replace with z in NaN matrix.
 P2B(i,j)=R2(i,j);
54
 else P2B(i,j)=0;
 end
 end
end
P2Bgraph=graph(P2B)
P2BAbsgraph=abs(P2B); %All edges in the graph should have a positive
weight.
P2BRgraph=graph(P2BAbsgraph) %P2BRgraph is the absolute value transformed
graph of filtered symmetry matrix P2B.
Weight=P2BAbsgraph;
LWidths=5*P2BRgraph.Edges.Weight/max(P2BRgraph.Edges.Weight);
node_names = {'CRP','IL-6','TNFa','Pugh','MELD','HR','SDNN','SD2','Glucose','Urea','Creatinine','Na','Ammo
nia','Tot Prot','Albumin','Tot
Bili','PT_pc','INR','Hb','Plt','ALT','Indolo'};
P2BRgraph = graph(P2BAbsgraph,node_names);
subplot (1,2,1)
plot(P2BRgraph,'LineWidth',LWidths,'Layout','force','UseGravity',true,'Node
Color','green','EdgeColor','black','MarkerSize',4)
%"'EdgeLabel',P2BRgraph.Edges.Weight," - option to label edges with weight
in the plot function code.
title ('Survivor 6M')
P2BRgraph.Nodes;
P2BR_degree = centrality(P2BRgraph,'Degree'); %Degree is the number of
edges connecting to each node.
P2BRgraph.Nodes.degree = P2BR_degree;
P2BR_betweenness = centrality(P2BRgraph,'Betweenness'); %Betweenness
measures how often each graph node appears on a shortest path between two
nodes in the graph.
P2BRgraph.Nodes.betweenness = P2BR_betweenness;
P2BR_closeness = centrality(P2BRgraph,'closeness'); %Closeness is the
inverse sum of the distance from a node to all other nodes in the graph
P2BRgraph.Nodes.closeness = P2BR_closeness;
P2BR_eigenvector = centrality(P2BRgraph,'eigenvector'); %Eigenvector is the
relative score assigned to the value of a node
P2BRgraph.Nodes.eigenvector = P2BR_eigenvector;
P2BRgraph.Nodes
P2BRshortestpaths=distances(P2BRgraph) %Shortest paths is the shortest path
possible between any 2 nodes.
%This section repeats network map generation for non-survivor data.
P3B=NaN(22);
for j=1:22
 for i=1:22
 if P3(i,j)<=0.05/K;
 P3B(i,j)=R3(i,j);
 else P3B(i,j)=0;
 end
 end
end
P3Bgraph=graph(P3B)
P3BAbsgraph=abs(P3B);
P3BRgraph=graph(P3BAbsgraph)
Weight=P3BAbsgraph;
LWidths=5*P3BRgraph.Edges.Weight/max(P3BRgraph.Edges.Weight);
node_names = {'CRP','IL-6','TNFa','Pugh','MELD','HR','SDNN','SD2','Glucose','Urea','Creatinine','Na','Ammo
nia','Tot Prot','Albumin','Tot
Bili','PT_pc','INR','Hb','Plt','ALT','Indolo'};
P3BRgraph = graph(P3BAbsgraph,node_names);
subplot (1,2,2)
55
plot(P3BRgraph,'LineWidth',LWidths,'Layout','force','UseGravity',true,'Node
Color','red','EdgeColor','black','MarkerSize',4)
%"'EdgeLabel',P3BRgraph.Edges.Weight," - option to label edges with weight
in the plot function code.
title ('Non survivor 6M')
P3BRgraph.Nodes;
P3BR_degree = centrality(P3BRgraph,'Degree'); %Degree is the number of
edges connecting to each node.d
P3BRgraph.Nodes.degree = P3BR_degree;
P3BR_betweenness = centrality(P3BRgraph,'Betweenness'); %Betweenness
measures how often each graph node appears on a shortest path between two
nodes in the graph.
P3BRgraph.Nodes.betweenness = P3BR_betweenness;
P3BR_closeness = centrality(P3BRgraph,'closeness');%Closeness is the
inverse sum of the distance from a node to all other nodes in the graph
P3BRgraph.Nodes.closeness = P3BR_closeness;
P3BR_eigenvector = centrality(P3BRgraph,'eigenvector');%Eigenvector is the
relative score assigned to the value of a node
P3BRgraph.Nodes.eigenvector = P3BR_eigenvector;
P3BRgraph.Nodes
P3BRshortestpaths=distances(P3BRgraph)
%The outputs of this code include:
%1. Network map with 2 subplots comparing survivors and non-survivors of
the indicated time period.
%2. Graph properties with number of edges and nodes.
%3. Network measure table with common measures used to quantify networks.
%4. Matrix containing shortest paths between nodes.
