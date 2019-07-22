%Code for Mutual Information network analysis
%This section clears all existing entities in the workspace, command window
lines and open menus and windows to prevent potential overlapping of data.
clc
close all
clear all
%This section loads in the data files: DatasetName =
importdata('filename.csv')
%This section generates Pearson's Correlation Coefficient symmetry
correlation matrices with corresponding p-values.
tic
th=1.05; %This is the selected significance threshold for the Mutual
Information analysis.
Surv6M = importdata ('22_6MData.csv');
A=Surv6M;
A=A'; %Information theory measures read code row-wise, hence the matrix
must be transformed.
toc
%This section filters and removes non-significant correlations from the
symmetry matrix and generates the network graph of the survivor and nonsurvivor group for a specific time period. The code shown is generating the
network graphs for the time period of 6 months.
X1=zeros(22);
for j=1:22
 for i=j:21
 X1(i+1,j)=kernelmi(A(j,:),A(i+1,:)); %kernelmi is an external mutual
information script. See methods for source.
 end
end
X2=zeros(22);
for p=1:22
 for q=p:22
56
 X2(p,q)=X1(q,p);
 end
end
X=X1+X2;
P4B=zeros(22);
for c=1:22;
 for d=1:22;
 if X(c,d)>=th;
 P4B(c,d)=X(c,d);
 else P4B(c,d)=0;
 end
 end
end
P4Bgraph=graph(P4B);
P4BAbsgraph=abs(P4B);
P4BRgraph=graph(P4BAbsgraph)
Weight=P4BAbsgraph;
LWidths=1*P4BRgraph.Edges.Weight/max(P4BRgraph.Edges.Weight);
node_names = {'CRP','IL-6','TNFa','Pugh','MELD','HR','SDNN','SD2','Glucose','Urea','Creatinine','Na','Ammo
nia','Tot Prot','Albumin','Tot
Bili','PT_pc','INR','Hb','Plt','ALT','Indolo'};
P4BRgraph = graph(P4BAbsgraph,node_names);
subplot (1,2,1)
plot(P4BRgraph,'LineWidth',LWidths,'Layout','force','UseGravity',true,'Node
Color','green','EdgeColor','black','MarkerSize',4)%'EdgeLabel',P2BRgraph.Ed
ges.Weight,
title ('Survivor 6M')
P4BRgraph.Nodes;
P4BR_degree = centrality(P4BRgraph,'Degree'); %Degree is the number of
edges connecting to each node.
P4BRgraph.Nodes.degree = P4BR_degree;
P4BR_betweenness = centrality(P4BRgraph,'Betweenness'); %Betweenness
measures how often each graph node appears on a shortest path between two
nodes in the graph.
P4BRgraph.Nodes.betweenness = P4BR_betweenness;
P4BR_closeness = centrality(P4BRgraph,'closeness');%Closeness is the
inverse sum of the distance from a node to all other nodes in the graph
P4BRgraph.Nodes.closeness = P4BR_closeness;
P4BR_eigenvector = centrality(P4BRgraph,'eigenvector');%Eigenvector is the
relative score assigned to the value of a node
P4BRgraph.Nodes.eigenvector = P4BR_eigenvector;
P4BRgraph.Nodes
P4BRshortestpaths=distances(P4BRgraph)
%This section repeats network map generation for non-survivor data.
NonSurv6M = importdata ('22_6MNData.csv');
B=NonSurv6M;
B=B';
X3=zeros(22);
for e=1:22
 for f=e:21
 X3(f+1,e)=kernelmi(B(e,:),B(f+1,:));
 end
end
X4=zeros(22);
for r=1:22
 for s=r:22
 X4(r,s)=X3(s,r);
 end
end
XT=X3+X4;
57
P5B=zeros(22);
for m=1:22;
 for n=1:22;
 if XT(m,n)>=th;
 P5B(m,n)=XT(m,n);
 else P5B(m,n)=0;
 end
 end
end
P5Bgraph=graph(P5B)
P5BAbsgraph=abs(P5B);
P5BRgraph=graph(P5BAbsgraph)
Weight=P5BAbsgraph;
LWidths=5*P5BRgraph.Edges.Weight/max(P5BRgraph.Edges.Weight);
node_names = {'CRP','IL-6','TNFa','Pugh','MELD','HR','SDNN','SD2','Glucose','Urea','Creatinine','Na','Ammo
nia','Tot Prot','Albumin','Tot
Bili','PT_pc','INR','Hb','Plt','ALT','Indolo'};
P5BRgraph = graph(P5BAbsgraph,node_names);
subplot (1,2,2)
plot(P5BRgraph,'LineWidth',LWidths,'Layout','force','UseGravity',true,'Node
Color','red','EdgeColor','black','MarkerSize',4)%'EdgeLabel',P2BRgraph.Edge
s.Weight,
title ('Non survivor 6M')
P5BRgraph.Nodes;
P5BR_degree = centrality(P5BRgraph,'Degree'); %Degree is the number of
edges connecting to each node.d
P5BRgraph.Nodes.degree = P5BR_degree;
P5BR_betweenness = centrality(P5BRgraph,'Betweenness'); %Betweenness
measures how often each graph node appears on a shortest path between two
nodes in the graph.
P5BRgraph.Nodes.betweenness = P5BR_betweenness;
P5BR_closeness = centrality(P5BRgraph,'closeness');%Closeness is the
inverse sum of the distance from a node to all other nodes in the graph
P5BRgraph.Nodes.closeness = P5BR_closeness;
P5BR_eigenvector = centrality(P5BRgraph,'eigenvector');%Eigenvector is the
relative score assigned to the value of a node
P5BRgraph.Nodes.eigenvector = P5BR_eigenvector;
P5BRgraph.Nodes
P5BRshortestpaths=distances(P5BRgraph)
%The outputs of this code include:
%1. Network map with 2 subplots comparing survivors and non-survivors of
the indicated time period.
%2. Graph properties with number of edges and nodes.
%3. Network measure table with common measures used to quantify networks.
%4. Matrix containing shortest paths between nodes.
