% 
% This function is to compare rich and poor nodes 
% in terms of betweenness centrality , degree, nodal efficiency,
% participation coefficient 
% 
% This idea is based on results in
% 
% Vértes, P. E., Alexander-Bloch, A. & Bullmore, E. T. 
% Philos. Trans. R. Soc. B Biol. Sci. 369, 20130531 (2014).
% 
% INPUTS: 
%   directory where you keep the file to be analysed
%   filename containing the variable called adjM which is a weighted
%   adjacency matrix
% 
% OUTPUTS:
%   pc is the mean participation coefficient for the rich nodes (column 1)
%   and poor nodes (column 2). bc, dk and ne are the same for betweenness
%   centrality, degree 'k', and LOCAL nodal efficiency 
%   number of nodes and edges in each network are also outputted
% 
% Author: AWE DUNN, JUNE 2020, U. of Cambridge
% 
function [pc, bc, dk, ne, num_nodes, num_edges] = richVSpoor_fcn(directory, filename)


load(filename);

if ~exist('channels')
    fprintf(2,'\n  WARNING: channel order not saved in spike matrix \n used default in line 50 of batch_getHeatMaps_fcn \n \n')
    channels = [47,48,46,45,38,37,28,36,27,17,26,16,35,25,15,14,24,34,13,23,12,22,33,21,32,31,44,43,41,42,52,51,53,54,61,62,71,63,72,82,73,83,64,74,84,85,75,65,86,76,87,77,66,78,67,68,55,56,58,57];
end

if ~exist('adjM') & contains(filename,'CTRL')
    adjM = adjM2s(:,:,1);
end

%% get RC and BC scores
count2 = 1;

RC_score = zeros(length(channels),1);
percentile_thresholds = [80] ; % i.e. vector of percentiles, 80 means top 20 % weights
for cutoff = percentile_thresholds
    
    threshold = prctile(adjM(:), cutoff);
    
    edges=adjM;
    edges = edges - eye(size(edges));
    edges(find(isnan(edges))) = 0;
    edges(find(edges < threshold)) = 0;
    edges(find(edges >= threshold))= 1;    
    
    if round(max(sum(edges))/5)*5>0 % check for enough nodes
        
                % rich club
        k               =   max(sum(edges));
        [Rc,Nk,Ek]      =   rich_club_bu(edges,k);
        RC              =   max(Rc);
        maxKrand        =   min(find(Rc==RC));
        RC_nodes_vec    =   find(sum(edges) >= maxKrand);
        Poor_nodes_vec  =   find(sum(edges) <  maxKrand);
        RC_node_IDs     =   channels(RC_nodes_vec);
        Poor_node_IDs   =   channels(Poor_nodes_vec);
        
        RC_score(RC_nodes_vec)      =    RC_score(RC_nodes_vec)+1;
        
        BC_vec(:,count2)          = betweenness_bin(edges);
        DegreeVec(:,count2)       = sum(edges);
        L_eff_Vec(:,count2)       = efficiency_bin(edges,1); % second input: 1 is local eff. 0 is global
            gamma = 1; % 1 is default; gamma>1, detects smaller modules; 0<=gamma<1 means large modules
            % schroeter et al method (blondel et al - see scriot for
            % reference
            [M,Q]=community_louvain(edges,gamma);
            % vertes et al method (newman method - see script for
            % reference)
%             [Ci,Q]=modularity_und(edges,gamma); % Q is modularity score, Ci is the vector of groups
        pcoef_Vec(:,count2)       = participation_coef(edges,M,0);%input 3: 0 means binary undirected graph

        pc = [  mean(pcoef_Vec(RC_nodes_vec))     mean(pcoef_Vec(Poor_nodes_vec)) ];
        bc = [  mean(BC_vec(RC_nodes_vec))        mean(BC_vec(Poor_nodes_vec)) ];
        dk = [  mean(DegreeVec(RC_nodes_vec))     mean(DegreeVec(Poor_nodes_vec)) ];
        ne = [  mean(L_eff_Vec(RC_nodes_vec))     mean(L_eff_Vec(Poor_nodes_vec)) ];
        num_nodes = length(find(sum(edges) >=1));
        num_edges = sum(sum(edges)) / 2; % only counts each connection once
        
    end
    
    count2 = count2 + 1;
    
end

end



