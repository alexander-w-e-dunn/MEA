function  [CIJtree] = mst_wu(CIJ)
%mst_wu        mimimum spanning tree
%   This script is based on backbone_wu but takes only the tree
%   adaption by AWE Dunn and D Akarca, Cambridge, 2020
%   Reference: Hidalgo et al. (2007) Science 317, 482.
%              Hagmann et al. (2008) PLoS Biol
%
%   Olaf Sporns, Indiana University, 2007/2008/2010/2012
%
%   [CIJtree] = mst_wu(CIJ)
%
%   The network backbone contains the dominant connections in the network
%   and may be used to aid network visualization. This function computes
%   the backbone of a given weighted and undirected connection matrix CIJ, 
%   using a minimum-spanning-tree based algorithm.
%
%   input:      CIJ,    connection/adjacency matrix (weighted, undirected)
%   output: 
%           CIJtree,    connection matrix of the minimum spanning tree of CIJ
%
%   NOTE: nodes with zero strength are discarded.

N = size(CIJ,1);
CIJtree = zeros(N);

% find strongest edge (note if multiple edges are tied, only use first one)
[i,j,s] = find(max(max(CIJ))==CIJ);                      %#ok<*ASGLU>
im = [i(1) i(2)];
jm = [j(1) j(2)];

% copy into tree graph
CIJtree(im,jm) = CIJ(im,jm);
in = im;
out = setdiff(1:N,in);

% repeat N-2 times
for n=1:N-2
    
    % find strongest link between 'in' and 'out',ignore tied ranks
    [i,j,s] = find(max(max(CIJ(in,out)))==CIJ(in,out)); 
    im = in(i(1));
    jm = out(j(1));
    
    % copy into tree graph
    CIJtree(im,jm) = CIJ(im,jm); CIJtree(jm,im) = CIJ(jm,im);
    in = [in jm];                                       %#ok<AGROW>
    out = setdiff(1:N,in);

end;



