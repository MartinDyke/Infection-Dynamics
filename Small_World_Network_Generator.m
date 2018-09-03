function [C, L, A] = Small_World_Network_Generator( k, n, p, r)
%This function generates a regular lattice with n individuals each
%connected to their nearest 2k neighbours. It then randomly rewires each
%edge with probability p under the condition that the graph remains
%connected, outputting the clustering coefficient, average path length and
%adjacency matrix for the network.

%Inputs are
%   n  = number of nodes
%   2k = average degree of the nodes
%   p  = rewiring probability
%	r  = number of iterations

%Outputs are
%   C = the average clustering coefficient over all iterations
%   L = the average path length over all iterations
%   A = the adjacency matrix, used for when other functions call this
%       function

%G is used to store the average path lengths and clustering coefficients
%for the r repetitions, and y is used as a counter for the repetitions

G=zeros(r,2);
y=1;

while y<=r
    
    %We define A to be the adjacency matrix for the regular lattice. 
    %We store it as a sparse matrix since most of the entries are zero, 
    %and since later on the function graphconncomp requires the input to 
    %be a sparse matrix.
    
A = sparse(gallery('circul', [ 0 ones(1,k) zeros(1,n - 2*k-1) ones(1,k)]));
    
    for i = 1:n
        
        
        %This for and if loop chooses the k nearest neighbours which are 
        %numbered higher than node i. We only need to run for the k larger 
        %neighbors since the graph is undirected. 
        %Hence, for a 10 node network, each node connected to its 4 nearest
        %neighbours, for node 1 it possibly rewires the links to nodes 2 
        %and 3. The links to node 10 will be covered when the for loop 
        %reaches 10, as it will rewire the link from 10 to 1, since 1 is 
        %'larger' than 10. This means we must use mod in order to make sure
        %it determines node 10+1 to be node 1, as opposed to the 
        %nonexistent node 11
        
        for m=1:k
            j= mod(i+m,n);
            if j==0
                j=n;
            end
            w=rand(1);
            
            %This if loop is for randomly rewiring edges. So if the 
            %random number w between 0 and 1 is less than p the if loop is 
            %executed.
            
            if w <=p
                q1=1;
                q2=1;
                S=0;
                while S ~=1
                    A(q1,q2) = 0;
                    A(q2,q1) = 0;
                    
                    %q1 and q2 are our randomly chosen nodes.
                    
                    q1=randi(n,1);
                    q2=randi(n,1);
                    
                    
                    %The algorithm then checks that q1 is not q2, and
                    %that there is currently no edge between q1 and q2, and
                    %keeps generating new q1 and q2 until this is the case.
                    
                    while A(q1,q2)==1 || q1==q2
                        q1=randi(n,1);
                        q2=randi(n,1);
                    end
                    
                    %The edges between nodes i and j are then deleted, and 
                    %edges between q1 and q2 are added in.
                    
                    A(i,j) = 0;
                    A(j,i) = 0;
                    A(q1,q2) = 1;
                    A(q2,q1) = 1;
                    
                    %The algorithm then checks to see if the graph is still
                    %connected. If S=1 then the graph is connected and the 
                    %while loop is terminated. If S~=1 then the graph has 
                    %multiple components. In this case the algorithm begins 
                    %the while loop again, deletes the edges just added in,
                    %and generates new edges. 
                    %It repeats this until the graph remains connected.
                    %We have to define S=0 before every random rewiring 
                    %since otherwise S=1 would still be 0 from the previous 
                    %rewiring and the while loop would not run.
                    %We also have to set q1=q2=1 before every rewiring. 
                    %Otherwise the algorithm would immediately delete the 
                    %last rewiring. Since there cannot be an edge between 
                    %node 1 and itself it is no issue if it attempts to 
                    %delete this nonexistent edge
                    
                    [S,~]=graphconncomp(A);
                end
                
            end
        end
        
    end
    
    %Calculate the clustering coefficient for the network
    %First we need to calculate the clustering coefficient for each node
    
    c=zeros(1,n);
    for i=1:n
        
        %Q returns the indexes of all the vertices connected to node i, 
        %which l then calculates the length of to give the number of 
        %vertices attached to node i
        
        Q=find(A(:,i));
        l=length(Q);
        
        %the if loop is since the degree of a vertex must be at least 2,
        %otherwise there cannot be any triangles
        
        if l >= 2
            
            %L is the adjacency matrix between the vertices attached to 
            %node i. We can then sum all the edges, each 2 of which would 
            %represent a triangle, and divide by the total possible number 
            %of edges
            
            K =A(Q,Q);
            c(i)=sum(K(:))/(l*(l-1));
        end
    end
    clustering_coefficient = sum(c)/n; %average clustering coefficient
    
    %shortest distance between every two vertices
    
    [dist] = graphallshortestpaths(A); 
    
    average_path_length=sum(dist(:))/(n*(n-1));
    G(y,1)=clustering_coefficient;
    G(y,2)=average_path_length;
    y=y+1;
end
C=sum(G(:,1))/r;
L=sum(G(:,2))/r;
end











