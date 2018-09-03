function [Dead] = quarantine(k, n, p, q, qu, r)
%Quarantine: This algorithm models the effect that quarantining individuals
%has on the spread of an infection through a network. The algorithm runs
%until all nodes are classed as 'dead', i.e. there are no more infected
%nodes. The algorithm outputs the average proportion of the network that
%has 'died'.

%Inputs are
%   n  = number of nodes
%   2k = average degree of the nodes
%   p  = rewiring probability
%   q  = infection probability
%   qu = number of units of time after an individual becomes infected that
%        they are quarantined
%	r  = number of iterations



%y is the counter for the repetition, T is the storage of the values for
%each repeat of the proportion of the network eventually infected.

y=1;
T=zeros(1,r);
while y<=r
    
    %For each repetition we create a new adjacency matrix using the
    %Small_World_Network_Generator algorithm
    
    [~,~,adjacency_matrix]=Small_World_Network_Generator(k, n, p,1);
    
    %We define the vectors infected, dead and counter. infected(i)=1 of 
    %node i is infected but not yet dead, infected(i)=0 if node i is not 
    %yet infected, or is 'dead'. dead(i)=1 if node i has 'died'. counter(i) 
    %counts how long the node i has been infected for.
    
    
    infected=zeros(1,n);
    dead=zeros(1,n);
    counter=zeros(1,n);
    
    %We choose an integer ran between 1 and the number of nodes randomly 
    %and define node ran to be infected
    
    ran=randi(n,1);
    infected(ran)=1;
    counter(ran)=1;
    
    %while loop continues until no more individuals are infected
    
    while sum(infected)~=0
        for i=1:n
            if infected(i)==1
                
                %sm is used as a counter. If node i has infected all its 
                %neighbours, then we may as well class it as 'dead'. 
                %This means that the algorithm does not then keep running 
                %the algorithm for node i even though it cannot possibly 
                %infect other nodes
                
                sm=0;
                for j=1:n
                    
                    %The if loop checks for each j if nodes i and j have 
                    %an edge between them. If there is, and node j is 
                    %neither infected or 'dead' then we infect node j with 
                    %probability q
                    
                    if adjacency_matrix(i,j)==1 && infected(j)==0 && dead(j)==0
                        sm=sm+1;
                        ran=rand(1);
                        if ran<=q
                            infected(j)=1;
                            counter(j)=1;
                        end
                    end
                end
                
                %If node i has been infected for the length of time 
                %qu then it is classed as dead. Else, if it has 
                %infected all connected nodes it is classed dead. 
                %Otherwise we add one to the counter
                
                if counter(i)==qu || sm==n
                    infected(i)=0;
                    dead(i)=1;
                end
                counter(i)=counter(i)+1;
            end
        end
    end
    T(y)=sum(dead);
    y=y+1;
    
end
Dead=sum(T)/(r*n);
end

