function [Infected, Dead, timetaken] = proportioninfection(k, n, p, q, e, r)
%This function generates a network using the Small_World_Network_Generator
%function and then simulates the spread of an infection through the
%network until a target proportion of the network has been infected. This
%function is intended to measure what effect the structure of the network
%has on the time taken for the infection to spread to the target
%proportion. If there are no more infected individuals which are not 'dead'
%then the algorithm redoes this iteration

%Inputs are
%   n  = number of nodes
%   2k = average degree of the nodes
%   p  = rewiring probability
%   q  = infection probability
%   e  = target proportion. Each iteration terminates once this proportion
%        of the network has been infected, or if there are no more infected
%        individuals
%	r  = number of iterations

%Outputs are
%   Infected  = Average proportion of the network infected at end of
%               iterations
%   Dead      = Average proportion of the network 'dead' at end of all 
%               iterations
%   timetaken = Average number of units of time taken to reach target
%               proportion or until there are no more infected individuals

%y is used as the counter to count the number of iterations completed. R, Q
%and T are used to store the proportions of infected, 'dead' and time
%taken for each iteration

y=1;
R=zeros(1,r);
Q=zeros(1,r);
T=zeros(1,r);
while y<=r
    
    %For each iteration we generate a new network using the
    %Small_World_Network_Generator algorithm
    
    [~,~, A] = Small_World_Network_Generator( k, n, p, 1);
    
    %We define the vectors infected, dead and counter. infected(i)=1 of 
    %node i is infected but not yet dead, infected(i)=0 if node i is not 
    %yet infected, or is 'dead'. dead(i)=1 if node i has 'died'. 
    %counter(i) counts how long the node i has been infected for.
    
    infected=zeros(1,n);
    dead=zeros(1,n);
    counter=zeros(1,n);
    
    %We choose an integer ran between 1 and the number of nodes randomly 
    %and define node ran to be infected. We also set the time to 0, and 
    %proportion_affected, %which will be the proportion of the network 
    %either infected or 'dead' to 0
    
    ran=randi(n,1);
    infected(ran)=1;
    counter(ran)=1;
    time=0;
    proportion_affected=0;
    
    while proportion_affected<e && sum(infected)~=0
        
        for i=1:n
            
            if infected(i)==1
                
                for j=1:n
                    
                    %The if loop checks for each j if nodes i and j have 
                    %an edge between them. If there is, and node j is 
                    %neither infected or 'dead' then we infect node j with 
                    %probability q
                    
                    if A(i,j)==1 && infected(j)==0 && dead(j)==0
                        ran=rand(1);
                        
                        if ran<=q
                            infected(j)=1;
                            counter(j)=1;
                        end
                    end
                end
                
                %We generate a new random number between 0 and 1, and use 
                %this to determine whether the infected node i 'dies' at 
                %this time. As counter(i) increases, the value of 
                %1/(exp(counter(i)/30)) decreases, and so it is more 
                %likely that the individual will 'die'
                
                ran=rand(1);
                if ran > 1/(exp(counter(i)/30))
                    infected(i)=0;
                    dead(i)=1;
                else
                    counter(i)=counter(i)+1;
                end
            end
        end
        
        proportion_affected=(sum(dead)+sum(infected))/(n);
        time=time+1;
        
    end
    
    R(y)=time;
    Q(y)=sum(infected);
    T(y)=sum(dead);
    
    %If the number of infected individuals is now zero, implying that the
    %target proportion was not reached, then y (the iteration counter) 
    %does not increase, and so the algorithm runs this iteration again
    
    if sum(infected)~=0
        y=y+1;
    end
    
end

Dead=sum(T)/(r*n);
Infected=sum(Q)/(r*n);
timetaken=sum(R)/(r);

end

