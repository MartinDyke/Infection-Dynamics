function [  Infected, Dead ] = timedinfection(k, n, p, q, t, r)
%This function generates a network using the Small_World_Network_Generator
%function and then simulates the spread of an infection through the
%network until t units of time have passed. The function then outputs the
%proportion of the network infected but not 'dead', and the proportion of
%the network 'dead'.


%Inputs are
%   n  = number of nodes
%   2k = average degree of the nodes
%   p  = rewiring probability
%   q  = infection probability
%   t  = number of units of time to run infection for
%	r  = number of iterations

%Outputs are
%   Infected = proportion of the network infected, but not
%                         'dead' at time t
%   Dead     = proportion of the network 'dead' at time t

y=1;
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
    %s is used as the timer
    
    infected=zeros(1,n);
    dead=zeros(1,n);
    counter=zeros(1,n);
    ran=randi(n,1);
    infected(ran)=1;
    counter(ran)=1;
    s=0;
    while s<t
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
        s=s+1;
    end
    Q(y)=sum(infected);
    T(y)=sum(dead);
    y=y+1;
end
Dead=sum(T)/(r*n);
Infected=sum(Q)/(r*n);
end

