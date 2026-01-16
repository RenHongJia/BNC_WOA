function [it_gb_k2score,DAG] =Learn_bnc_woa_struct(data,bnetinfo,nPop,MaxIt,params)



max_fan_in= params.max_fan_in;

nNodes=size(data,1);
attNodes = nNodes-1; 
nAssignments = bnetinfo.node_sizes;


it_gb_k2score=zeros(1,MaxIt);

firefly.Position = [];
firefly.Cost = [];

pop = repmat(firefly, nPop, 1);
GlobalBest.Cost = -Inf;

for i = 1:nPop  

    Position_temp = mk_rnd_dag(attNodes,max_fan_in);

     Position_temp = transform_BNC(Position_temp,nNodes);

    pop(i).Position = Position_temp(:)'; 

    pop(i).Cost = score_ren(data,pop(i).Position,nAssignments);
    
end



for idx=1:length(pop)
    if pop(idx).Cost > GlobalBest.Cost
        GlobalBest.Cost=pop(idx).Cost;
        GlobalBest=pop(idx);
    end
end

  
 for it = 1:MaxIt
   
    newpop=repmat(firefly,nPop,1);
    a = 1 - (1 / (1 + exp(-10 * ((it/MaxIt) - 0.5))));

    for i = 1:nPop
        A=a;

        l=rand();
        p=rand();
        
        newpop(i).Cost = -inf;
        new_pos=zeros(1,length(pop(i).Position));
 
        rnum1=rand();
       
        if p<0.5  
            if rnum1<A   
                rand_leader_index = randi([1,nPop]);
                X_rand = pop(rand_leader_index).Position;
                
                diff=X_rand-pop(i).Position;   

                C=rand(1,length(pop(i).Position));
              
                val=0.5;     
                for ii=1:length(pop(i).Position)
                    if C(ii) > val && diff(ii)~=0
                        new_pos(ii)=~X_rand(ii);
                    else
                        new_pos(ii)=X_rand(ii);
                    end
                end
                
                newpop(i).Position =  new_pos;

                newpop(i).Position =check_solution(newpop(i).Position);
                newpop(i).Position = illegal_structure_resive(newpop(i).Position,nNodes);
 
                newpop(i).Cost = score_ren(data,newpop(i).Position,nAssignments);
            
            else
               
                
                diff=GlobalBest.Position-pop(i).Position;  
                bestPos=GlobalBest.Position;

                C=rand(1,length(pop(i).Position));
                
                val=0.5;  
                for ii=1:length(pop(i).Position)
                    if C(ii) > val && diff(ii)~=0
                        new_pos(ii)=~bestPos(ii);
                    else
                        new_pos(ii)=bestPos(ii);

                    end
                end   
                newpop(i).Position = new_pos;

                newpop(i).Position =check_solution(newpop(i).Position);
                newpop(i).Position = illegal_structure_resive(newpop(i).Position,nNodes);

                newpop(i).Cost = score_ren(data,newpop(i).Position,nAssignments);

            end
        else

            D=GlobalBest.Position-pop(i).Position;
            bestPos=GlobalBest.Position;
            Hanming_dis=sum(D~=0);
            if Hanming_dis==0
                newpop(i)= GlobalBest;
                continue;
            end
         
            l=0.5*(it/MaxIt);
            k=0.6*(exp(l)).*cos(pi*l);
            beta=Hanming_dis*k;
            if beta<=1
                beta=1;
                
            elseif beta >1        
                beta=round(beta); 
            end
            if beta==1
                beta1=0;
                beta2=beta-beta1;
            end
            if beta>=2
                beta1=ceil(beta/2);
                beta2=beta-beta1;
            end
            mutation_pos=randperm(length(pop(i).Position));
            bt1=0;
            bt2=0;
            jj=1;
            while(bt1~=beta1)
                if mod(mutation_pos(jj),nNodes)==0 
                    jj=jj+1;
                    continue;
                end
                if D(mutation_pos(jj))~=0     
                    bestPos(mutation_pos(jj))=~bestPos(mutation_pos(jj));
                    bt1=bt1+1;
                    jj=jj+1;
                    continue;
                end
                jj=jj+1;
            end
            jj=1;
             while(bt2~=beta2)
                if mod(mutation_pos(jj),nNodes)==0 
                    jj=jj+1;
                    continue;
                end
                if D(mutation_pos(jj))==0 || beta==1    
                    bestPos(mutation_pos(jj))=~bestPos(mutation_pos(jj));
                    bt2=bt2+1;
                    jj=jj+1;
                    continue;
                end
                jj=jj+1;
             end
            if bt1+bt2~=beta
                error('    ');
            end

            newpop(i).Position = bestPos;

            newpop(i).Position =check_solution(newpop(i).Position);
            newpop(i).Position = illegal_structure_resive(newpop(i).Position,nNodes);

            newpop(i).Cost = score_ren(data,newpop(i).Position,nAssignments);
 
        end
    end

    pop=newpop;
    for idx=1:length(pop)
        if pop(idx).Cost > GlobalBest.Cost
            GlobalBest.Cost=pop(idx).Cost;
            GlobalBest=pop(idx);
        end
    end

    it_gb_k2score(it)=GlobalBest.Cost;
    if  mod(it,20)==0
        disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(GlobalBest.Cost)]);
        
    end    
 end

DAG = GlobalBest.Position;
DAG = reshape(DAG,nNodes,nNodes);

end

