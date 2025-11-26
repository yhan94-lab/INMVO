function [Best_universe_Inflation_rate,Best_universe,Convergence_curve]=INMVO(N,MaxFEs,lb,ub,dim,fobj)
fitcount=0;
Best_universe=zeros(1,dim);
Best_universe_Inflation_rate=inf;
Universes=initialization(N,dim,ub,lb);
WEP_Max=1;
WEP_Min=0.2;
FEs=0;
Convergence_curve=[];
while FEs<MaxFEs
    WEP=WEP_Min+FEs*((WEP_Max-WEP_Min)/MaxFEs);
    TDR=1-((FEs)^(1/6)/(MaxFEs)^(1/6));
    Inflation_rates=zeros(1,size(Universes,1));
    for i=1:size(Universes,1)
        Flag4ub=Universes(i,:)>ub;
        Flag4lb=Universes(i,:)<lb;
        Universes(i,:)=(Universes(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        Inflation_rates(1,i)=fobj(Universes(i,:));
        FEs=FEs+1;
        if Inflation_rates(1,i)<Best_universe_Inflation_rate
            Best_universe_Inflation_rate=Inflation_rates(1,i);
            Best_universe=Universes(i,:);
        end
        fitcount=fitcount+1;
        Convergence_curve(1,fitcount)= Best_universe_Inflation_rate;
    end
    [sorted_Inflation_rates,sorted_indexes]=sort(Inflation_rates);
    for newindex=1:N
        Sorted_universes(newindex,:)=Universes(sorted_indexes(newindex),:);
    end
    normalized_sorted_Inflation_rates=normr(sorted_Inflation_rates);
    Universes(1,:)= Sorted_universes(1,:);
    for i=2:size(Universes,1)
        Back_hole_index=i;
        for j=1:size(Universes,2)

            r1 = sin(0.7 * pi / rand);
            if r1<normalized_sorted_Inflation_rates(i)
                White_hole_index=RouletteWheelSelection(-sorted_Inflation_rates);
                if White_hole_index==-1
                    White_hole_index=1;
                end
                Universes(Back_hole_index,j)=Sorted_universes(White_hole_index,j);
            end
            if (size(lb,2)==1)
                r2 = sin(0.7 * pi / rand);
                if r2<WEP
                    r3 = sin(0.7 * pi / rand);
                    if r3<0.5
                        Universes(i,j)=Best_universe(1,j)+TDR*((ub-lb)*sin(0.7 * pi / rand)+lb);
                    end
                    if r3>0.5
                        Universes(i,j)=Best_universe(1,j)-TDR*((ub-lb)*sin(0.7 * pi / rand)+lb);
                    end
                end
            end
            if (size(lb,2)~=1)
                r2 = sin(0.7 * pi / rand);
                if r2<WEP
                    r3 = sin(0.7 * pi / rand);
                    if r3<0.5
                        Universes(i,j)=Best_universe(1,j)+TDR*((ub(j)-lb(j))*sin(0.7 * pi / rand)+lb(j));
                    end
                    if r3>0.5
                        Universes(i,j)=Best_universe(1,j)-TDR*((ub(j)-lb(j))*sin(0.7 * pi / rand)+lb(j));
                    end
                end
            end

        end
    end
    options = optimset('MaxFunEvals', floor(MaxFEs * 0.1));
    [pos, fval, ~, output]  = fminsearchbnd(fobj,Best_universe,lb,ub,options);
    if fval < Best_universe_Inflation_rate
        Best_universe_Inflation_rate = fval;
        Best_universe = pos;
    end
    FEs = FEs + output.funcCount;

    Convergence_curve(1,fitcount+1:fitcount+output.funcCount)= Best_universe_Inflation_rate;
    fitcount=fitcount+output.funcCount;
end
Convergence_curve=Convergence_curve(:,(fitcount-MaxFEs+1):end);
end
