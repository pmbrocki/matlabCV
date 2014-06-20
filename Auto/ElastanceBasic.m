function Elh = ElastanceBasic(t,T,Trf,Ed,Es,Tsf)

Tr  = Trf*T;
Ts  = Tsf*T;

for i = 1:length(t)
    if t(i)<=Ts
        Elh(i) = Ed + ((Es-Ed)./2).*(1-cos(pi.*t(i)./Ts(i)));
    elseif (t(i)<=Ts+Tr) 
        Elh(i) = Ed + ((Es-Ed)./2).*(cos((pi./Tr(i)).*(t(i)-Ts(i))) + 1);
    else
        Elh(i)=Ed;
    end
end
