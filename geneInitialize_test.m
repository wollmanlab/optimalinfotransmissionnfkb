function delg = geneInitialize_test(t,g,TFt,TF,params0)
% Model for Down strem gene regulated by NFkB as a transcription factor
if t > params0(6)
   TF1 = interp1(TFt,TF,t-params0(6)); % interpolate the data set (TFt TF) at time t for g1
else
   TF1 = TF(1);
end
%
delg = params0(1) + ((params0(2)*((TF1)^params0(5)))/((params0(4)^params0(5))+((TF1)^params0(5))))...
         - (params0(3)* g(1)); 
end





