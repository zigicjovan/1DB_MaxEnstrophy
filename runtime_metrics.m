
name = 'CGPR';
% name = 'RCGPR';
% name = 'CGRMIL';
% name = 'RCGRMIL';
ensvec = logspace(3,6,31);
timepts = 25;
E0levels = 18;

max = NaN(timepts,3*E0levels);

for i = 1:E0levels

        ens = ensvec(i);

        filemax = [ pwd '/runtime/runtime_2048_LuMax_t1_' name '_E0(' num2str(ens) ')_lambda(1).dat'];
        max( : , (i*3)-2 : i*3 ) = readmatrix(filemax);

end

%%% find average runtime for all E0 levels by timept for each gradient method
%%% matrix = [ timept , average runtime ]

avemax = zeros(timepts,2);
avemax(:,1) = max(:,1);
for i = 1: E0levels
    avemax(:,2) = avemax(:,2) + max(:,i*3)  ;
end
avemax(:,2) = avemax(:,2)/E0levels;
fileavemax = [ pwd '/runtime/runtimeave_2048_LuMax_t1_' name '_E0(' num2str(ens) ')_lambda(1).dat'];
writematrix(avemax,fileavemax,'Delimiter','tab')

%%% find total runtime (timept 25, col 3) for all E0 levels for each gradient method
%%% matrix = [ E0 , total runtime ]

totalmax = zeros(E0levels,2);
for i = 1: E0levels
    totalmax(i,1) = ensvec(i);
    totalmax(i,2) = max(25,i*3);
end
filetotalmax = [ pwd '/runtime/runtimetotal_2048_LuMax_t1_' name '_E0(' num2str(ens) ')_lambda(1).dat'];
writematrix(totalmax,filetotalmax,'Delimiter','tab')

