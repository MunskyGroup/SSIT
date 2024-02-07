%%
function P = parseFspVector(states,fsp)
P = zeros(1,size(states,2));
for i = 1:size(states,2)
    P(i) = fsp(states(1,i)+1,states(2,i)+1,states(3,i)+1);
end
end
