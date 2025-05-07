function isbetter = IsRbetter(ind1,ind2,M)
isbetter = true;
for i = 1 : M
    if ind1 > ind2
        isbetter = false;
        break;
    end
end