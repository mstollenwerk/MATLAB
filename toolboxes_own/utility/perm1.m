function perms = perm1(perm_,number)

p = size(perm_,2);
perms = NaN(p,p);

rest_perm = perm_(perm_ ~= number);
perms(1,:) = [number, rest_perm];
perms(end,:) = [rest_perm, number];
for jj = 2:(p-1)
    perms(jj,:) = [rest_perm(1:jj-1) number rest_perm(jj:end)];
end

end

