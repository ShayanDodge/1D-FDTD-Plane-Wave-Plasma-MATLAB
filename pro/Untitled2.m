n=1
for n=1:10000
dn=dt*(k_5ea.*n_e.*n_Ar-k_6ea.*n_Ari.*n_e);
n_e=n_e+dn;
n_Ari=n_Ari+dt.*(k_5ea.*n_e.*n_Ar-k_6ea.*n_Ari.*n_e);
n_Ar=n_Ar-dt.*(k_5ea.*n_e.*n_Ar-k_6ea.*n_Ari.*n_e);

A(floor(n))=n_e(1);
h=figure(1)
h=plot(zcoorp,n_e)
h_2=figure(2)
h_2=plot(A)
end