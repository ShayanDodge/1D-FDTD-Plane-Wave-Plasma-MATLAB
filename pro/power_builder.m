power_total=zeros(1,28);
for i=10:10:1000
load(['C:\Users\Shayan-pc\Desktop\10KHZ\Jet_' num2str(i) '.mat'])
power_total=power_total+P_total;
end
power_total=power_total/100;
power_diss=zeros(1,28);
for i=10:10:1000
load(['C:\Users\Shayan-pc\Desktop\10KHZ\Jet_' num2str(i) '.mat'])
power_diss=power_diss+p_Diss;
end
power_diss=power_diss/100;
temp=zeros(1,28);
for i=10:10:1000
load(['C:\Users\Shayan-pc\Desktop\10KHZ\Jet_' num2str(i) '.mat'])
temp=temp+T_e_ev;
end
temp=temp/100;