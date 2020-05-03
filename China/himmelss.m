function ss = himmelss(k,data)

time = data.ydata(:,1);
Iqdata = data.ydata(:,2);
R1data = data.ydata(:,3);

Iq0 = Iqdata(1); E0 = k(13); A0 = k(14); 
I0 = k(15); R10 = R1data(1); Eq0 = k(16); R20 = k(17); Aq0= k(18);
y0 = [Iq0;E0;A0;I0;R10;Eq0;R20;Aq0];

[t,y] = ode45(@himmelode,time,y0,[],k);
Iqmodel = y(:,1);
R1model = y(:,5);

ss = sum((Iqdata-Iqmodel).^2+(R1data-R1model).^2);