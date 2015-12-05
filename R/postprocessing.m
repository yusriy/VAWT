%NOTE: THE USER NEED TO CHANGE THE FILE NAME ACCORDINGLY

tic;    %Start run timer
u0 = rawdata(:,1);  %'rawdata' is the txt file imported to Matlab
v0 = rawdata(:,2);
w0 = rawdata(:,3);
t0 = rawdata(:,4) + 273.15; %in the units Kelvin (K)

avgu0 = mean(u0(:)); %averages of u0, v0, w0, and t0 for record
avgv0 = mean(v0(:));
avgw0 = mean(w0(:));
avgt0 = mean(t0(:));

stdu0 = std(u0(:)); %standard deviations of u0, v0, w0, and t0 for record
stdv0 = std(v0(:));
stdw0 = std(w0(:));
stdt0 = std(t0(:));

upperboun_u0 = avgu0 + 3.5*stdu0;   %upper and lower boundary of u0, v0, w0, and t0 for record
lowerboun_u0 = avgu0 - 3.5*stdu0;

upperboun_v0 = avgv0 + 3.5*stdv0;
lowerboun_v0 = avgv0 - 3.5*stdv0;

upperboun_w0 = avgw0 + 3.5*stdw0;
lowerboun_w0 = avgw0 - 3.5*stdw0;

upperboun_t0 = avgt0 + 3.5*stdt0;
lowerboun_t0 = avgt0 - 3.5*stdt0;

%to remove spikes of 3.5 stdev from the mean, positive and negative
u01 = u0(u0 < upperboun_u0 & u0 > lowerboun_u0 & v0 < upperboun_v0 & v0 > lowerboun_v0 & w0 < upperboun_w0 & w0 > lowerboun_w0 & t0 < upperboun_t0 & t0 > lowerboun_t0);
v01 = v0(u0 < upperboun_u0 & u0 > lowerboun_u0 & v0 < upperboun_v0 & v0 > lowerboun_v0 & w0 < upperboun_w0 & w0 > lowerboun_w0 & t0 < upperboun_t0 & t0 > lowerboun_t0);
w01 = w0(u0 < upperboun_u0 & u0 > lowerboun_u0 & v0 < upperboun_v0 & v0 > lowerboun_v0 & w0 < upperboun_w0 & w0 > lowerboun_w0 & t0 < upperboun_t0 & t0 > lowerboun_t0);
t01 = t0(u0 < upperboun_u0 & u0 > lowerboun_u0 & v0 < upperboun_v0 & v0 > lowerboun_v0 & w0 < upperboun_w0 & w0 > lowerboun_w0 & t0 < upperboun_t0 & t0 > lowerboun_t0);

%To determine number of data in record
record = length(u01); %An abitrary choice of u01 is chosen

%To determine the number of datasets in the record for different block
%averages
n5 = record / 3000;
n10 = record / 6000;
n15 = record / 9000;
n30 = record / 18000;
n60 = record / 36000;
n120 = record / 72000;

%Cutting the data into 5 min sets
n = 0;

% Compute size and allocate data, added by Eric Johnson.  Need to apply to
% other variables.

numCols = n5;
numRows = 1 + n:3000 + n;
numRows = max(numRows) - n;

u05 = zeros(numRows,numCols);
v05 = zeros(numRows,numCols);
w05 = zeros(numRows,numCols);
t05 = zeros(numRows,numCols);

% The K loop can be a PARFOR loop

for k = 1:n5     

for i = 1 + n:3000 + n  %3000 data points for 5 mins
    
    u05(i - n,k) = u01(i);
    v05(i - n,k) = v01(i);
    w05(i - n,k) = w01(i);
    t05(i - n,k) = t01(i);
    
end

n = i;

end

%Cutting the data into 10 min sets
n =0;
% Compute size and allocate data, added by Eric Johnson.  Need to apply to
% other variables.

numCols = n10;
numRows = 1 + n:6000 + n;
numRows = max(numRows) - n;

u010 = zeros(numRows,numCols);
v010 = zeros(numRows,numCols);
w010 = zeros(numRows,numCols);
t010 = zeros(numRows,numCols);

for k = 1:n10     

for i = 1 + n:6000 + n  %6000 data points for 10 mins
    
    u010(i - n,k) = u01(i);
    v010(i - n,k) = v01(i);
    w010(i - n,k) = w01(i);
    t010(i - n,k) = t01(i);
    
end

n = i;

end

%Cutting the data into 15 min sets

n = 0;

% Compute size and allocate data, added by Eric Johnson.  Need to apply to
% other variables.

numCols = n15;
numRows = 1 + n:9000 + n;
numRows = max(numRows) - n;

u015 = zeros(numRows,numCols);
v015 = zeros(numRows,numCols);
w015 = zeros(numRows,numCols);
t015 = zeros(numRows,numCols);

for k = 1:n15     

for i = 1 + n:9000 + n  %9000 data points for 15 mins
    
    u015(i - n,k) = u01(i);
    v015(i - n,k) = v01(i);
    w015(i - n,k) = w01(i);
    t015(i - n,k) = t01(i);
    
end

n = i;

end


%Cutting the data into 30 min sets

n = 0;
% Compute size and allocate data, added by Eric Johnson.  Need to apply to
% other variables.

numCols = n30;
numRows = 1 + n:18000 + n;
numRows = max(numRows) - n;

u030 = zeros(numRows,numCols);
v030 = zeros(numRows,numCols);
w030 = zeros(numRows,numCols);
t030 = zeros(numRows,numCols);

for k = 1:n30     

for i = 1 + n : 18000 + n  %18000 data points for 30 min
    
    u030(i - n,k) = u01(i);
    v030(i - n,k) = v01(i);
    w030(i - n,k) = w01(i);
    t030(i - n,k) = t01(i);
    
end

n = i;

end


%Cutting the data into 60 min sets
n = 0;
% Compute size and allocate data, added by Eric Johnson.  Need to apply to
% other variables.

numCols = n60;
numRows = 1 + n:36000 + n;
numRows = max(numRows) - n;

u060 = zeros(numRows,numCols);
v060 = zeros(numRows,numCols);
w060 = zeros(numRows,numCols);
t060 = zeros(numRows,numCols);


for k = 1:n60    

for i = 1 + n : 36000 + n  %36000 data points for 60 min
    
    u060(i - n,k) = u01(i);
    v060(i - n,k) = v01(i);
    w060(i - n,k) = w01(i);
    t060(i - n,k) = t01(i);
    
end

n = i;

end

%Cutting the data into 120 min sets
n = 0;
% Compute size and allocate data, added by Eric Johnson.  Need to apply to
% other variables.

numCols = n120;
numRows = 1 + n:72000 + n;
numRows = max(numRows) - n;

u0120 = zeros(numRows,numCols);
v0120 = zeros(numRows,numCols);
w0120 = zeros(numRows,numCols);
t0120 = zeros(numRows,numCols);

for k = 1:n120     

for i = 1 + n : 72000 + n  %72000 data points for 120 min
    
    u0120(i - n,k) = u01(i);
    v0120(i - n,k) = v01(i);
    w0120(i - n,k) = w01(i);
    t0120(i - n,k) = t01(i);
    
end

n = i;

end

%horizontal rotation on all datasets and other calculations

%5 min sets
if n5 >= 1

% Compute size and allocate data, added by Eric Johnson.  Need to apply to
% other variables.

numCols = n5;
numRows = 3000;

u05bar = zeros(numCols);
v05bar = zeros(numCols);
w05bar = zeros(numCols);
thetabar5 = zeros(numCols);

V5 = zeros(numRows,numCols);
thetai5 = zeros(numRows,numCols);
thetaprimei5 = zeros(numRows,numCols);
ui5 = zeros(numRows,numCols);
vi5 = zeros(numRows,numCols);
wi5 = zeros(numRows,numCols);
ti5 = zeros(numRows,numCols);
witi5 = zeros(numRows,numCols);
ui5bar = zeros(numCols);
uiwi5 = zeros(numRows,numCols);
viwi5 = zeros(numRows,numCols);

skewui5 = zeros(numCols);              
skewvi5 = zeros(numCols);             
skewwi5 = zeros(numCols);
skewti5 = zeros(numCols);

kurtui5 = zeros(numCols);               
kurtvi5 = zeros(numCols);              
kurtwi5 = zeros(numCols);
kurtti5 = zeros(numCols);

haarui5 = zeros(numCols);            
haarvi5 = zeros(numCols);
haarwi5 = zeros(numCols);
haarti5 = zeros(numCols);

witi5bar = zeros(numCols); 
ustar5 = zeros(numCols);
L5 = zeros(numCols);
FRICT5 = zeros(numCols);
SIGMAU5 = zeros(numCols);
SIGMAV5 = zeros(numCols);
SIGMAW5 = zeros(numCols);
SIGMAT5 = zeros(numCols);

fistar = zeros(numRows);
fir = zeros(numRows);


for k = 1:n5      
    
u05bar(k) = mean(u05(:,k));
v05bar(k)= mean(v05(:,k));
w05bar(k) = mean(w05(:,k)); %added recently
thetabar5(k) = atand(u05bar(k)/v05bar(k));


for i = 1:3000

V5(i,k) = sqrt((u05(i,k)^2) + (v05(i,k)^2));
thetai5(i,k) = atand(u05(i,k)/((v05(i,k)+0.001)));   %0.001 is added to remove ZERO values in v(:)
thetaprimei5(i,k) = thetai5(i,k) - thetabar5(k) ;
ui5(i,k) = V5(i,k) * cosd(thetaprimei5(i,k));
vi5(i,k) = V5(i,k) * sind(thetaprimei5(i,k));
wi5(i,k) = w05(i,k) - w05bar(k); 
ti5(i,k) = t05(i,k) - mean(t05(:,k));
witi5(i,k) = wi5(i,k) * ti5(i,k); 



end

ui5bar(k) = mean(ui5(:,k));
    
for i = 1:3000
    
uiwi5(i,k) = (ui5(i,k) - ui5bar(k)) * wi5(i,k); %to calculate friction velocity
viwi5(i,k) = (vi5(i,k) - mean(vi5(:,k))) * wi5(i,k); %to calculate friction velocity %ADDED

end



skewui5(k) = skewness(ui5(:,k));              %skewness needs to be in the range of -2 to 2 (hard flagged) or -1 to 1 (soft flagged)
skewvi5(k) = skewness(vi5(:,k));              %ideal skewness = 0
skewwi5(k) = skewness(wi5(:,k));
skewti5(k) = skewness(ti5(:,k));

kurtui5(k) = kurtosis(ui5(:,k));              %kurtosis needs to be in the range of 1 to 8 (hard flagged) or 2 to 5 (soft flagged) 
kurtvi5(k) = kurtosis(vi5(:,k));              %ideal kurtosis = 3
kurtwi5(k) = kurtosis(wi5(:,k));
kurtti5(k) = kurtosis(ti5(:,k));

%Haar transform calculation for 2 half windows of one set
%u
variance1ui5 = var(ui5(1:1500,k));
variance2ui5 = var(ui5(1501:3000,k));
variancetotalui5 = var(ui5(:,k));
%v
variance1vi5 = var(vi5(1:1500,k));
variance2vi5 = var(vi5(1501:3000,k));
variancetotalvi5 = var(vi5(:,k));
%w
variance1wi5 = var(wi5(1:1500,k));
variance2wi5 = var(wi5(1501:3000,k));
variancetotalwi5 = var(wi5(:,k));
%t
variance1ti5 = var(ti5(1:1500,k));
variance2ti5 = var(ti5(1501:3000,k));
variancetotalti5 = var(ti5(:,k));



haarui5(k)= abs((variance1ui5 - variance2ui5)/variancetotalui5);            %Haar transform must be less than 2 (soft) and 3 (hard ideal value is 0
haarvi5(k) = abs((variance1vi5 - variance2vi5)/variancetotalvi5);
haarwi5(k) = abs((variance1wi5 - variance2wi5)/variancetotalwi5);
haarti5(k) = abs((variance1ti5 - variance2ti5)/variancetotalti5);
%To calculate the average of w'T' for each k
witi5bar(k) = mean(witi5(:,k)); 
%To calculate friction velocity
ustar5(k) = power((mean(uiwi5(:,k)))^2 + (mean(viwi5(:,k)))^2,0.25);
%To calculate the Monin-Obukhov length, L
L5(k) = (-(power(ustar5(k),3)) * mean(t05(:,k)))/(9.81*0.4*witi5bar(k));
%To calculate friction temperature
FRICT5(k) = witi5bar(k)/ustar5(k);
%To calculate the turbulent characteristics
SIGMAU5(k) = std(ui5(:,k));
SIGMAV5(k) = std(vi5(:,k));
SIGMAW5(k) = std(wi5(:,k));
SIGMAT5(k) = std(ti5(:,k));

end



%Final w'T' for record
witi5bar1 = mean(witi5bar(:));  %Newly added


%Using linear regression function

temp = (0:k-1);
x = temp(:);
y = witi5bar(:,1) - witi5bar1;

[a0 a1] = linear_regression(x,y); %Linear correlation coefficients are obtained

end

%Calculation of RFE (random flux error) using L = 5 min
for i = 1:k
    
fistar(i) = witi5bar(i) - witi5bar1 - a0 - a1*(i-1);      %a0 must be substracted from record mean witi5bar1

end

stdfistar = std(fistar(:));

N = k;

RFE = stdfistar / (abs(witi5bar1)*N^0.5); %Nonstationarity of turbulent fluxes test for record. Must be less than 0.25

for i = 1:k
    
fir(i) = fistar(i) + witi5bar(i) - witi5bar1;

end

stdfir = std(fir(:));

RN = stdfir / (abs(witi5bar1)*N^0.5);   %Nonstationarity of turbulent fluxes due to mesoscale motion. Must be less than 0.25


%10 min sets

if n10 >= 1

% Compute size and allocate data, added by Eric Johnson.  Need to apply to
% other variables.

numCols = n10;
numRows = 6000;

u010bar = zeros(numCols);
v010bar = zeros(numCols);
w010bar = zeros(numCols);
thetabar10 = zeros(numCols);

V10 = zeros(numRows,numCols);
thetai10 = zeros(numRows,numCols);
thetaprimei10 = zeros(numRows,numCols);
ui10 = zeros(numRows,numCols);
vi10 = zeros(numRows,numCols);
wi10 = zeros(numRows,numCols);
ti10 = zeros(numRows,numCols);
witi10 = zeros(numRows,numCols);
ui10bar = zeros(numCols);
uiwi10 = zeros(numRows,numCols);
viwi10 = zeros(numRows,numCols);

skewui10 = zeros(numCols);              
skewvi10 = zeros(numCols);             
skewwi10 = zeros(numCols);
skewti10 = zeros(numCols);

kurtui10 = zeros(numCols);               
kurtvi10 = zeros(numCols);              
kurtwi10 = zeros(numCols);
kurtti10 = zeros(numCols);

haarui10 = zeros(numCols);            
haarvi10 = zeros(numCols);
haarwi10 = zeros(numCols);
haarti10 = zeros(numCols);

witi10bar = zeros(numCols); 
ustar10 = zeros(numCols);
L10 = zeros(numCols);
FRICT10 = zeros(numCols);
SIGMAU10 = zeros(numCols);
SIGMAV10 = zeros(numCols);
SIGMAW10 = zeros(numCols);
SIGMAT10 = zeros(numCols);

 
for k = 1:n10      
    
u010bar(k) = mean(u010(:,k));
v010bar(k)= mean(v010(:,k));
w010bar(k) = mean(w010(:,k)); %added recently
thetabar10(k) = atand(u010bar(k)/v010bar(k));


for i = 1:6000

V10(i,k) = sqrt((u010(i,k)^2) + (v010(i,k)^2));
thetai10(i,k) = atand(u010(i,k)/((v010(i,k)+0.001)));   %0.001 is added to remove ZERO values in v(:)
thetaprimei10(i,k) = thetai10(i,k) - thetabar10(k) ;
ui10(i,k) = V10(i,k) * cosd(thetaprimei10(i,k));
vi10(i,k) = V10(i,k) * sind(thetaprimei10(i,k));
wi10(i,k) = w010(i,k) - w010bar(k); %added recently
ti10(i,k) = t010(i,k) - mean(t010(:,k));
witi10(i,k) = wi10(i,k) * ti10(i,k); %newly added



end

ui10bar(k) = mean(ui10(:,k));
    
for i = 1:6000
    
uiwi10(i,k) = (ui10(i,k) - ui10bar(k)) * wi10(i,k); %to calculate friction velocity
viwi10(i,k) = (vi10(i,k) - mean(vi10(:,k))) * wi10(i,k); %to calculate friction velocity %ADDED
end

skewui10(k) = skewness(ui10(:,k));              %skewness needs to be in the range of -2 to 2 (hard flagged) or -1 to 1 (soft flagged)
skewvi10(k) = skewness(vi10(:,k));              %ideal skewness = 0
skewwi10(k) = skewness(wi10(:,k));
skewti10(k) = skewness(ti10(:,k));

kurtui10(k) = kurtosis(ui10(:,k));              %kurtosis needs to be in the range of 1 to 8 (hard flagged) or 2 to 5 (soft flagged) 
kurtvi10(k) = kurtosis(vi10(:,k));              %ideal kurtosis = 3
kurtwi10(k) = kurtosis(wi10(:,k));
kurtti10(k) = kurtosis(ti10(:,k));

%Haar transform calculation for 2 half windows of one set
%u
variance1ui10 = var(ui10(1:3000,k));
variance2ui10 = var(ui10(3001:6000,k));
variancetotalui10 = var(ui10(:,k));
%v
variance1vi10 = var(vi10(1:3000,k));
variance2vi10 = var(vi10(3001:6000,k));
variancetotalvi10 = var(vi10(:,k));
%w
variance1wi10 = var(wi10(1:3000,k));
variance2wi10 = var(wi10(3001:6000,k));
variancetotalwi10 = var(wi10(:,k));
%t
variance1ti10 = var(ti10(1:3000,k));
variance2ti10 = var(ti10(3001:6000,k));
variancetotalti10 = var(ti10(:,k));



haarui10(k)= abs((variance1ui10 - variance2ui10)/variancetotalui10);            %Haar transform must be less than 2 (soft) and 3 (hard ideal value is 0
haarvi10(k) = abs((variance1vi10 - variance2vi10)/variancetotalvi10);
haarwi10(k) = abs((variance1wi10 - variance2wi10)/variancetotalwi10);
haarti10(k) = abs((variance1ti10 - variance2ti10)/variancetotalti10);
%To calculate the average of w'T' for each k
witi10bar(k) = mean(witi10(:,k)); %Newly added
%To calculate friction velocity
ustar10(k) = power((mean(uiwi10(:,k)))^2 + (mean(viwi10(:,k)))^2,0.25);
%To calculate the Monin-Obukhov length, L
L10(k) = (-(power(ustar10(k),3)) * mean(t010(:,k)))/(9.81*0.4*witi10bar(k));
%To calculate friction temperature
FRICT10(k) = witi10bar(k)/ustar10(k);
%To calculate the turbulent characteristics
SIGMAU10(k) = std(ui10(:,k));
SIGMAV10(k) = std(vi10(:,k));
SIGMAW10(k) = std(wi10(:,k));
SIGMAT10(k) = std(ti10(:,k));
end
%Final w'T' for record
witi10bar1 = mean(witi10bar(:));  %Newly added

%15 min sets

% Compute size and allocate data, added by Eric Johnson.  Need to apply to
% other variables.

numCols = n15;
numRows = 9000;

u015bar = zeros(numCols);
v015bar = zeros(numCols);
w015bar = zeros(numCols);
thetabar15 = zeros(numCols);

V15 = zeros(numRows,numCols);
thetai15 = zeros(numRows,numCols);
thetaprimei15 = zeros(numRows,numCols);
ui15 = zeros(numRows,numCols);
vi15 = zeros(numRows,numCols);
wi15 = zeros(numRows,numCols);
ti15 = zeros(numRows,numCols);
witi15 = zeros(numRows,numCols);
ui15bar = zeros(numCols);
uiwi15 = zeros(numRows,numCols);
viwi15 = zeros(numRows,numCols);

skewui15 = zeros(numCols);              
skewvi15 = zeros(numCols);             
skewwi15 = zeros(numCols);
skewti15 = zeros(numCols);

kurtui15 = zeros(numCols);               
kurtvi15 = zeros(numCols);              
kurtwi15 = zeros(numCols);
kurtti15 = zeros(numCols);

haarui15 = zeros(numCols);            
haarvi15 = zeros(numCols);
haarwi15 = zeros(numCols);
haarti15 = zeros(numCols);

witi15bar = zeros(numCols); 
ustar15 = zeros(numCols);
L15 = zeros(numCols);
FRICT15 = zeros(numCols);
SIGMAU15 = zeros(numCols);
SIGMAV15 = zeros(numCols);
SIGMAW15 = zeros(numCols);
SIGMAT15 = zeros(numCols);

for k = 1:n15      
    
u015bar(k) = mean(u015(:,k));
v015bar(k)= mean(v015(:,k));
w015bar(k) = mean(w015(:,k)); %added recently
thetabar15(k) = atand(u015bar(k)/v015bar(k));


for i = 1:9000

V15(i,k) = sqrt((u015(i,k)^2) + (v015(i,k)^2));
thetai15(i,k) = atand(u015(i,k)/((v015(i,k)+0.001)));   %0.001 is added to remove ZERO values in v(:)
thetaprimei15(i,k) = thetai15(i,k) - thetabar15(k) ;
ui15(i,k) = V15(i,k) * cosd(thetaprimei15(i,k));
vi15(i,k) = V15(i,k) * sind(thetaprimei15(i,k));
wi15(i,k) = w015(i,k) - w015bar(k); %added recently
ti15(i,k) = t015(i,k) - mean(t015(:,k));
witi15(i,k) = wi15(i,k) * ti15(i,k); %newly added


end

ui15bar(k) = mean(ui15(:,k));
    
for i = 1:9000
    
uiwi15(i,k) = (ui15(i,k) - ui15bar(k)) * wi15(i,k); %to calculate friction velocity
viwi15(i,k) = (vi15(i,k) - mean(vi15(:,k))) * wi15(i,k); %to calculate friction velocity %ADDED

end

skewui15(k) = skewness(ui15(:,k));              %skewness needs to be in the range of -2 to 2 (hard flagged) or -1 to 1 (soft flagged)
skewvi15(k) = skewness(vi15(:,k));              %ideal skewness = 0
skewwi15(k) = skewness(wi15(:,k));
skewti15(k) = skewness(ti15(:,k));

kurtui15(k) = kurtosis(ui15(:,k));              %kurtosis needs to be in the range of 1 to 8 (hard flagged) or 2 to 5 (soft flagged) 
kurtvi15(k) = kurtosis(vi15(:,k));              %ideal kurtosis = 3
kurtwi15(k) = kurtosis(wi15(:,k));
kurtti15(k) = kurtosis(ti15(:,k));

%Haar transform calculation for 2 half windows of one set
%u
variance1ui15 = var(ui15(1:4500,k));
variance2ui15 = var(ui15(4501:9000,k));
variancetotalui15 = var(ui15(:,k));
%v
variance1vi15 = var(vi15(1:4500,k));
variance2vi15 = var(vi15(4501:9000,k));
variancetotalvi15 = var(vi15(:,k));
%w
variance1wi15 = var(wi15(1:4500,k));
variance2wi15 = var(wi15(4501:9000,k));
variancetotalwi15 = var(wi15(:,k));
%t
variance1ti15 = var(ti15(1:4500,k));
variance2ti15 = var(ti15(4501:9000,k));
variancetotalti15 = var(ti15(:,k));



haarui15(k)= abs((variance1ui15 - variance2ui15)/variancetotalui15);            %Haar transform must be less than 2 (soft) and 3 (hard ideal value is 0
haarvi15(k) = abs((variance1vi15 - variance2vi15)/variancetotalvi15);
haarwi15(k) = abs((variance1wi15 - variance2wi15)/variancetotalwi15);
haarti15(k) = abs((variance1ti15 - variance2ti15)/variancetotalti15);
%To calculate the average of w'T' for each k
witi15bar(k) = mean(witi15(:,k)); %Newly added
%To calculate friction velocity
ustar15(k) = power((mean(uiwi15(:,k)))^2 + (mean(viwi15(:,k)))^2,0.25);
%To calculate the Monin-Obukhov length, L
L15(k) = (-(power(ustar15(k),3)) * mean(t015(:,k)))/(9.81*0.4*witi15bar(k));
%To calculate friction temperature
FRICT15(k) = witi15bar(k)/ustar15(k);
%To calculate the turbulent characteristics
SIGMAU15(k) = std(ui15(:,k));
SIGMAV15(k) = std(vi15(:,k));
SIGMAW15(k) = std(wi15(:,k));
SIGMAT15(k) = std(ti15(:,k));
end
%Final w'T' for record
witi15bar1 = mean(witi15bar(:));  %Newly added

end

%30 min sets

if n30 >= 1
 % Compute size and allocate data, added by Eric Johnson.  Need to apply to
% other variables.

numCols = n30;
numRows = 18000;

u030bar = zeros(numCols);
v030bar = zeros(numCols);
w030bar = zeros(numCols);
thetabar30 = zeros(numCols);

V30 = zeros(numRows,numCols);
thetai30 = zeros(numRows,numCols);
thetaprimei30 = zeros(numRows,numCols);
ui30 = zeros(numRows,numCols);
vi30 = zeros(numRows,numCols);
wi30 = zeros(numRows,numCols);
ti30 = zeros(numRows,numCols);
witi30 = zeros(numRows,numCols);
ui30bar = zeros(numCols);
vi30bar = zeros(numCols);
uiwi30 = zeros(numRows,numCols);
viwi30 = zeros(numRows,numCols);

skewui30 = zeros(numCols);              
skewvi30 = zeros(numCols);             
skewwi30 = zeros(numCols);
skewti30 = zeros(numCols);

kurtui30 = zeros(numCols);               
kurtvi30 = zeros(numCols);              
kurtwi30 = zeros(numCols);
kurtti30 = zeros(numCols);

haarui30 = zeros(numCols);            
haarvi30 = zeros(numCols);
haarwi30 = zeros(numCols);
haarti30 = zeros(numCols);

witi30bar = zeros(numCols); 
ustar30 = zeros(numCols);
L30 = zeros(numCols);
FRICT30 = zeros(numCols);
SIGMAU30 = zeros(numCols);
SIGMAV30 = zeros(numCols);
SIGMAW30 = zeros(numCols);
SIGMAT30 = zeros(numCols);
uii=zeros(numRows,numCols);
vii=zeros(numRows,numCols);

for k = 1:n30      
    
u030bar(k) = mean(u030(:,k));
v030bar(k)= mean(v030(:,k));
w030bar(k) = mean(w030(:,k)); %added recently
thetabar30(k) = atand(u030bar(k)/v030bar(k));
%V30bar(k) = sqrt((u030bar(k))^2 + (v030bar(k)^2)); %!!!

for i = 1:18000

V30(i,k) = sqrt((u030(i,k)^2) + (v030(i,k)^2)); %!!!
thetai30(i,k) = atand(u030(i,k)/((v030(i,k)+0.001))); %0.001 is added to remove ZERO values in v(:)
thetaprimei30(i,k) = thetai30(i,k) - thetabar30(k) ;
ui30(i,k) = V30(i,k) * cosd(thetaprimei30(i,k)); %!!!
vi30(i,k) = V30(i,k) * sind(thetaprimei30(i,k));%!!!
wi30(i,k) = w030(i,k) - w030bar(k); %added recently
ti30(i,k) = t030(i,k) - mean(t030(:,k));
witi30(i,k) = wi30(i,k) * ti30(i,k); %newly added



end


ui30bar(k) = mean(ui30(:,k)); %!!!
vi30bar(k) = mean(ui30(:,k));%!!!

for i = 1:18000

    uii(i,k) = ui30(i,k) - ui30bar(k); %!!!
    vii(i,k) = vi30(i,k) - vi30bar(k); %!!!
    uiwi30(i,k) = uii(i,k) * wi30(i,k); %to calculate friction velocity
    viwi30(i,k) = vii(i,k) * wi30(i,k); %to calculate friction velocity %ADDED

end

skewui30(k) = skewness(ui30(:,k));              %skewness needs to be in the range of -2 to 2 (hard flagged) or -1 to 1 (soft flagged)
skewvi30(k) = skewness(vi30(:,k));              %ideal skewness = 0
skewwi30(k) = skewness(wi30(:,k));
skewti30(k) = skewness(ti30(:,k));

kurtui30(k) = kurtosis(ui30(:,k));              %kurtosis needs to be in the range of 1 to 8 (hard flagged) or 2 to 5 (soft flagged) 
kurtvi30(k) = kurtosis(vi30(:,k));              %ideal kurtosis = 3
kurtwi30(k) = kurtosis(wi30(:,k));
kurtti30(k) = kurtosis(ti30(:,k));


%Haar transform calculation for 2 half windows of one set
%u
variance1ui30 = var(ui30(1:9000,k));
variance2ui30 = var(ui30(9001:18000,k));
variancetotalui30 = var(ui30(:,k));
%v
variance1vi30 = var(vi30(1:9000,k));
variance2vi30 = var(vi30(9001:18000,k));
variancetotalvi30 = var(vi30(:,k));
%w
variance1wi30 = var(wi30(1:9000,k));
variance2wi30 = var(wi30(9001:18000,k));
variancetotalwi30 = var(wi30(:,k));
%t
variance1ti30 = var(ti30(1:9000,k));
variance2ti30 = var(ti30(9001:18000,k));
variancetotalti30 = var(ti30(:,k));



haarui30(k)= abs((variance1ui30 - variance2ui30)/variancetotalui30);            %Haar transform must be less than 2 (soft) and 3 (hard ideal value is 0
haarvi30(k) = abs((variance1vi30 - variance2vi30)/variancetotalvi30);
haarwi30(k) = abs((variance1wi30 - variance2wi30)/variancetotalwi30);
haarti30(k) = abs((variance1ti30 - variance2ti30)/variancetotalti30);
%To calculate the average of w'T' for each k
witi30bar(k) = mean(witi30(:,k)); %Newly added
%To calculate friction velocity
ustar30(k) = power((mean(uiwi30(:,k)))^2 + (mean(viwi30(:,k)))^2,0.25);
%To calculate the Monin-Obukhov length, L
L30(k) = (-(power(ustar30(k),3)) * mean(t030(:,k)))/(9.81*0.4*witi30bar(k));
%To calculate friction temperature
FRICT30(k) = witi30bar(k)/ustar30(k);
%To calculate the turbulent characteristics
SIGMAU30(k) = std(ui30(:,k));
SIGMAV30(k) = std(vi30(:,k));
SIGMAW30(k) = std(wi30(:,k));
SIGMAT30(k) = std(ti30(:,k));
end
%Final w'T' for record
witi30bar1 = mean(witi30bar(:));  %Newly added

end 

%60 min sets

if n60 >= 1
% Compute size and allocate data, added by Eric Johnson.  Need to apply to
% other variables.

numCols = n60;
numRows = 36000;

u060bar = zeros(numCols);
v060bar = zeros(numCols);
w060bar = zeros(numCols);
thetabar60 = zeros(numCols);

V60 = zeros(numRows,numCols);
thetai60 = zeros(numRows,numCols);
thetaprimei60 = zeros(numRows,numCols);
ui60 = zeros(numRows,numCols);
vi60 = zeros(numRows,numCols);
wi60 = zeros(numRows,numCols);
ti60 = zeros(numRows,numCols);
witi60 = zeros(numRows,numCols);
ui60bar = zeros(numCols);
uiwi60 = zeros(numRows,numCols);
viwi60 = zeros(numRows,numCols);

skewui60 = zeros(numCols);              
skewvi60 = zeros(numCols);             
skewwi60 = zeros(numCols);
skewti60 = zeros(numCols);

kurtui60 = zeros(numCols);               
kurtvi60 = zeros(numCols);              
kurtwi60 = zeros(numCols);
kurtti60 = zeros(numCols);

haarui60 = zeros(numCols);            
haarvi60 = zeros(numCols);
haarwi60 = zeros(numCols);
haarti60 = zeros(numCols);

witi60bar = zeros(numCols); 
ustar60 = zeros(numCols);
L60 = zeros(numCols);
FRICT60 = zeros(numCols);
SIGMAU60 = zeros(numCols);
SIGMAV60 = zeros(numCols);
SIGMAW60 = zeros(numCols);
SIGMAT60 = zeros(numCols);
  
for k = 1:n60      
    
u060bar(k) = mean(u060(:,k));
v060bar(k)= mean(v060(:,k));
w060bar(k) = mean(w060(:,k)); %added recently
thetabar60(k) = atand(u060bar(k)/v060bar(k));


for i = 1:36000

V60(i,k) = sqrt((u060(i,k)^2) + (v060(i,k)^2));
thetai60(i,k) = atand(u060(i,k)/((v060(i,k)+0.001)));       %0.001 is added to remove ZERO values in v(:)
thetaprimei60(i,k) = thetai60(i,k) - thetabar60(k) ;
ui60(i,k) = V60(i,k) * cosd(thetaprimei60(i,k));
vi60(i,k) = V60(i,k) * sind(thetaprimei60(i,k));
wi60(i,k) = w060(i,k) - w060bar(k); %added recently
ti60(i,k) = t060(i,k) - mean(t060(:,k));
witi60(i,k) = wi60(i,k) * ti60(i,k); %newly added

viwi60(i,k) = vi60(i,k) * wi60(i,k);    %to calculate friction velocity

end

ui60bar(k) = mean(ui60(:,k));
    
for i = 1:36000
    
uiwi60(i,k) = (ui60(i,k) - ui60bar(k)) * wi60(i,k); %to calculate friction velocity

end

skewui60(k) = skewness(ui60(:,k));              %skewness needs to be in the range of -2 to 2 (hard flagged) or -1 to 1 (soft flagged)
skewvi60(k) = skewness(vi60(:,k));              %ideal skewness = 0
skewwi60(k) = skewness(wi60(:,k));
skewti60(k) = skewness(ti60(:,k));

kurtui60(k) = kurtosis(ui60(:,k));              %kurtosis needs to be in the range of 1 to 8 (hard flagged) or 2 to 5 (soft flagged) 
kurtvi60(k) = kurtosis(vi60(:,k));              %ideal kurtosis = 3
kurtwi60(k) = kurtosis(wi60(:,k));
kurtti60(k) = kurtosis(ti60(:,k));


%Haar transform calculation for 2 half windows of one set
%u
variance1ui60 = var(ui60(1:18000,k));
variance2ui60 = var(ui60(18001:36000,k));
variancetotalui60 = var(ui60(:,k));
%v
variance1vi60 = var(vi60(1:18000,k));
variance2vi60 = var(vi60(18001:36000,k));
variancetotalvi60 = var(vi60(:,k));
%w
variance1wi60 = var(wi60(1:18000,k));
variance2wi60 = var(wi60(18001:36000,k));
variancetotalwi60 = var(wi60(:,k));
%t
variance1ti60 = var(ti60(1:18000,k));
variance2ti60 = var(ti60(18001:36000,k));
variancetotalti60 = var(ti60(:,k));



haarui60(k)= abs((variance1ui60 - variance2ui60)/variancetotalui60);            %Haar transform must be less than 2 (soft) and 3 (hard ideal value is 0
haarvi60(k) = abs((variance1vi60 - variance2vi60)/variancetotalvi60);
haarwi60(k) = abs((variance1wi60 - variance2wi60)/variancetotalwi60);
haarti60(k) = abs((variance1ti60 - variance2ti60)/variancetotalti60);
%To calculate the average of w'T' for each k
witi60bar(k) = mean(witi60(:,k)); %Newly added
%To calculate friction velocity
ustar60(k) = power((mean(uiwi60(:,k)))^2 + (mean(viwi60(:,k)))^2,0.25);
%To calculate the Monin-Obukhov length, L
L60(k) = (-(power(ustar60(k),3)) * mean(t060(:,k)))/(9.81*0.4*witi60bar(k));
%To calculate friction temperature
FRICT60(k) = witi60bar(k)/ustar60(k);
%To calculate the turbulent characteristics
SIGMAU60(k) = std(ui60(:,k));
SIGMAV60(k) = std(vi60(:,k));
SIGMAW60(k) = std(wi60(:,k));
SIGMAT60(k) = std(ti60(:,k));
end

%Final w'T' for record
witi60bar1 = mean(witi60bar(:));  %Newly added

end

%120 min sets

if n120 >= 1
% Compute size and allocate data, added by Eric Johnson.  Need to apply to
% other variables.

numCols = n120;
numRows = 72000;

u0120bar = zeros(numCols);
v0120bar = zeros(numCols);
w0120bar = zeros(numCols);
thetabar120 = zeros(numCols);

V120 = zeros(numRows,numCols);
thetai120 = zeros(numRows,numCols);
thetaprimei120 = zeros(numRows,numCols);
ui120 = zeros(numRows,numCols);
vi120 = zeros(numRows,numCols);
wi120 = zeros(numRows,numCols);
ti120 = zeros(numRows,numCols);
witi120 = zeros(numRows,numCols);
ui120bar = zeros(numCols);
uiwi120 = zeros(numRows,numCols);
viwi120 = zeros(numRows,numCols);

skewui120 = zeros(numCols);              
skewvi120 = zeros(numCols);             
skewwi120 = zeros(numCols);
skewti120 = zeros(numCols);

kurtui120 = zeros(numCols);               
kurtvi120 = zeros(numCols);              
kurtwi120 = zeros(numCols);
kurtti120 = zeros(numCols);

haarui120 = zeros(numCols);            
haarvi120 = zeros(numCols);
haarwi120 = zeros(numCols);
haarti120 = zeros(numCols);

witi120bar = zeros(numCols); 
ustar120 = zeros(numCols);
L120 = zeros(numCols);
FRICT120 = zeros(numCols);
SIGMAU120 = zeros(numCols);
SIGMAV120 = zeros(numCols);
SIGMAW120 = zeros(numCols);
SIGMAT120 = zeros(numCols);
 
for k = 1:n120      
    
u0120bar(k) = mean(u0120(:,k));
v0120bar(k)= mean(v0120(:,k));
w0120bar(k) = mean(w0120(:,k)); %added recently
thetabar120(k) = atand(u0120bar(k)/v0120bar(k));

for i = 1:72000

V120(i,k) = sqrt((u0120(i,k)^2) + (v0120(i,k)^2));
thetai120(i,k) = atand(u0120(i,k)/((v0120(i,k)+0.001)));    %0.001 is added to remove ZERO values in v(:)
thetaprimei120(i,k) = thetai120(i,k) - thetabar120(k) ;
ui120(i,k) = V120(i,k) * cosd(thetaprimei120(i,k));
vi120(i,k) = V120(i,k) * sind(thetaprimei120(i,k));
wi120(i,k) = w0120(i,k) - w0120bar(k); %added recently
ti120(i,k) = t0120(i,k) - mean(t0120(:,k));
witi120(i,k) = wi120(i,k) * ti120(i,k); %newly added



end

ui120bar(k) = mean(ui120(:,k));
    
for i = 1:72000
    
uiwi120(i,k) = (ui120(i,k) - ui120bar(k)) * wi120(i,k); %to calculate friction velocity
viwi120(i,k) = (vi120(i,k) - mean(vi120(:,k))) * wi120(i,k); %to calculate friction velocity %ADDED

end

skewui120(k) = skewness(ui120(:,k));              %skewness needs to be in the range of -2 to 2 (hard flagged) or -1 to 1 (soft flagged)
skewvi120(k) = skewness(vi120(:,k));              %ideal skewness = 0
skewwi120(k) = skewness(wi120(:,k));
skewti120(k) = skewness(ti120(:,k));

kurtui120(k) = kurtosis(ui120(:,k));              %kurtosis needs to be in the range of 1 to 8 (hard flagged) or 2 to 5 (soft flagged) 
kurtvi120(k) = kurtosis(vi120(:,k));              %ideal kurtosis = 3
kurtwi120(k) = kurtosis(wi120(:,k));
kurtti120(k) = kurtosis(ti120(:,k));


%Haar transform calculation for 2 half windows of one set
%u
variance1ui120 = var(ui120(1:36000,k));
variance2ui120 = var(ui120(36001:72000,k));
variancetotalui120 = var(ui120(:,k));
%v
variance1vi120 = var(vi120(1:36000,k));
variance2vi120 = var(vi120(36001:72000,k));
variancetotalvi120 = var(vi120(:,k));
%w
variance1wi120 = var(wi120(1:36000,k));
variance2wi120 = var(wi120(36001:72000,k));
variancetotalwi120 = var(wi120(:,k));
%t
variance1ti120 = var(ti120(1:36000,k));
variance2ti120 = var(ti120(36001:72000,k));
variancetotalti120 = var(ti120(:,k));



haarui120(k)= abs((variance1ui120 - variance2ui120)/variancetotalui120);            %Haar transform must be less than 2 (soft) and 3 (hard ideal value is 0
haarvi120(k) = abs((variance1vi120 - variance2vi120)/variancetotalvi120);
haarwi120(k) = abs((variance1wi120 - variance2wi120)/variancetotalwi120);
haarti120(k) = abs((variance1ti120 - variance2ti120)/variancetotalti120);
%To calculate the average of w'T' for each k
witi120bar(k) = mean(witi120(:,k)); %Newly added
%To calculate friction velocity
ustar120(k) = power((mean(uiwi120(:,k)))^2 + (mean(viwi120(:,k)))^2,0.25);
%To calculate the Monin-Obukhov length, L
L120(k) = (-(power(ustar120(k),3)) * mean(t0120(:,k)))/(9.81*0.4*witi120bar(k));
%To calculate friction temperature
FRICT120(k) = witi120bar(k)/ustar120(k);
%To calculate the turbulent characteristics
SIGMAU120(k) = std(ui120(:,k));
SIGMAV120(k) = std(vi120(:,k));
SIGMAW120(k) = std(wi120(:,k));
SIGMAT120(k) = std(ti120(:,k));
end
%Final w'T' for record
witi120bar1 = mean(witi120bar(:));  %Newly added

end

%Writing the computed values into a text file

%Creating a text file for 5 min

file_id = fopen('ProcessedData_5_min.txt', 'w');

%Writing code must be written in between creating and closing of text file

fprintf(file_id, 'Date:\n\n'); 
fprintf(file_id, 'Initial time of measurement:\n\n');
fprintf(file_id, 'Height of measurement:\n\n');
fprintf(file_id, 'Number of data in the entire original record: %g\n\n', length(u0));
fprintf(file_id, 'Number of spikes in the entire record: %g (%g percent)\n\n\n\n',length(u0) - length(u01),(((length(u0) - length(u01))/length(u0))*100));

fprintf(file_id, 'L = 5 min \n\n');
fprintf(file_id, 'Nonstationarity tests of turbulent fluxes\n (only available for L = 5 min)\n\n');
fprintf(file_id, 'RFE = %g   (Note: Categorized as soft flagged if > 0.25)\n', RFE);
fprintf(file_id, 'RN = %g   (Note: Categorized as soft flagged if > 0.25)\n\n\n', RN);

%kurtosis
fprintf(file_id, 'KURTOSIS\n\n');
fprintf(file_id, 'no\t kurtU\t kurtV\t kurtW\t kurtT\n\n');

for b = 1:n5
    fprintf(file_id, '%g\t%g\t%g\t%g\t%g\n\n', b, kurtui5(b), kurtvi5(b), kurtwi5(b), kurtti5(b));
end

fprintf(file_id, 'NOTE: Kurtosis needs to be in the range of 1 to 8 (hard flagged) or 2 to 5 (soft flagged). Ideal value of kurtosis is 3\n\n\n');

%skewness
fprintf(file_id, 'SKEWNESS\n\n');
fprintf(file_id, 'no\t skewU\t skewV\t skewW\t skewT\n\n');

for b = 1:n5
    fprintf(file_id, '%g\t%g\t%g\t%g\t%g\n\n', b, skewui5(b), skewvi5(b), skewwi5(b), skewti5(b));
end

fprintf(file_id, 'NOTE: Skewness needs to be in the range of -2 to 2 (hard flagged) or -1 to 1 (soft flagged). Ideal value of skewness is 0\n\n\n');


%Haar transform
fprintf(file_id, 'HAAR TRANSFORM\n\n');
fprintf(file_id, 'no\t HaarU\t HaarV\t HaarW\t HaarT\n\n');

for b = 1:n5
    fprintf(file_id, '%g\t%g\t%g\t%g\t%g\n\n', b, haarui5(b), haarvi5(b), haarwi5(b), haarti5(b));
end

fprintf(file_id, 'NOTE: Haar transform must be less than 2 (soft) and 3 (hard), ideal value is 0\n\n\n');

%Turbulence statistics
fprintf(file_id, 'Turbulence statistics\n\n');
fprintf(file_id, 'no\t SigU\t\t SigV\t\t SigW\t\t SigT\t\t u*\n\n');

for b = 1:n5
    fprintf(file_id, '%g\t%6.3g\t%6.3g\t%6.3g\t%6.3g\t%6.3g\n\n', b, SIGMAU5(b), SIGMAV5(b), SIGMAW5(b), SIGMAT5(b), ustar5(b));
end

fprintf(file_id, 'CONTINUED...\n\n');

fprintf(file_id, 'no\t wT\t\t ui\t\t vi\t\t V\t\t L\t\t CD\n\n');

for b = 1:n5
    fprintf(file_id, '%g\t%6.3g\t%6.3g\t%6.3g\t%6.3g\t%6.3g  %6.3g\n\n', b, witi5bar(b), u05bar(b), v05bar(b), mean(V5(:,b)), L5(b),(ustar5(b)/mean(V5(:,b)))^2);
end

fprintf(file_id, 'CONTINUED...\n\n');

fprintf(file_id, 'no\t T*\t\t uw\n\n');

for b = 1:n5
    fprintf(file_id, '%g\t%6.3g\t%6.3g\n\n', b, FRICT5(b), mean(uiwi5(:,b)));
end

%Closing the created text file
fclose(file_id);

%Creating a text file for 10 min

file_id1 = fopen('ProcessedData_10_min.txt', 'w');

%Writing code must be written in between creating and closing of text file

fprintf(file_id1, 'Date:\n\n'); 
fprintf(file_id1, 'Initial time of measurement:\n\n');
fprintf(file_id1, 'Height of measurement:\n\n');
fprintf(file_id1, 'Number of data in the entire original record: %g\n\n', length(u0));
fprintf(file_id1, 'Number of spikes in the entire record: %g (%g percent)\n\n\n\n',length(u0) - length(u01),(((length(u0) - length(u01))/length(u0))*100));

fprintf(file_id1, 'L = 10 min \n\n');


%kurtosis
fprintf(file_id1, 'KURTOSIS\n\n');
fprintf(file_id1, 'no\t kurtU\t kurtV\t kurtW\t kurtT\n\n');

for b = 1:n10
    fprintf(file_id1, '%g\t%g\t%g\t%g\t%g\n\n', b, kurtui10(b), kurtvi10(b), kurtwi10(b), kurtti10(b));
end

fprintf(file_id1, 'NOTE: Kurtosis needs to be in the range of 1 to 8 (hard flagged) or 2 to 5 (soft flagged). Ideal value of kurtosis is 3\n\n\n');

%skewness
fprintf(file_id1, 'SKEWNESS\n\n');
fprintf(file_id1, 'no\t skewU\t skewV\t skewW\t skewT\n\n');

for b = 1:n10
    fprintf(file_id1, '%g\t%g\t%g\t%g\t%g\n\n', b, skewui10(b), skewvi10(b), skewwi10(b), skewti10(b));
end

fprintf(file_id1, 'NOTE: Skewness needs to be in the range of -2 to 2 (hard flagged) or -1 to 1 (soft flagged). Ideal value of skewness is 0\n\n\n');


%Haar transform
fprintf(file_id1, 'HAAR TRANSFORM\n\n');
fprintf(file_id1, 'no\t HaarU\t HaarV\t HaarW\t HaarT\n\n');

for b = 1:n10
    fprintf(file_id1, '%g\t%g\t%g\t%g\t%g\n\n', b, haarui10(b), haarvi10(b), haarwi10(b), haarti10(b));
end

fprintf(file_id1, 'NOTE: Haar transform must be less than 2 (soft) and 3 (hard), ideal value is 0\n\n\n');

%Turbulence statistics
fprintf(file_id1, 'Turbulence statistics\n\n');
fprintf(file_id1, 'no\t SigU\t\t SigV\t\t SigW\t\t SigT\t\t u*\n\n');

for b = 1:n10
    fprintf(file_id1, '%g\t%6.3g\t%6.3g\t%6.3g\t%6.3g\t%6.3g\n\n', b, SIGMAU10(b), SIGMAV10(b), SIGMAW10(b), SIGMAT10(b), ustar10(b));
end

fprintf(file_id1, 'CONTINUED...\n\n');

fprintf(file_id1, 'no\t wT\t\t ui\t\t vi\t\t V\t\t L\t\t CD\n\n');

for b = 1:n10
    fprintf(file_id1, '%g\t%6.3g\t%6.3g\t%6.3g\t%6.3g\t%6.3g  %6.3g\n\n', b, witi10bar(b), u010bar(b), v010bar(b), mean(V10(:,b)), L10(b),(ustar10(b)/mean(V10(:,b)))^2);
end

fprintf(file_id1, 'CONTINUED...\n\n');

fprintf(file_id1, 'no\t T*\t\t uw\n\n');

for b = 1:n10
    fprintf(file_id1, '%g\t%6.3g\t%6.3g\n\n', b, FRICT10(b), mean(uiwi10(:,b)));
end

%Closing the created text file
fclose(file_id1);


%Creating a text file for 15 min

file_id2 = fopen('ProcessedData_15_min.txt', 'w');

%Writing code must be written in between creating and closing of text file

fprintf(file_id2, 'Date:\n\n'); 
fprintf(file_id2, 'Initial time of measurement:\n\n');
fprintf(file_id2, 'Height of measurement:\n\n');
fprintf(file_id2, 'Number of data in the entire original record: %g\n\n', length(u0));
fprintf(file_id2, 'Number of spikes in the entire record: %g (%g percent)\n\n\n\n',length(u0) - length(u01),(((length(u0) - length(u01))/length(u0))*100));

fprintf(file_id2, 'L = 15 min \n\n');


%kurtosis
fprintf(file_id2, 'KURTOSIS\n\n');
fprintf(file_id2, 'no\t kurtU\t kurtV\t kurtW\t kurtT\n\n');

for b = 1:n15
    fprintf(file_id2, '%g\t%g\t%g\t%g\t%g\n\n', b, kurtui15(b), kurtvi15(b), kurtwi15(b), kurtti15(b));
end

fprintf(file_id2, 'NOTE: Kurtosis needs to be in the range of 1 to 8 (hard flagged) or 2 to 5 (soft flagged). Ideal value of kurtosis is 3\n\n\n');

%skewness
fprintf(file_id2, 'SKEWNESS\n\n');
fprintf(file_id2, 'no\t skewU\t skewV\t skewW\t skewT\n\n');

for b = 1:n15
    fprintf(file_id2, '%g\t%g\t%g\t%g\t%g\n\n', b, skewui15(b), skewvi15(b), skewwi15(b), skewti15(b));
end

fprintf(file_id2, 'NOTE: Skewness needs to be in the range of -2 to 2 (hard flagged) or -1 to 1 (soft flagged). Ideal value of skewness is 0\n\n\n');


%Haar transform
fprintf(file_id2, 'HAAR TRANSFORM\n\n');
fprintf(file_id2, 'no\t HaarU\t HaarV\t HaarW\t HaarT\n\n');

for b = 1:n15
    fprintf(file_id2, '%g\t%g\t%g\t%g\t%g\n\n', b, haarui15(b), haarvi15(b), haarwi15(b), haarti15(b));
end

fprintf(file_id2, 'NOTE: Haar transform must be less than 2 (soft) and 3 (hard), ideal value is 0\n\n\n');

%Turbulence statistics
fprintf(file_id2, 'Turbulence statistics\n\n');
fprintf(file_id2, 'no\t SigU\t\t SigV\t\t SigW\t\t SigT\t\t u*\n\n');

for b = 1:n15
    fprintf(file_id2, '%g\t%6.3g\t%6.3g\t%6.3g\t%6.3g\t%6.3g\n\n', b, SIGMAU15(b), SIGMAV15(b), SIGMAW15(b), SIGMAT15(b), ustar15(b));
end

fprintf(file_id2, 'CONTINUED...\n\n');

fprintf(file_id2, 'no\t wT\t\t ui\t\t vi\t\t V\t\t L\t\t CD\n\n');

for b = 1:n15
    fprintf(file_id2, '%g\t%6.3g\t%6.3g\t%6.3g\t%6.3g\t%6.3g  %6.3g\n\n', b, witi15bar(b), u015bar(b), v015bar(b), mean(V15(:,b)), L15(b),(ustar15(b)/mean(V15(:,b)))^2);
end

fprintf(file_id2, 'CONTINUED...\n\n');

fprintf(file_id2, 'no\t T*\t\t uw\n\n');

for b = 1:n15
    fprintf(file_id2, '%g\t%6.3g\t%6.3g\n\n', b, FRICT15(b), mean(uiwi15(:,b)));
end

%Closing the created text file
fclose(file_id2);


%Creating a text file for 30 min

file_id3 = fopen('ProcessedData_30_min.txt', 'w');

%Writing code must be written in between creating and closing of text file

fprintf(file_id3, 'Date:\n\n'); 
fprintf(file_id3, 'Initial time of measurement:\n\n');
fprintf(file_id3, 'Height of measurement:\n\n');
fprintf(file_id3, 'Number of data in the entire original record: %g\n\n', length(u0));
fprintf(file_id3, 'Number of spikes in the entire record: %g (%g percent)\n\n\n\n',length(u0) - length(u01),(((length(u0) - length(u01))/length(u0))*100));

fprintf(file_id3, 'L = 30 min \n\n');


%kurtosis
fprintf(file_id3, 'KURTOSIS\n\n');
fprintf(file_id3, 'no\t kurtU\t kurtV\t kurtW\t kurtT\n\n');

for b = 1:n30
    fprintf(file_id3, '%g\t%g\t%g\t%g\t%g\n\n', b, kurtui30(b), kurtvi30(b), kurtwi30(b), kurtti30(b));
end

fprintf(file_id3, 'NOTE: Kurtosis needs to be in the range of 1 to 8 (hard flagged) or 2 to 5 (soft flagged). Ideal value of kurtosis is 3\n\n\n');

%skewness
fprintf(file_id3, 'SKEWNESS\n\n');
fprintf(file_id3, 'no\t skewU\t skewV\t skewW\t skewT\n\n');

for b = 1:n30
    fprintf(file_id3, '%g\t%g\t%g\t%g\t%g\n\n', b, skewui30(b), skewvi30(b), skewwi30(b), skewti30(b));
end

fprintf(file_id3, 'NOTE: Skewness needs to be in the range of -2 to 2 (hard flagged) or -1 to 1 (soft flagged). Ideal value of skewness is 0\n\n\n');


%Haar transform
fprintf(file_id3, 'HAAR TRANSFORM\n\n');
fprintf(file_id3, 'no\t HaarU\t HaarV\t HaarW\t HaarT\n\n');

for b = 1:n30
    fprintf(file_id3, '%g\t%g\t%g\t%g\t%g\n\n', b, haarui30(b), haarvi30(b), haarwi30(b), haarti30(b));
end

fprintf(file_id3, 'NOTE: Haar transform must be less than 2 (soft) and 3 (hard), ideal value is 0\n\n\n');

%Turbulence statistics
fprintf(file_id3, 'Turbulence statistics\n\n');
fprintf(file_id3, 'no\t SigU\t\t SigV\t\t SigW\t\t SigT\t\t u*\n\n');

for b = 1:n30
    fprintf(file_id3, '%g\t%6.3g\t%6.3g\t%6.3g\t%6.3g\t%6.3g\n\n', b, SIGMAU30(b), SIGMAV30(b), SIGMAW30(b), SIGMAT30(b), ustar30(b));
end

fprintf(file_id3, 'CONTINUED...\n\n');

fprintf(file_id3, 'no\t wT\t\t ui\t\t vi\t\t V\t\t L\t\t CD\n\n');

for b = 1:n30
    fprintf(file_id3, '%g\t%6.3g\t%6.3g\t%6.3g\t%6.3g\t%6.3g  %6.3g\n\n', b, witi30bar(b), u030bar(b), v030bar(b), mean(V30(:,b)), L30(b),(ustar30(b)/mean(V30(:,b)))^2);
end

fprintf(file_id3, 'CONTINUED...\n\n');

fprintf(file_id3, 'no\t T*\t\t uw\n\n');

for b = 1:n30
    fprintf(file_id3, '%g\t%6.3g\t%6.3g\n\n', b, FRICT30(b), mean(uiwi30(:,b)));
end
%Closing the created text file
fclose(file_id3);


%Creating a text file for 60 min

file_id4 = fopen('ProcessedData_60_min.txt', 'w');

%Writing code must be written in between creating and closing of text file

fprintf(file_id4, 'Date:\n\n'); 
fprintf(file_id4, 'Initial time of measurement:\n\n');
fprintf(file_id4, 'Height of measurement:\n\n');
fprintf(file_id4, 'Number of data in the entire original record: %g\n\n', length(u0));
fprintf(file_id4, 'Number of spikes in the entire record: %g (%g percent)\n\n\n\n',length(u0) - length(u01),(((length(u0) - length(u01))/length(u0))*100));

fprintf(file_id4, 'L = 60 min \n\n');


%kurtosis
fprintf(file_id4, 'KURTOSIS\n\n');
fprintf(file_id4, 'no\t kurtU\t kurtV\t kurtW\t kurtT\n\n');

for b = 1:n60
    fprintf(file_id4, '%g\t%g\t%g\t%g\t%g\n\n', b, kurtui60(b), kurtvi60(b), kurtwi60(b), kurtti60(b));
end

fprintf(file_id4, 'NOTE: Kurtosis needs to be in the range of 1 to 8 (hard flagged) or 2 to 5 (soft flagged). Ideal value of kurtosis is 3\n\n\n');

%skewness
fprintf(file_id4, 'SKEWNESS\n\n');
fprintf(file_id4, 'no\t skewU\t skewV\t skewW\t skewT\n\n');

for b = 1:n60
    fprintf(file_id4, '%g\t%g\t%g\t%g\t%g\n\n', b, skewui60(b), skewvi60(b), skewwi60(b), skewti60(b));
end

fprintf(file_id4, 'NOTE: Skewness needs to be in the range of -2 to 2 (hard flagged) or -1 to 1 (soft flagged). Ideal value of skewness is 0\n\n\n');


%Haar transform
fprintf(file_id4, 'HAAR TRANSFORM\n\n');
fprintf(file_id4, 'no\t HaarU\t HaarV\t HaarW\t HaarT\n\n');

for b = 1:n60
    fprintf(file_id4, '%g\t%g\t%g\t%g\t%g\n\n', b, haarui60(b), haarvi60(b), haarwi60(b), haarti60(b));
end

fprintf(file_id4, 'NOTE: Haar transform must be less than 2 (soft) and 3 (hard), ideal value is 0\n\n\n');

%Turbulence statistics
fprintf(file_id4, 'Turbulence statistics\n\n');
fprintf(file_id4, 'no\t SigU\t\t SigV\t\t SigW\t\t SigT\t\t u*\n\n');

for b = 1:n60
    fprintf(file_id4, '%g\t%6.3g\t%6.3g\t%6.3g\t%6.3g\t%6.3g\n\n', b, SIGMAU60(b), SIGMAV60(b), SIGMAW60(b), SIGMAT60(b), ustar60(b));
end

fprintf(file_id4, 'CONTINUED...\n\n');

fprintf(file_id4, 'no\t wT\t\t ui\t\t vi\t\t V\t\t L\t\t CD\n\n');

for b = 1:n60
    fprintf(file_id4, '%g\t%6.3g\t%6.3g\t%6.3g\t%6.3g\t%6.3g  %6.3g\n\n', b, witi60bar(b), u060bar(b), v060bar(b), mean(V60(:,b)), L60(b),(ustar60(b)/mean(V60(:,b)))^2);
end

fprintf(file_id4, 'CONTINUED...\n\n');

fprintf(file_id4, 'no\t T*\t\t uw\n\n');

for b = 1:n60
    fprintf(file_id4, '%g\t%6.3g\t%6.3g\n\n', b, FRICT60(b), mean(uiwi60(:,b)));
end

%Closing the created text file
fclose(file_id4);

%Creating a text file for 120 min

file_id5 = fopen('ProcessedData_120_min.txt', 'w');

%Writing code must be written in between creating and closing of text file

fprintf(file_id5, 'Date:\n\n'); 
fprintf(file_id5, 'Initial time of measurement:\n\n');
fprintf(file_id5, 'Height of measurement:\n\n');
fprintf(file_id5, 'Number of data in the entire original record: %g\n\n', length(u0));
fprintf(file_id5, 'Number of spikes in the entire record: %g (%g percent)\n\n\n\n',length(u0) - length(u01),(((length(u0) - length(u01))/length(u0))*100));

fprintf(file_id5, 'L = 120 min \n\n');


%kurtosis
fprintf(file_id5, 'KURTOSIS\n\n');
fprintf(file_id5, 'no\t kurtU\t kurtV\t kurtW\t kurtT\n\n');

for b = 1:n120
    fprintf(file_id5, '%g\t%g\t%g\t%g\t%g\n\n', b, kurtui120(b), kurtvi120(b), kurtwi120(b), kurtti120(b));
end

fprintf(file_id5, 'NOTE: Kurtosis needs to be in the range of 1 to 8 (hard flagged) or 2 to 5 (soft flagged). Ideal value of kurtosis is 3\n\n\n');

%skewness
fprintf(file_id5, 'SKEWNESS\n\n');
fprintf(file_id5, 'no\t skewU\t skewV\t skewW\t skewT\n\n');

for b = 1:n120
    fprintf(file_id5, '%g\t%g\t%g\t%g\t%g\n\n', b, skewui120(b), skewvi120(b), skewwi120(b), skewti120(b));
end

fprintf(file_id5, 'NOTE: Skewness needs to be in the range of -2 to 2 (hard flagged) or -1 to 1 (soft flagged). Ideal value of skewness is 0\n\n\n');


%Haar transform
fprintf(file_id5, 'HAAR TRANSFORM\n\n');
fprintf(file_id5, 'no\t HaarU\t HaarV\t HaarW\t HaarT\n\n');

for b = 1:n120
    fprintf(file_id5, '%g\t%g\t%g\t%g\t%g\n\n', b, haarui120(b), haarvi120(b), haarwi120(b), haarti120(b));
end

fprintf(file_id5, 'NOTE: Haar transform must be less than 2 (soft) and 3 (hard), ideal value is 0\n\n\n');

%Turbulence statistics
fprintf(file_id5, 'Turbulence statistics\n\n');
fprintf(file_id5, 'no\t SigU\t\t SigV\t\t SigW\t\t SigT\t\t u*\n\n');

for b = 1:n120
    fprintf(file_id5, '%g\t%6.3g\t%6.3g\t%6.3g\t%6.3g\t%6.3g\n\n', b, SIGMAU120(b), SIGMAV120(b), SIGMAW120(b), SIGMAT120(b), ustar120(b));
end

fprintf(file_id5, 'CONTINUED...\n\n');

fprintf(file_id5, 'no\t wT\t\t ui\t\t vi\t\t V\t\t L\t\t CD\n\n');

for b = 1:n120
    fprintf(file_id5, '%g\t%6.3g\t%6.3g\t%6.3g\t%6.3g\t%6.3g  %6.3g\n\n', b, witi120bar(b), u0120bar(b), v0120bar(b), mean(V120(:,b)), L120(b),(ustar120(b)/mean(V120(:,b)))^2);
end

fprintf(file_id5, 'CONTINUED...\n\n');

fprintf(file_id5, 'no\t T*\t\t uw\n\n');

for b = 1:n120
    fprintf(file_id5, '%g\t%6.3g\t%6.3g\n\n', b, FRICT120(b), mean(uiwi120(:,b)));
end
%Closing the created text file
fclose(file_id5);

toc;    %End run timer