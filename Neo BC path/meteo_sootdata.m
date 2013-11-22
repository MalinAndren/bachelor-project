close all 
clear all


daymat=[31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
meteodata1=dlmread(['meteo.dat'], '', 0, 0);

%meteodata1=meteodata1(meteodata1(:, size(meteodata1, 2))== -999, :);

year=meteodata1(:, 3);
month=meteodata1(:, 2);
day=meteodata1(:, 1);
hour=meteodata1(:, 4);
min=meteodata1(:, 5);

juldate=datenum(year, month, day, hour, min, 0)-datenum(year-1, 12, 31, 0, 0, 0);

wind_speed=meteodata1(:, 7);
wind_dir=meteodata1(:, 8);
air_temp=meteodata1(:, 12);
humidity=meteodata1(:, 13);
dew_point=meteodata1(:, 14);


soot2011_13=load('soot11_13.dat');

soot2011_13=soot2011_13(soot2011_13(:, size(soot2011_13, 2))== 0 | soot2011_13(:, size(soot2011_13, 2))== 0.777, :);

year_s=soot2011_13(:, 1);
month_s=soot2011_13(:,2);
day_s=soot2011_13(:, 3);

hour_s=soot2011_13(:, 4);
min_s=soot2011_13(:, 5);
juldate=datenum(year, month, day, hour, min, 0)-datenum(year-1, 12, 31, 0, 0, 0);

abs=soot2011_13(:, 6);




%timmes samt månadsvariation/mean vindstyrka nedan

hourly_windspeed=[];


for timme=0:23
   
   ind=timme+1
   hourly_windspeed(ind, 1:3)=prctile(wind_speed(hour==timme), [25, 50, 75])
end
figure(10+ind) 
errorbar(1:24, hourly_windspeed(:, 2), hourly_windspeed(:, 2)-hourly_windspeed(:, 1), hourly_windspeed(:, 3) -hourly_windspeed(:, 2), 'r')

title('HOURLY VARIATION IN WINDSPEED') 
xlabel('Variations m s^-1')
ylabel('Time')


monthly_windspeed=[]
for month2=1:12
     
    if sum(month==month2)>0
        monthly_windspeed(month2, 1:3)=prctile(wind_speed(month==month2), [25, 50, 75])
    else
        monthly_windspeed(month2, 1:3)=NaN
   
    end
end
figure(1+ind)
errorbar(1:12, monthly_windspeed(:, 2), monthly_windspeed(:, 2)-monthly_windspeed(:, 1), monthly_windspeed(:, 3)-monthly_windspeed(:, 2), 'k') 
title('MONTHLY VARIATION IN WINDSPEED')
xlabel('month')
ylabel('Wind speed m s^{-1}')


hourly_winddir=[]

for timme=0:23
    ind=timme+1
    hourly_winddir(ind, 1:3)=prctile(wind_dir(hour==timme), [25, 50, 75]);
end
figure(10+ind)
errorbar(1:24, hourly_winddir(:, 1), hourly_winddir(:, 2)-hourly_winddir(:, 1), hourly_winddir(:, 3)-hourly_winddir(:, 2), 'r')
title('HOURLY VARIATION IN WIND DIRECION')
xlabel('time')
ylabel('direction in degrees')



monthly_winddir=[]

for month2=1:12
if sum (month==month2)>0
    monthly_winddir(month2, 1:3)=prctile(wind_dir(month==month2), [25, 50, 75])
else
    monthly_winddir(month2, 1:3)=NaN
end
end
figure(1+ind)
errorbar(1:12, monthly_winddir(:, 1), monthly_winddir(:, 2)-monthly_winddir(:, 1), monthly_winddir(:, 3)-monthly_winddir(:, 2), 'k')
title('MONTHLY VARIATION IN WIND DIRECTION')
xlabel('month')
ylabel('directions in degrees')

%Timmes/månadsvariation temperatur

hourly_temp=[];
    for timme=0:23
        ind=timme+1
        hourly_temp(ind, 1:3)=prctile(air_temp(hour==timme), [25, 50, 75])
    end
    
    figure(10+ind)
   
    errorbar(1:24, hourly_temp(:, 2), hourly_temp(:, 2)-hourly_temp(:, 1), hourly_temp(:, 3)-hourly_temp(:, 2), 'r')
    
    title('HOURLY VARIATION IN TEMPERATURE') 
    xlabel('time')
    ylabel('temperature celcius')

    
    monthly_temp=[];

for month2=1:12
    if sum(month==month2)>0 
    
         monthly_temp(month2, 1:3)=prctile(air_temp(month==month2), [25, 50, 75])
    else
        monthly_temp(month2, 1:3)=NaN
       
   
     end 
end
figure(1+ind)
errorbar(1:12,...
    monthly_temp(:, 2),...
    monthly_temp(:, 2) - monthly_temp(:, 1),...
    monthly_temp(:, 3)-monthly_temp(:, 2), 'k')
title('MONTHLY VARIATION IN TEMPERATURE')
xlabel('month')
ylabel('temperature celcius')



%timme samt månads mean fuktighet här under
hourly_humidity=[];
for timme=0:23
    
    ind=timme+1
    hourly_humidity(ind, 1:3)=prctile(humidity(hour==timme), [25, 50, 75])
end
figure(10+ind)
errorbar(1:24, hourly_humidity(:, 2), hourly_humidity(:, 2)-hourly_humidity(:, 1), hourly_humidity(:, 3)-hourly_humidity(:, 2), 'r')

title('HOURLY VARIATION IN HUMIDITY')
xlabel('time')
ylabel('humidity %')

%semilogx(x,y)
%semilogy(x,y)
%loglog(x, y)
%set(gca, 'xlim', [0, 23], 'ylim', [40, 100])
%xlabel('Hour of day')

monthly_humidity=[];
for month2=1:12
if sum (month==month2)>0
    monthly_humidity(month2, 1:3)=prctile(humidity(month==month2),[25, 50, 75])
else
    monthly_humidity(month2, 1:3)=NaN
end 
end 
errorbar(1:12, monthly_humidity(:, 2), monthly_humidity(:, 2)-monthly_humidity(:, 1), monthly_humidity(:, 3)-monthly_humidity(:, 2), 'k')
title('MONTHLY VARIATION IN HUMIDITY')
xlabel('month')
ylabel('humidity %')




%timme samt månads mean daggpunkt nedan

hourly_dewpoint=[];
for timme=0:23
    ind=timme+1
    hourly_dewpoint(ind, 1:3)= prctile(dew_point(hour==timme), [25, 50, 75])
end
figure(10+ind) 

errorbar(1:24, hourly_dewpoint(:, 2),  hourly_dewpoint(:, 2)- hourly_dewpoint(:, 1), hourly_dewpoint(:, 3)- hourly_dewpoint(:, 2), 'r')
title('HOURLY VARIATION IN DEW POINT') 
xlabel('time')
ylabel('dew point')



monthly_dewpoint=[];
for month2=1:12
    if sum (month==month2)>0 
        monthly_dewpoint(month2, 1:3)=prctile(dew_point(month==month2), [25, 50, 75])
    else
        monthly_dewpoint(month2, 1:3)=NaN
    end 
end
errorbar(1:12, monthly_dewpoint(:, 2), monthly_dewpoint(:, 2)-monthly_dewpoint(:, 1), monthly_dewpoint(:, 3)-monthly_dewpoint(:, 2), 'k')
title('MONTHLY VARIATION IN DEW POINT')
xlabel('month')
ylabel('dew point')




%PERCENTIL/MEDIAN TIMME SOT

hourly_abs=[]
for timme=0:23
    ind=timme+1
    hourly_abs(ind, 1:3)=prctile(abs(hour_s==timme), [25, 50, 75])
    
    end
figure(10+ind)
errorbar(1:24, hourly_abs(:, 2), hourly_abs(:, 2)-hourly_abs(:, 1), hourly_abs(:, 3)-hourly_abs(:, 2), 'r')
title('HOURLY SOOT ABSORBTION')
xlabel('time')
ylabel('Absorption m^{-1}')




    
    %MÅNADS PECENTIL/MEDIAN SOT
    
    monthly_abs=[];
    for month2=2:12
        if sum(month_s==month2)>0
            monthly_abs(month2, 1:3)=prctile(abs(month_s==month2), [25, 50, 75])
        else
            monthly_abs(month2, 1:3)=NaN
        end 
    end
    errorbar(1:12, monthly_abs(:, 2), monthly_abs(:, 2)-monthly_abs(:, 1), monthly_abs(:, 3)-monthly_abs(:, 2), 'r')
    title('MONTHLY ABSORPTION')
    xlabel('decimal month of year')
    ylabel('absorption m^{-1}')





%PERCENTIL/MEDIAN VECKA SOT

weekly_abs=[]
 for days_week=1:7
     if sum (day_s==days_week)>0
         weekly_abs(days_week, 1:3)=prctile(abs(day_s==days_week), [25, 50, 75])
     else
         monthly_abs(days_week,1:3)=NaN
     end 
 end
 
errorbar(1:7,...
    weekly_abs(:, 2), ...
    weekly_abs(:, 2) - weekly_abs(:, 1),...
    weekly_abs(:, 3)-weekly_abs(:, 2), 'g')
title('WEEKLY ABSORPTION')
xlabel('decimal day of week')
ylabel('absorption')

%korrelation mellan abs_hourly och temperatur MÅNADSVIS

figure(88)

 errorbar(1:12, monthly_abs(:, 2), monthly_abs(:, 2)-monthly_abs(:, 1), monthly_abs(:, 3)-monthly_abs(:, 2), 'r')
 legend('abs')
 

   


 
 
 
%KORRELATION MELLAN TEMPERATUR OCH SOT MÅNADSVIS
 




[haxes, hline1, hline2]=plotyy(1:12, monthly_temp(:, 2), 1:12, monthly_abs(:, 2),'plot')

xlabel('month')

axes(haxes(1))
ylabel('temperature')

axes(haxes(2))
ylabel('absorption')

title('correlation between temperature and absorption') 


%% Korrelation mellan sot/temp/daggpunkt månadsvis

[ax,hlines] = plotyyy(1:12, monthly_temp(:, 2), 1:12, monthly_abs(:, 2), 1:12, monthly_dewpoint(:, 2), ylabels)

xlabel('decimal moth of year')

ylabels{1}='temp';
ylabels{2}='abs';
ylabels{3}='dewpoint';
title('correlation between temp/abs/dewpoint')


[ax, hlines] = plotyyy(1:12, monthly_abs(:, 2), 1:12, monthly_winddir(:, 2), 1:12, monthly_temp(:, 2), ylabels)

xlabel('decimal month of year')
ylabels{1}='abs';
ylabels{2}='wind direction degrees';
ylabels{3}='temp';


%KORRELATION MELLAN WINDRIKTNING; TEMP OCH SOT

[ax, hlines] = plotyyy(1:24, hourly_abs(:, 2),1:24, hourly_winddir (:, 2), 1:24, hourly_temp(:, 2), ylabels)

xlabel=('time')
ylabels{1}='abs';
ylabels{2}= 'wind direction';
ylabels{3}= 'temperature';
title('hourly correaltion between soot/temp/wind direction')


%% Korrelation timvis mellan vindriktning och absorption.



[haxes, hline1, hline2]=plotyy(1:24, hourly_winddir(:, 2), 1:24, hourly_abs(:, 2), 'plot')

xlabel('time')

axes(haxes(1))
ylabel('wind direction')

axes(haxes(2))
ylabel('absorption')


title('hourly correlation between wind direction and absorption') 

%% Korrelation månadsvis mellan vindriktning och absorption

[haxes, hline1, hline2,]=plotyy(1:12, monthly_winddir(:, 2), 1:12, monthly_abs(:, 2), 'plot')
xlabel('decimal month of year')

axes(haxes(1))

ylabel('wind direction')

axes(haxes(2))
ylabel('absorption')


title('correlation between wind direction and absorption NEO')

%% Korrelation timvis mellan temperatur och absorption 
[haxes, hline1, hline2]=plotyy(1:24, hourly_temp(:, 2), 1:24, hourly_abs(:, 2), 'plot')

xlabel('time')

axes(haxes(1))
ylabel('temperature')

axes(haxes(2))
ylabel('absorption')

title('hourly correlation between temperature and absorption')

%% Korrealtion månadsvis mellan temperatur och sot 

[haxes, hline1, hline2]=plotyy(1:24, monthly_temp(:, 2),1:24, monthly_abs(:, 2), 'plot')

xlabel('decimal month of year')
axes(haxes(1))
ylabel('temperature')

axes(haxes(2))
ylabel('absorption')

title('monthly correlation between temperature and absorption')

%% Korrelation mellan daggpunkt och soot

[haxes, hline1, hline2]=plotyy(1:12, monthly_dewpoint(:, 2), 1:12, monthly_abs(:, 2), 'plot')

xlabel('decimal month of year')

axes(haxes(1))
ylabel('dewpoint')

axes(haxes(2))
ylabel('absorption')

title('monthly correlation between soot and dewpoint NEO')
%% 

[haxes, hline1, hline2]=plotyy(1:24, hourly_dewpoint(:, 2), 1:24, hourly_abs(:, 2), 'plot')

xlabel('time')

axes(haxes(1))
ylabel('dew point')

axes(haxes(2))
ylabel('absorption')

title('hourly correlation between dewpoint and absorption NEO')

%% korelation mellan abs, wind speed and wind direction 
%%h = figure(3)
[ax,hlines] = plotyyy(1:24, hourly_abs(:, 2), 1:24, hourly_windspeed(:, 2), 1:24, hourly_winddir(:, 2), ylabels)

xlabel('decimal month of year')

ylabels{1}='absorption m{^-1}';
ylabels{2}='windspeed';
ylabels{3}='wind direction';
title('correlation between abs/windspeed/wind direction')

%%saveas(h,'C:\Users\itm\Desktop\NEO_BC\fig\blub.png')
%%


t=hourly_winddir(:, 2)

polar(t,sin(2*t).*cos(2*t), '--go');


hold on 

t2=hourly_windspeed(:, 2)
polar(t2,sin(2*t2).*cos(2*t2), '--ro');
legend('windspeed')

hold on

sortrows

t3=hourly_abs(:, 2)
polar(t3,sin(2*t3).*cos(2*t3), 'b')


%%% BINDATA

x=hourly_abs(:, 2);
y=hourly_winddir(:,1);

TopEdge=1;
BotEdge=0;
nbin=2;

BinEdges=linspace(BotEdge, TopEdge, nbin+1)

[h, WhichBin]=histc(x, BinEdges);

for i=1:nbin
    flagBinMember = (WhichBin==i);
    BinMember = y(flagBinMembers);
    BinMean(i)= mean(BinMember);
end





%timmedelvärde för sot och vindriktning 
ftp=0

for year1=2011
    for month1=1:12
        for day1=1:daymat(month1)
            for hour1=0:23
                [a]=find(year_s==year1 & month_s==month1 & day_s==day1 & hour_s==hour1);
                [b]=find(year==year1 & month==month1 & day==day1 & hour==hour1);
                 if ~isempty(a) & ~isempty (b)
                     ftp=ftp+1
                flag=0
              outmatrix(ftp, 1:7)=[year1, month1, day1, hour1, median(abs(a, :)), median(wind_dir(b, :)), flag];
                 elseif ~isempty(a) & isempty(b)
                     flag=0.666
                     ftp=ftp+1
                     outmatrix(ftp, 1:7)=[year1, month1, day1, hour1, median(abs(a, :)), -999, flag];
                     
                 elseif isempty(a) & ~isempty(b)
                     
                    flag=0.777
                     ftp=ftp+1
                     outmatrix(ftp, 1:7)=[year1, month1, day1, hour1, -999, -999, flag];
                     
                
                outmatrix(outmatrix(:, 7)==0, :);
                P=outmatrix(outmatrix(:, 7)==0, :);
                plot(P(:, 6), P(:, 5), 'xr')
                xlabel('wind')
                ylabel('soot')
                 end
            end
        end
    end
end
 
%timmesmedelvärde sot och windhastighet
kkl=0
for year2=2011
    for month2=1:12
        for day2=1:daymat(month2)
            for hour2=0:23
                [a]=find(year_s==year2 & month_s==month2 & day_s==day2 & hour_s==hour2);
                [b]=find(year==year & month==month2 & day==day2 & hour==hour2);
                 if ~isempty(a) & ~isempty (b)
                     flag=0
                
                     kkl=kkl+1
                     outmatrix(kkl, 1:7)=[year2, month2, day2, hour2, median(abs(a, :)), median(wind_speed(b, :)), flag];
                     
                 elseif ~isempty(a) & isempty(b)
                     flag=0.666
                     kkl=kkl+1
                     
                    outmatrix(kkl, 1:7)=[year2, month2, day2, hour2, median(abs(a, :)), -999, flag];
                    
                 elseif isempty(a) & ~isempty(b) 
                     flag=0.777
                     kkl=kkl+1
                     
                     outmatrix(kkl, 1:7)=[year2, month2, day2, hour2, -999, -999, flag];
                       
                outmatrix(outmatrix(:, 7)==0, :);
                P=outmatrix(outmatrix(:, 7)==0, :);
                plot(P(:, 6), P(:, 5), 'xr')
                xlabel('wind')
                ylabel('soot')
                xlim('0:100')
                
               
                 end
            end
        end
    end
end

                 
                
               
              
%timmedel för sot och temp

kk=0
for year3=2011
    for month3=1:12
        for day3=1:daymat(month3)
        for hour3=0:23
            [a]=find(year_s==year3 & month_s==month3 & day_s==day3 & hour_s==hour3);
            [b]=find(year==year3 & month==month3 & day==day3 & hour==hour3);
            if ~isempty(a) & ~isempty(b)
                flag=0
                kk=kk+1
                outmatrix(kk, 1:7)=[year3, month3, day3, hour3, median(abs(a, :)), median(air_temp(b, :)), flag];
            elseif ~isempty(a) & isempty (b)
                flag=0.666
                kk=kk+1
                outmatrix(kk, 1:7)=[year3, month3, day3, hour3, median(abs(a, :)), -999, flag];
                
            elseif isempty (a) & ~isempty (b) 
                
                flag=0.888
                
                kk=kk+1
                outmatrix(kk, 1:7)=[year3, month3, day3, hour3, -999, -999, flag];
                
                outmatrix(outmatrix(:, 7)==0, :);
                P=outmatrix(outmatrix(:, 7)==0, :);
                plot(P(:, 6), P(:, 5), 'xr')
                xlabel('temp')
                ylabel('soot')
                set(gca, 'xlim', [0, 60])
            end
           end
        end
    end
end

ftw=0
for wd=0:5:365
    [a]=find(wind_dir>=wd)
    [b]=find(wind_dir<wd) %& wind_dir(:, size(wind_dir))==0); 
    if~isempty (a) & ~isempty(b)
    ftw=ftw+1;
    outmatrix(ftw, 1:4)=wd
    elseif is empty(a) & ~isempty(b)
        flag=0.666
        ftw=ftw+1
    else if ~isempty(a) & isempty(b)
            flag=0.777
            ftw=ftw+1
            
        
    
    
   % out_matrix(ftw, 1:2)=[wd+10/2]
  
    end
end

  y=sin(wind_dir)              
           
figure 
t=y
    polar(t, sin(2*t).*cos(2*t), 'or')


     
            


polar








     