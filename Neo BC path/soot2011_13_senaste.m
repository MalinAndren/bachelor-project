close all
clear all
%A=A(A(:, size(A, 2))==0, :)


soot2011_13=load('correctedsootdata2011_13.dat')

soot2011_13=soot2011_13(soot2011_13(:, size(soot2011_13, 2))== 0 | soot2011_13(:, size(soot2011_13, 2))== 0.777, :);

year=soot2011_13(:, 1);
mon=soot2011_13(:,2);
day=soot2011_13(:, 3);

hour=soot2011_13(:, 4);
min=soot2011_13(:, 5);
juldate=datenum(year, mon, day, hour, min, 0)-datenum(year-1, 12, 31, 0, 0, 0);

abs=soot2011_13(:, 6);

plot(juldate, soot2011_13(:, 6), 'g');
hold on
%close all
%clear all

%flaggad2012=dlmread('screened_data2013.dat')
%plot(flaggad2012(:,1), flaggad2012(:,3), 'b')
%hold on
%plot(flaggad2012(flaggad2012(:,6)==0.999,1), flaggad2012(flaggad2012(:,6)==0.999,3), '.k')


 

ylabel('m^{-1}')
xlabel('decimal day of year')
    
%timmedelvärde    
hourly_abs=[];
    for timme=0:23
        ind=timme+1
        hourly_abs(ind)=mean(abs(hour==timme))
    end
    
    figure(10+ind)
    plot(1:24, hourly_abs, 'g')
title('hourly mean 2011 NEO')
xlabel('hour')
ylabel('absorption')

%PERCENTIL/MEDIAN

hourly_abs=[]
for timme=0:23
    ind=timme+1
    hourly_abs(ind, 1:3)=prctile(abs(hour==timme), [25, 50, 75])
    
    end
figure(10+ind)
errorbar(1:24, hourly_abs(:, 2), hourly_abs(:, 2)-hourly_abs(:, 1), hourly_abs(:, 3)-hourly_abs(:, 2), 'r')
title('HOURLY ABSORBTION')
xlabel('time')
ylabel('Absorption m^{-1}')




%månads medel-absorption under 
month_abs=[]
for month2=1:12
    if sum(mon==month2)>0
        month_abs(month2)=mean(abs(mon==month2))
    else month_abs(month2)=NaN
    end
end

figure(20+ind)
plot(1:12, month_abs, 'k')
    title('monthly mean 2011 NEO')
    xlabel('month')
    ylabel('absorption')
    
    %MÅNADS PECENTIL/MEDIAN
    
    monthly_abs=[];
    for month2=2:12
        if sum(mon==month2)>0
            monthly_abs(month2, 1:3)=prctile(abs(mon==month2), [25, 50, 75])
        else
            monthly_abs(month2, 1:3)=NaN
        end 
    end
    errorbar(1:12, monthly_abs(:, 2), monthly_abs(:, 2)-monthly_abs(:, 1), monthly_abs(:, 3)-monthly_abs(:, 2), 'r')
    title('MONTHLY ABSORPTION')
    xlabel('decimal month of year')
    ylabel('absorption m^{-1}')





%veckomedelvärde?

weekly_abs=[];

for days_week=1:7
    if sum(day==days_week)>0
        weekly_abs(days_week)=mean(abs(day==days_week))
    else weekly_abs(days_week)=NaN
    end
    %weekly_absorption=mean(abs(day==days_week))
end

plot(1:7, weekly_abs, 'r')
title('weekly mean 2011 NEO')
ylabel('absorption')
xlabel('day of week')


%PERCENTIL/MEDIAN VECKA

weekly_abs=[]
 for days_week=1:7
     if sum (day==days_week)>0
         weekly_abs(days_week, 1:3)=prctile(abs(day==days_week), [25, 50, 75])
     else
         monthly_abs(days_week,1:3)=NaN
     end 
 end
 
errorbar(1:7, weekly_abs(:, 2), weekly_abs(:, 2)-weekly_abs(:, 1), weekly_abs(:, 3)-weekly_abs(:, 2), 'g')
title('WEEKLY ABSORPTION')
xlabel('decimal day of week')
ylabel('absorption')





