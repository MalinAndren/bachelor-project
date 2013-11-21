close all
clear all
%A=A(A(:, size(A, 2))==0, :)


soot2013=load('correctedsootdata2013.dat')

soot2013=soot2013(soot2013(:, size(soot2013, 2))== 0, :);

year=soot2013(:, 1);
mon=soot2013(:,2);
day=soot2013(:, 3);

hour=soot2013(:, 4);
min=soot2013(:, 5);
juldate=datenum(year, mon, day, hour, min, 0)-datenum(year-1, 12, 31, 0, 0, 0);

abs=soot2013(:, 6);

plot(juldate, soot2013(:, 6), 'g');
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
title('hourly mean 2013 NEO')
xlabel('hour')
ylabel('absorption')




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
    title('monthly mean 2013 NEO')
    xlabel('month')
    ylabel('absorption')





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
title('weekly mean 2013 NEO')
ylabel('absorption')
xlabel('day of week')
