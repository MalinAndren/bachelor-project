close all 
clear all
soot_raw=dlmread(['pastedsootrawfile_2013.dat'], '', 0, 0);
soot_raw=soot_raw(1:size(soot_raw, 1), :);
datum=soot_raw(:, 1);

datum=num2str(datum);

year=str2num(datum(:, 1:4));
mon=str2num(datum(:, 5:6));
dag=str2num(datum(:, 7:8));
tid=floor(soot_raw(:, 2));

tid=1e6+tid;
tid_raw=num2str(tid);
tid=num2str(tid);

hour=str2num(tid(:, 2:3));
min=str2num(tid(:, 4:5));
sek=str2num(tid(:, 6:7));

I0=soot_raw(:, 3);
I=soot_raw(:, 4);
qobs=soot_raw(:, 5);
filter_change=soot_raw(:, 11);
juldate=datenum(year, mon, dag, hour, min, sek)-datenum(year-1, 12, 31, 0, 0, 0);

diff_in_I=diff(I);


diff_crit=diff_in_I>0.3;

step_qual=0;
aa=0;
invalid=0

period=60
timespan=60

kk=0;

for step=period:1:(size(I, 1)-period)
    kk=kk+1
    Iout(kk, 1:5)=[juldate(step), median(I((step-period/2):(step+period/2))), median(I0((step-period/2):(step+period/2))), diff_crit(step), qobs(step)];
end
%slutar här
%below->the actual calculation of absorption using corr for loading only.

filterarea=pi*(1.5e-3)^2;
filterconst=2
orgratio=0.83

indo=0;
out_inter=[];
flow_tot=0;
Out_absorption=[]
dd=0

for a=1:size(Iout, 1)
    indo=indo+1
    flow_tot=flow_tot+1e-3*Iout(a, 5);
        if Iout(a, 4)==1
        indo=0;
        out_inter=[];
        end
        
       
            
            
            if indo==1 %first value in selected timeseries
         I_out_inter(1)=Iout(a, 2);
         I0_out_inter(1)=Iout(a, 3);
        end
    
        if indo==floor(timespan/2) %%%%Take the center value of I and I0
                                   ...during the intergration period to 
                                   ...estimate loading
                                   
              I_center(1:2)=[Iout(a, 2), Iout(a, 3)]   
        end
        
            
        if indo==timespan %last value in selected timeseries
        dd=dd+1 %index for outputmatrix "Out_absorption"
          
          I_out_inter(2)=Iout(a, 2);
          I0_out_inter(2)=Iout(a, 3);
            
          %%CALCULATES THE log((I0(1)/I(1))/(I0(last)/I(last))))
          xa=I0_out_inter(1)/I_out_inter(1)
          xb=I0_out_inter(2)/I_out_inter(2)
          lndata=log(xb/xa)
          %%%%%%%%%%%%%%%%%%%%%%%%%%%
          
        %%%CALCULATION FOR LOADIUNG CORRECTION BELOW%%%%%%%  
        loading=(I_center(1)/I_center(2))/orgratio;
        R=(0.355+0.5398*loading);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%CALCULATION FOR LOADIUNG CORRECTED ABSORBTION COEFFICIENT BELWO%%%%%%%  
        babs_uncorrected=lndata*filterarea/(flow_tot)/filterconst %m-1
        babs_corrected=lndata*filterarea/(flow_tot)/filterconst/R %m-1
         
   
      %  plot(Iout(a, 1), babs_corrected,'.r', Iout(a, 1), babs_uncorrected, '.g')
        
        %Below is the outputmatrix. Column 1 is decimal day of year, rest
        %should be intuative
        Out_absorption(dd, 1:5)=[Iout(a, 1),babs_corrected, babs_uncorrected, loading, flow_tot];
        hold on
  
        
        flow_tot=0; %resets the total flow when a new period starts
        
          indo=0;
        end
 
end
 
 
           %absorbtion data plotted below
          plot(Out_absorption(:, 1), Out_absorption(:, 2), 'r',Out_absorption(:, 1), Out_absorption(:, 3), 'g' )
          legend('corrected', 'uncorrected')
          
          ylabel('m^{-1}')
          xlabel('decimal day of year')
          datan=[];
          doy_int=floor(Out_absorption(:, 1));
          d_out=[];
         close all
          for a=350:370
              
              aa=find(doy_int==a);
              flag=[]
              datan=[]
       
              if ~isempty(aa)
                  
                  datan=Out_absorption(aa, 1:5);
                  
            
         plot(datan(:, 1), datan(:, 2), 'r',datan(:, 1), datan(:, 3), 'g' )
        flag(1:size(datan, 1), 1)=0;
              end
 
        button=0
         
         if ~isempty(datan)
             close all
         while button<3
        plot(datan(:, 1), datan(:, 2), 'r.-',datan(:, 1), datan(:, 3), 'g.-' )
        
        
            [x, y, button]=ginput(2)
            
            if ~isempty(datan(:,1)>x(1) & datan(:, 1)<x(2))
            
            flag(datan(:,1)>x(1) & datan(:, 1)<x(2), 1)=0.999 
         
            end
            
                % datan(flag>0, 2:5)=NaN;
                 
                 datan2=[datan, flag];
                 plot(datan2(datan2(:, 6)==0.999,1),datan2(datan2(:, 6)==0.999,3), 'k+', 'markersize', 20, 'linewidth', 3)
 
         hold on
         datan_out=[datan, flag]
 
        
        
       
         end
        
         d_out(size(d_out, 1)+1:size(d_out, 1)+size(datan_out, 1), 1:6)=datan_out;
         end
          end
          
save screened_data2011.dat d_out -ascii


    
    

