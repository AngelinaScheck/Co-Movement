
%%Manual time trace analysis of hiFRET DNA stds in vivo.  Data from
%%selections of movies taken from Fiji and just graphed here
%%RC 6/7/12
clear all; close all;

ALEX = 0;
firstframeG = 1;
timeres = 0.010;
frame_after_Gbleach = 1791;
startframe = 1;
frame_after_Rbleach = 1791;



if (ALEX == 0)

    filename_DA = uigetfile('*.txt','Select DA data');
    filename_DD = uigetfile('*.txt','Select DD data');
    %BKGD_DA = uigetfile('*.*','Select BKGD DA data');
    %BKGD_DD = uigetfile('*.*','Select BKGD DD data');

    DA_data = importdata(filename_DA);
    frames = DA_data.data(:,6);
    DA_data = DA_data.data(:,2);
    DD_data = importdata(filename_DD);
    DD_data = DD_data.data(:,2);
    
    endframe = length(DA_data);
%     BG_DA = importdata(BKGD_DA);
%     BG_DA = BG_DA.data(:,2);
%     BG_DD = importdata(BKGD_DD);
%     BG_DD = BG_DD.data(:,2);

    starttime = startframe*timeres;
    times = frames.*timeres;
    
    xmax = max(times)+1;

    BG_DA_mean = mean(DA_data(frame_after_Gbleach:endframe));
    BG_DD_mean = mean(DD_data(frame_after_Gbleach:endframe));
    

    %subplot(3,1,1,'Fontsize',20)
    rect = [500, 500, 500, 350]
    figure('OuterPosition',rect)
    subplot(3,1,[2 3])
    plot(times(startframe:endframe),DA_data(startframe:endframe),'r-', 'LineWidth',2);
    hold on
    plot(times(startframe:endframe),DD_data(startframe:endframe),'g-', 'LineWidth',2);xlim([starttime xmax])
    hold on

    %bkgd subtract
    DA_data_BS = DA_data - BG_DA_mean;
    DD_data_BS = DD_data - BG_DD_mean;

    %subplot(3,1,2,'Fontsize',20)
    %plot(time,DA_data_BS,'r-', 'LineWidth',2);
    %hold on
    %plot(time,DD_data_BS,'g-', 'LineWidth',2);xlim('auto'); title('BKGD subtracted','Fontsize',20);xlabel('','Fontsize',20), ylabel('Intensity','Fontsize',20);

    %calculate FRET
    E = DA_data_BS./(DA_data_BS + DD_data_BS);
    %subplot(3,1,3,'Fontsize',20)
    subplot(3,1,1)
    plot(times(startframe:frame_after_Rbleach-1),E(startframe:frame_after_Rbleach-1),'b-', 'LineWidth',2);ylim([0 1]);xlim([starttime xmax]);
    E_b4_bleach = E(startframe:frame_after_Rbleach-1);
   
    save([filename_DA 'Edata.mat']);
    mean(E_b4_bleach)
    
else
    
    
    filename_Aem = uigetfile('*.*','Select Aem data');
    filename_Dem = uigetfile('*.*','Select Dem data');
    BKGD_R = uigetfile('*.*','Select BKGD Aem data');
    BKGD_G = uigetfile('*.*','Select BKGD Dem data');

    %import data and sort out ALEX frames
    Aem_data = importdata(filename_Aem);
    frames = Aem_data.data(:,6);
    Aem_data = Aem_data.data(:,2);
    Dem_data = importdata(filename_Dem);
    Dem_data = Dem_data.data(:,2);
    BG_Aem = importdata(BKGD_R);
    BG_Aem = BG_Aem.data(:,2);
    BG_Dem = importdata(BKGD_G);
    BG_Dem = BG_Dem.data(:,2);
    
    times = frames.*timeres;

    if (firstframeG  == 1);
        
        DA_data = Aem_data(1:2:end);
        AA_data = Aem_data(2:2:end);
        DD_data = Dem_data(1:2:end);
        
        BG_DA = BG_Aem(1:2:end); BG_DA_mean = mean(BG_DA);
        BG_AA = BG_Aem(2:2:end); BG_AA_mean = mean(BG_AA);
        BG_DD = BG_Dem(1:2:end); BG_DD_mean = mean(BG_DD);
        
    else
        DA_data = Aem_data(2:2:end);
        AA_data = Aem_data(3:2:end);
        DD_data = Dem_data(2:2:end);
        
        BG_DA = BG_Aem(2:2:end,2); BG_DA_mean = mean(BG_DA);
        BG_AA = BG_Aem(3:2:end,2); BG_AA_mean = mean(BG_AA);
        BG_DD = BG_Dem(2:2:end,2); BG_DD_mean = mean(BG_DD);
            
    end
    
    subplot(3,1,1,'Fontsize',20)
    plot(times,DA_data,'r-', 'LineWidth',2);
    hold on
    plot(times,DD_data,'g-', 'LineWidth',2);
    hold on
    plot(times,AA_data,'k-', 'LineWidth',2); xlim('auto'); title('Raw intensities','Fontsize',20); xlabel('','Fontsize',20);ylabel('Intensity','Fontsize',20);


    %bkgd subtract
    DA_data_BS = DA_data - BG_DA_mean;
    DD_data_BS = DD_data - BG_DD_mean;
    AA_data_BS = AA_data - BG_AA_mean;

    subplot(3,1,2,'Fontsize',20)
    plot(times,DA_data_BS,'r-', 'LineWidth',2);
    hold on
    plot(times,DD_data_BS,'g-', 'LineWidth',2);
    hold on
    plot(times,AA_data_BS,'k-', 'LineWidth',2); xlim('auto'); title('Raw intensities','Fontsize',20); xlabel('','Fontsize',20);ylabel('Intensity','Fontsize',20);


    %calculate FRET & Stoichiometry. Careful with no. of frames as
    %can end on Dex *or* Aex
    E = zeros(length(AA_data_BS),1);
    S = zeros(length(AA_data_BS),1);
    
    for i = 1:length(AA_data_BS);
        E(i) = DA_data_BS(i)./(DA_data_BS(i) + DD_data_BS(i));
        S(i) = (DA_data_BS(i)+DD_data_BS(i))./(DA_data_BS(i)+DD_data_BS(i)+AA_data_BS(i));
    end
    
    subplot(3,1,3,'Fontsize',20)
    plot(times,E,'r-', 'LineWidth',2);ylim([0 1]);xlim('auto'); xlabel('Time (s)','Fontsize',20);
    hold on
    plot(times,S,'k-', 'LineWidth',2);ylim([0 1]);xlim('auto'); title('Sraw (black) / Eraw (red)','Fontsize',20); xlabel('Frame','Fontsize',20);
end
