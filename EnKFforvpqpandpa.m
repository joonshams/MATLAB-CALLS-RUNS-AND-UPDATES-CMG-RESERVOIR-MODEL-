clear all
close all
clc
filename='sandtank2D.out';
td=364;
sgx=40;
sgz=40;
%% Load True Model
load('pres_truth_enkf2020');
load('sg_truth_enkf2020');
sg_truth(sg_truth<=0)=0.01;
pres_truth=pres_truth*0.00689476; %psia to mPA

%% plotting the time-lapse results of the truth
% for i=1:1:td
%     subplot(1,2,1)
%    imagesc(sg_truth(:,:,i));
%    colorbar
%    xlabel('distance, x (ft)')
%    ylabel('depth, z (ft)')
%    title('True Gas Saturation');
%    
%    subplot(1,2,2)
%    imagesc(pres_truth(:,:,i));
%    colorbar
%    xlabel('distance, x (ft)')
%    ylabel('depth, z (ft)')
%    title('True Pressure (Pa)');
%    pause(0.1);
% end
%store the pressure monitoring data at the inj/prod site
% x_inj=36;
% z_inj=38;
% x_prod=4;
% z_prod=38;
% pres_pm=pres_truth(x_prod,z_prod,:)./0.00689476; %psia
% pres_im=pres_truth(x_inj,z_inj,:)./0.00689476; %psia

%% Generate Observations of Vp Qp 
%things to change for new models: grid perm, phi, b, and pressure
for i=1
Ks=34; % GPa Grain Bulk Modulus
Gs=35.3; %GPa Grain Shear Modulus
rhos=2585; %kg/m^3 Grain density

%Matrix (Kdry and Km, Gm and Gdry are the same)
Km=8.67;Kdry=Km; %Gpa Matrix Bulk Modulus 
Gm=6.61; Gdry=Gm;%Gpa Matrix Shear Modulus
phi(1:5,1:40)=0.1;%mD CHANGE THESE FOR NEW MODELS
phi(6:40,1:40)=0.3;%mD
% Kdry = Ks*(1 - phi)^(4/(1-phi)); %Carcione and Picotti 2006f
% Gdry = Kdry*Gs/Ks;
%T=2.5; %Tortuosity
%perm(1:40,1:40)=0.01;
perm(1:5,1:40)=0.001;%mD CHANGE THESE FOR NEW MODELS
perm(6:40,1:40)=0.1;%mD
perm=perm*1e-3;%D
perm=perm*9.869233e-13;%D to m2

%Brine
%Kw=2.4;%GPa Brine Bulk Modulus
%rhow=1040; %kg/m^3 brine density
nuw=0.00018; %Pa.s brine visc

%Gas @ 1-km depth
%Kg=0.01; %Gpa
%rhog=100; %kg/m^3
nug=0.00002;%Pa.s visc of gas

%patchwmodified for loop for varying kg and f values to reproduce Figure 2
freq=[50000 100000 250000 500000]; %Hz
%dx=0.000886; %grid spacing in meters
b=0.077*0.3048; % m (r=sqrt((z*x)/pi))

%flprop properties
method=2;
sal=0; %NaCl ppm
og=25; % oil gravity
gg=0.8; %Specific Gravity of gas
gor=0; %L/L
giib=1; %gas index in brine
giio=0; %gas index in oil
%P=PRES./1000000; %pore pressure (MPa)
T=25; %degC Rock temp
So=0; %saturation of oil
% Sg=0.8; %saturation of gas

end%%All the petrophysical properties used for patchw and flprop

%% OBSERVATIONS FWD WHITE MODEL 

for k=1:1:td
    for j=1:1:sgz
        for i=1:1:sgx
            a(i,j,k)=sg_truth(i,j,k)^(1/3)*b;
            [~,~,~,~,rhow(i,j,k),Kb(i,j,k),~,~,~,...
                ~,rhog(i,j,k),Kg(i,j,k),~]=flprop(method,sal,og,gg,gor,giib,giio,pres_truth(i,j,k),T,So,sg_truth(i,j,k));
            fl=[Kg(i,j,k).*1e9 Kb(i,j,k).*1e9; rhog(i,j,k) rhow(i,j,k); nug nuw];
            [vp(i,j,k),~,atn(i,j,k),~,~,~]= ...
                patchw(Kdry*1e9,Gdry*1e9,Ks*1e9,Gs*1e9,rhos,phi(i,j),perm(i,j),fl,sg_truth(i,j,k),a(i,j,k) ...
                ,freq(3));
            pres_o(i,j,k)=pres_truth(i,j,k)/0.00689476+normrnd(0,2);
            while(pres_o(i,j,k)<=5)
                pres_o(i,j,k)=pres_truth(i,j,k)/0.00689476+normrnd(0,2);
            end
             while(pres_o(i,j,k)>=40)
                pres_o(i,j,k)=pres_truth(i,j,k)/0.00689476+normrnd(0,2);
            end
            while(atn(i,j,k)<0.01)
                atn(i,j,k)=0.01+abs(normrnd(0,0.01));
                
            end
            %patchw(kdry,mudry,k0,mu0,ro0,phi,perm,fl,sg,a,f)
        end
    end
    atn(isnan(atn))=0.01+abs(normrnd(0,0.01));
    
    vp(isnan(vp))=3600+normrnd(0,20);
end
Xo(:,1,:)=reshape(vp,sgx*sgz,1,td);
Xo(:,2,:)=reshape(atn,sgx*sgz,1,td);
pobs=zeros(40,40,td);
pobs(38,36,:)=1;
pobs(38,4,:)=1;
pobs=reshape(pobs,sgx*sgz,1,td);
Xo(:,3,:)=reshape(pres_o,sgx*sgz,1,td).*pobs;

clear vp atn

%% 
%Xa:analysis, Xb: background, Xf:forecast, and Xo: observations
ens=20;
for i=1:ens
Xa(1,:,1,i)=reshape(sg_truth(:,:,1),1,40*40)+abs(normrnd(0,0.05,[1,sgx*sgz]));%gas saturation
Xa(2,:,1,i)=reshape(pres_truth(:,:,1)/0.00689476,1,1600)+abs(normrnd(0,max(max(pres_truth(:,:,1)))*0.1,[1,sgx*sgz]));%pressure MPa
Xa(3,:,1,i)=reshape(phi(:,:),1,40*40);
Xa(4,:,1,i)=reshape(perm,1,sgx*sgz);
Xa(5,:,1,i)=Xo(:,1,1)';
Xa(6,:,1,i)=Xo(:,2,1)';
end
perm_idx=reshape(perm,1,sgx*sgz);
phi_idx=reshape(phi,1,sgx*sgz);
%Xf(:,:,2)=Xa;
for k=2:td
    for en=1:ens
        %% Use CMG to generate forecast
        SGINT_cmg=reshape(Xa(1,:,k-1,en),sgx,sgz,1)';
        SGINT_cmg(SGINT_cmg<=0)=0.05+abs(normrnd(0,0.02));
        SGINT_cmg(SGINT_cmg>=1)=1-abs(normrnd(0,0.1));
        pres_cmg=reshape(Xa(2,:,k-1,en),sgx,sgz,1)';%pres_truth(:,:,k)/0.00689476%change to psia for cmg
        pres_cmg(pres_cmg<=3)=5+abs(normrnd(0,3));
        pres_cmg(pres_cmg>=40)=35+abs(normrnd(0,3));
        [sg_f(:,:,k,en),pres_f(:,:,k,en)]=SWINTcmg(1-SGINT_cmg,pres_cmg,k,en);
        Xf(1,:,k,en)=reshape(sg_f(:,:,k,en),1,sgx*sgz)+normrnd(0,0.05,[1,sgx*sgz]);
        Xf(2,:,k,en)=reshape(pres_f(:,:,k,en),1,sgx*sgz)+normrnd(0,1,[1,sgx*sgz]);
        Xf(3,:,k,en)=Xa(3,:,k-1,en);
        Xf(4,:,k,en)=Xa(4,:,k-1,en);
        Xf(5,:,k,en)=Xa(5,:,k-1,en);
        Xf(6,:,k,en)=Xa(6,:,k-1,en);
              
            %% H-operator to generate Vp and Qp from Sg
                  
        for idx=1:sgx*sgz
            while (Xf(1,idx,k,en)>=1)
                Xf(1,idx,k,en)=normrnd(0.95,0.05);
            end
            while (Xf(1,idx,k,en)<=0)
                Xf(1,idx,k,en)=normrnd(0.05,0.05);
            end
%             while (Xf(2,idx,k,en)>=30)
%                 Xf(2,idx,k,en)=normrnd(29,3);
%             end
            while (Xf(2,idx,k,en)<=0)
                Xf(2,idx,k,en)=normrnd(10,3);
            end
            a=Xf(1,idx,k,en)^(1/3)*b;
            [~,~,~,~,rhow,Kb,~,~,~,...
                ~,rhog,Kg,~]=flprop(method,sal,og,gg,gor,giib,giio,Xf(2,idx,k,en)*0.00689476,T,So,Xf(1,idx,k,en));
            fl=[Kg.*1e9 Kb.*1e9; rhog rhow; nug nuw];
            [vp,~,atn,~,~,~]= ...
                patchw(Kdry*1e9,Gdry*1e9,Ks*1e9,Gs*1e9,rhos,phi_idx(1,idx),perm_idx(1,idx),fl,Xf(1,idx,k,en),a ...
                ,freq(4));
            atn(isnan(atn))=0.01+abs(normrnd(0,0.1));
            atn(atn<0.01)=0.01+abs(normrnd(0,0.01));
            vp(isnan(vp))=3600+normrnd(0,100);
            HXf(1,idx,k,en)= Xf(1,idx,k,en);
            HXf(2,idx,k,en)= Xf(2,idx,k,en);
            HXf(3,idx,k,en)= Xf(3,idx,k,en);
            HXf(4,idx,k,en)= Xf(4,idx,k,en);
            HXf(5,idx,k,en)=vp;
            HXf(6,idx,k,en)=atn; 
            Xf(5,idx,k,en)=vp;
            Xf(6,idx,k,en)=atn;
        end
    end
    
    Pm=[pres_truth(4,38,k);pres_truth(36,38,k)];
%     [Xf(1,:,k,:),mu_xf1(:,k),stdev_xf1(:,k)]=logit(Xf(1,:,k,:));
%     [Xf(2,:,k,:),mu_xf2(:,k),stdev_xf2(:,k)]=logtrans(Xf(2,:,k,:));
%     
%     [Xo(:,1,k),mu_xo1(:,k),stdev_xo1(:,k)]=logtrans(Xo(:,1,k));
%     [Xo(:,2,k),mu_xo2(:,k),stdev_xo2(:,k)]=logtrans(Xo(:,2,k));
%     
%     [HXf(1,:,k,:),mu_Hxf1(:,k),stdev_Hxf1(:,k)]=logtrans(HXf(1,:,k,:));
%     [HXf(2,:,k,:),mu_Hxf2(:,k),stdev_Hxf2(:,k)]=logtrans(HXf(2,:,k,:));
    
    
    [Xa2]=modelspaceENKF(Xf(:,:,k,:),Xo(:,:,k),HXf(:,:,k,:),Pm,sgx,sgz,ens);
    %Xa2 from nonlinearENKf has dimensions 6X1600
    Xa(:,:,k,:)=permute(Xa2,[1,2,4,3]);
    Xfm2=permute(Xf(:,:,k,:),[2,4,1,3]);
    
    for n=1:6%change this if number of state vector ~= 6
        Xfmean(:,n,k)=mean(Xfm2(:,:,n)')';
    end
    

    
    Xam2=permute(Xa(:,:,k,:),[2,4,1,3]);
    for n=1:6%change this if number of state vector ~= 6
        Xamean(:,n,k)=sum(Xam2(:,:,n)')';
    end
    
    plotsg=reshape(Xa(1,:,k,:),sgx,sgz,ens);
    savesg2=mean(permute(plotsg(:,:,:),[3,1,2]));
    savesg(:,:,k)=permute(savesg2,[2,3,1]);
    
    plotpa=reshape(Xa(2,:,k,:),sgx,sgz,ens);
    savepa2=mean(permute(plotpa(:,:,:),[3,1,2]));
    savepa(:,:,k)=permute(savepa2,[2,3,1]);
    
    
    plotk=reshape(Xa(4,:,k,1),sgx,sgz);
    savek(:,:,k)=plotk;
    plotphi=reshape(Xa(3,:,k,1),sgx,sgz);
    savephi(:,:,k)=plotphi;
    plotvp=reshape(Xa(5,:,k,1),sgx,sgz);
    savevp(:,:,k)=plotvp;
    plotqp=reshape(Xa(6,:,k,1),sgx,sgz);
    saveqp(:,:,k)=plotqp; 
    
    plotXfsg=reshape(Xfmean(:,1,k),sgx,sgz);
    saveXfsg(:,:,k)=plotXfsg;
    plotXfpa=reshape(Xfmean(:,2,k),sgx,sgz);
    saveXfpa(:,:,k)=plotXfpa;
    
    
    figure (1)

    subplot(3,3,1)
    imagesc(savesg(:,:,k))
    colorbar
    caxis([0 1])
    title(['Analysis: SG at TIME= ',num2str(k),' days']);
    
    subplot(3,3,4)
    imagesc(savepa(:,:,k))
    colorbar
    caxis([0 30])
    title(['Analysis: Pressure at TIME= ',num2str(k),' days']);
    
    subplot(3,3,7)
    imagesc(log(1./plotqp))
    colorbar
    title(['Analysis: log(Q_p) at TIME= ',num2str(k),' days']);
    
    subplot(3,3,8)
    imagesc(reshape(Xo(:,3,k),40,40,1))
    colorbar
    title(['Observed Pressure at TIME= ',num2str(k),' days']);
    
    subplot(3,3,9)
    imagesc(plotvp)
    colorbar
    title(['Analysis: Vp at TIME= ',num2str(k),' days']);
    
    subplot(3,3,2)
    imagesc(plotXfsg)
    colorbar  
    caxis([0 1])
    title(['Forecast: SG at TIME= ',num2str(k),' days']);
    
    subplot(3,3,5)
    imagesc(plotXfpa)
    colorbar 
    caxis([0 30])
    title(['Forecast: Pressure at TIME= ',num2str(k),' days']);
    
    subplot(3,3,3)
    imagesc(abs(savesg(:,:,k)-sg_truth(:,:,k)))
    colorbar 
    caxis([0 1])
    title(['Error: SG at TIME= ',num2str(k),' days']);
    
    subplot(3,3,6)
    imagesc(abs(savepa(:,:,k)-pres_truth(:,:,k)./0.00689476))
    colorbar 
    caxis([0 30])
    title(['Error: Pressure at TIME= ',num2str(k),' days']);
end
%% Assimilation Error
sgdiff2=((savesg-sg_truth).^2).^(1/2);
padiff2=((savepa-pres_truth./0.00689476).^2).^(1/2);
sgdiff_mean2=mean(mean(sgdiff2));
padiff_mean2=mean(mean(padiff2));
for i=1:td
    sgdiffmean2(1,i)=sgdiff_mean2(:,:,i);
    padiffmean2(1,i)=padiff_mean2(:,:,i);
end
%% Forecast Error
for i=1:ens
sgfe=sum(sg_f(:,:,:,1:20))
end
%% Plot Error Statistics 
figure (2)
hold on
subplot(2,1,1)
plot(1:1:td,sgdiffmean2)
xlabel('time (days)');
ylabel('RMSE error');
title('Assimilated Gas Saturation RMSE ');

subplot(2,1,2)
hold on
plot(1:1:td,padiffmean2)
xlabel('time (days)');
ylabel('RMSE error');
title('Assimilated Pressure RMSE  (Mpa)');

figure (3)
hold on
subplot(2,1,1)
plot(1:1:td,permute(sgdiff2(15,15,1:td),[3,2,1]))
xlabel('time (days)');
ylabel('RMSE error');
title('Assimilated Gas Saturation RMSE ');

subplot(2,1,2)
hold on
plot(1:1:td,permute(padiff2(15,15,1:td),[3,2,1]))
xlabel('time (days)');
ylabel('RMSE error');
title('Assimilated Pressure RMSE  (Mpa)');

