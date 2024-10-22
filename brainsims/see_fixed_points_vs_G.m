clear; clc 
files = dir('NFxps_vs_Gs_N*'); 
figure(1); clf; hold on 
Nf = []; 
alf = .5; 
colores = lines(numel(files)); 
for n = 1:numel(files)
    load(files(n).name)
    plot(Gs,Nfxps,'LineWidth',2,'color',[colores(n,:), alf])
    scatter(Gs,Nfxps,60,colores(n,:),'filled')
    alpha(alf)
    Nf = [Nf; Nfxps]; 
end 
N = max(Nf); 
plot(Gs,N,'k-o')
set(gca,'YScale','linear')
ylim([0 40])
grid on 
%%
clc 
disp([Gs(Gs==0.08), N(Gs==0.08)]) 
disp([Gs(Gs==0.55), N(Gs==0.55)]) 