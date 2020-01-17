
%%
clear all
restoredefaultpath

addpath('~/Documents/MATLAB/fieldtrip/')
ft_defaults

%%
allfiles=dir('/Users/b1019548/.CMVolumes/Obob/obob/staff/gdemarchi/data/markov/ERFs/*.mat');

%%
for ii=1:length(allfiles)
    clear avg*
    load(fullfile('/Users/b1019548/.CMVolumes/Obob/obob/staff/gdemarchi/data/markov/ERFs/', allfiles(ii).name))
    
    cfg=[];
    cfg.channel = 'MEGMAG';
    SNDmat{ii}=ft_selectdata(cfg, avgSNDALL_RD);
    OMmat_RD{ii}=ft_selectdata(cfg, avgOMALL_RD);
    OMmat_MM{ii}=ft_selectdata(cfg, avgOMALL_MM);
    OMmat_MP{ii}=ft_selectdata(cfg, avgOMALL_MP);
    OMmat_OR{ii}=ft_selectdata(cfg, avgOMALL_OR);
end

clear avg*

save('/Users/b1019548/Desktop/MessingAround/Markov_Orig_ERFs/AllERFs.mat', 'SNDmat', 'OMmat*')

%% Normal Grand Average

GA_SND = ft_timelockgrandaverage([], SNDmat{:});
GA_OMrd = ft_timelockgrandaverage([], OMmat_RD{:});

%%
cfg=[];
cfg.layout='neuromag306mag.lay';
ft_multiplotER(cfg, GA_OMrd)

%%
X=SNDmat{1}.avg;
SNDmat_p = SNDmat;

for ii=1:length(SNDmat)
    Y{ii}=SNDmat{ii}.avg;
    [d,SNDmat_p{ii}.avg,tr{ii}] = procrustes(X,Y{ii});
end

%%

GA_SND_p = ft_timelockgrandaverage([], SNDmat_p{:});

%%
cfg=[];
cfg.layout='neuromag306mag.lay';
ft_multiplotER(cfg, GA_SND ,GA_SND_p)

%% transform omissions
%Z = b*Y*T + c;

OMmat_RD_p = OMmat_RD;
OMmat_MM_p = OMmat_MM;
OMmat_MP_p = OMmat_MP;
OMmat_OR_p = OMmat_OR;

for ii=1:length(OMmat_RD)
    OMmat_RD_p{ii}.avg = tr{ii}.b * OMmat_RD{ii}.avg * tr{ii}.T + tr{ii}.c;
    OMmat_MM_p{ii}.avg = tr{ii}.b * OMmat_MM{ii}.avg * tr{ii}.T + tr{ii}.c;
    OMmat_MP_p{ii}.avg = tr{ii}.b * OMmat_MP{ii}.avg * tr{ii}.T + tr{ii}.c;
    OMmat_OR_p{ii}.avg = tr{ii}.b * OMmat_OR{ii}.avg * tr{ii}.T + tr{ii}.c;
end

%%

GA_OMrd_p = ft_timelockgrandaverage([], OMmat_RD_p{:});
GA_OMmm_p = ft_timelockgrandaverage([], OMmat_MM_p{:});
GA_OMmp_p = ft_timelockgrandaverage([], OMmat_MP_p{:});
GA_OMor_p = ft_timelockgrandaverage([], OMmat_OR_p{:});

GA_OMrd = ft_timelockgrandaverage([], OMmat_RD{:});
GA_OMmm = ft_timelockgrandaverage([], OMmat_MM{:});
GA_OMmp = ft_timelockgrandaverage([], OMmat_MP{:});
GA_OMor = ft_timelockgrandaverage([], OMmat_OR{:});

%%
cfg=[];
cfg.xlim=[-.1 .3];
cfg.layout='neuromag306mag.lay';
ft_multiplotER(cfg, GA_OMrd ,GA_OMrd_p)


%%
cfg=[];
cfg.xlim=[-.1 .3];
cfg.layout='neuromag306mag.lay';
ft_multiplotER(cfg, GA_OMrd_p ,GA_OMor_p)


