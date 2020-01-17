%%
clear all
restoredefaultpath

addpath('~/Documents/MATLAB/fieldtrip/')
ft_defaults

addpath('~/Documents/MATLAB/meeg_hyperalignment/')

load('/Users/b1019548/Desktop/MessingAround/Markov_Orig_ERFs/AllERFs.mat')

%%
cfg=[];
cfg.lpfilter='yes';
cfg.lpfreq=25;
    
for ii=1:length(SNDmat)
    SNDmat{ii}=ft_preprocessing(cfg, SNDmat{ii});
end

%%
cfg=[];
cfg.scale=false;

for ii=1:length(SNDmat)
    tmp=nw_procrustes_calctr(cfg, SNDmat{1}, SNDmat{ii});
    SNDmat_p{ii}=tmp.dataP;
    tr{ii}=tmp.tr;
    clear tmp
end


%%
GA_SND_p = ft_timelockgrandaverage([], SNDmat_p{:});
GA_SND = ft_timelockgrandaverage([], SNDmat{:});

%%
cfg=[];
cfg.layout='neuromag306mag.lay';
ft_multiplotER(cfg, GA_SND ,GA_SND_p)

%%

OMmat_RD_p = OMmat_RD;
OMmat_MM_p = OMmat_MM;
OMmat_MP_p = OMmat_MP;
OMmat_OR_p = OMmat_OR;

for ii=1:length(OMmat_RD)
    cfg=[];
    cfg.tr=tr{ii};
    OMmat_RD_p{ii}=nw_procrustes_applytr(cfg, OMmat_RD{ii});
    OMmat_MM_p{ii}=nw_procrustes_applytr(cfg, OMmat_MM{ii});
    OMmat_MP_p{ii}=nw_procrustes_applytr(cfg, OMmat_MP{ii});
    OMmat_OR_p{ii}=nw_procrustes_applytr(cfg, OMmat_OR{ii});
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

%%
cfg=[];
cfg.xlim=[-.1 .3];
cfg.layout='neuromag306mag.lay';
ft_multiplotER(cfg, GA_OMrd ,GA_OMor)

%% LET'S TRY INVERSE TRANSFORM
% NOW TRANSFORM IN SPACE OF SUBJ 1 ... LET'S MOVE TRANSFORMED DATA OF
% SUBJECT 2 BACK TO ORIGINAL SPACE

cfg=[];
cfg.tr=tr{2}; %how to get data from subject 2 into subject 1
cfg.inverse=true;
SNDsub2_inv=nw_procrustes_applytr(cfg, SNDmat_p{2});

%%
cfg=[];
cfg.xlim=[-.1 .3];
cfg.layout='neuromag306mag.lay';
ft_multiplotER(cfg, SNDsub2_inv)
figure;ft_multiplotER(cfg, SNDmat{2})

%%
cfg=[];
cfg.xlim=[-.1 .3];
cfg.layout='neuromag306mag.lay';
ft_multiplotER(cfg, GA_OMrd_p, GA_OMmm_p, GA_OMmp_p ,GA_OMor_p)
figure; ft_multiplotER(cfg, GA_OMrd, GA_OMmm, GA_OMmp,GA_OMor)


%%

cfg = [];
cfg.method      = 'template';                         % try 'distance' as well
cfg.template    = 'neuromag306mag_neighb.mat';                % specify type of template
cfg.layout      = 'neuromag306mag_helmet.mat';                % specify layout of channels
cfg.feedback    = 'yes';                              % show a neighbour plot
neighbours      = ft_prepare_neighbours(cfg, GA_OMmm_p); %


%%
cfg=[];
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.correctm    = 'cluster';
cfg.correcttail = 'prob';
cfg.latency=[0 .3];
cfg.numrandomization = 1000;
cfg.neighbours  = neighbours;
cfg.spmversion = 'spm12';
Nsub = length(SNDmat);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2;

stat=ft_timelockstatistics(cfg, OMmat_RD{:}, OMmat_OR{:});
statP=ft_timelockstatistics(cfg, OMmat_RD_p{:}, OMmat_OR_p{:});

%%
cfg=[];
cfg.parameter='stat';
cfg.xlim=[0 .3];
cfg.layout='neuromag306mag_helmet.lay';
cfg.zlim='absmax';
ft_multiplotER(cfg, stat, statP) % tatP.posclusters(1) --> prob: 0.0020

%%

statP.stat2=statP.stat.*(statP.posclusterslabelmat==1);
stat.stat2=stat.stat.*(statP.posclusterslabelmat==1);

cfg=[];
cfg.parameter='stat2';
cfg.xlim=[0 .3];
cfg.layout='neuromag306mag_helmet.lay';
cfg.zlim='absmax';
ft_multiplotER(cfg, statP, stat) 

%% TEST 1 Subj single trials

load('data.mat')

%% data from MRGU --> vp 27

cfg=[];
cfg.tr=tr{27}; 
data_P=nw_procrustes_applytr(cfg, data);

%%
test_P=ft_timelockanalysis([], data_P); 
test=ft_timelockanalysis([], data);

cfg=[];
cfg.layout='neuromag306mag_helmet.lay';
cfg.zlim='absmax';
ft_multiplotER(cfg, test, test_P) 

%% data from MRGU --> backwards

cfg=[];
cfg.tr=tr{27}; 
cfg.inverse=true;
data_P_back=nw_procrustes_applytr(cfg, data_P);

%%

test_P_back=ft_timelockanalysis([], data_P_back); 

%%

cfg=[];
cfg.layout='neuromag306mag_helmet.lay';
cfg.zlim='absmax';
ft_multiplotER(cfg, test) 
figure; ft_multiplotER(cfg, test_P_back) 


