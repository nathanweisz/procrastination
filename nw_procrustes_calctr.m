function dat_trans=nw_procrustes_calctr(cfg, targetavg, sourceavg)
%Convenience function to obtain linear transformation data to align MEG
%data from a "source"-participant to a "target"-participant (latter could
%also be a grand average).
%
%Rationale the transformation information is done using averages or mTRFs
%contained in a *.avg field. The transformation information can be applied in
%a separate function also to single trials.
%
%Input:
%       -cfg.scale = per default false (i.e. scaling component = 1)
%       -sourceavg = timelock-like data i.e. must contain avg-field
% 
%Output: dat_trans structure with following fields
%       - dataP = The Procustes transformed data structure
%       - targetgrad = the grad structure of the target data
%       - d = dissimilarity metric
%       - tr = actual transformation info (translation, rotation, scaling)
%
%See also nw_procrustes_applytr.m
%
%Jan 2020: First Implementation NW

cfg.scale = ft_getopt(cfg, 'scale', false, 1);


if istrue(cfg.scale)
    [d,tmp,tr] = procrustes(targetavg.avg', sourceavg.avg', 'scaling', true);    
elseif ~istrue(cfg.scale)
    [d,tmp,tr] = procrustes(targetavg.avg', sourceavg.avg', 'scaling', false);
end

sourceavg=rmfield(sourceavg, 'grad');

dat_trans=[];
dat_trans.dataP=sourceavg;
dat_trans.dataP.avg=tmp';

if isfield(targetavg, 'grad')
    dat_trans.targetgrad=targetavg.grad;
end

dat_trans.d=d;
dat_trans.tr=tr;
dat_trans.tr.c=mean(tr.c);



