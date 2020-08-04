function [Info_tot, Info_thr_tot] = compute_SpectralInformation(opts, varargin)

[~,ncond]= size(varargin);
S = [];
for ii = 1 : ncond
    [nbchan, nfreq,m] = size(varargin{ii});
    S = [S;ii*ones(m,1)];
end
    
Info_tot = zeros(nbchan,nfreq);
Info_thr_tot = zeros(nbchan,nfreq);    


for chi = 1: nbchan 
    for kk = 1: nfreq
        R_raw = [];
        for ii = 1 : ncond
            R_raw = [R_raw; squeeze(varargin{ii}(chi,kk,:))];
        end       
        [R,nT] = buildr(S,R_raw);
        opts.nt = nT;
        R_bin = binr(R,nT,ncond,'eqpop');
        I = information(R_bin,opts,'I');
        I(I<0) = 0;
        Info_tot(chi,kk) = I(1);
        Info_thr_tot(chi,kk) = prctile(I,95);
    end  
end
