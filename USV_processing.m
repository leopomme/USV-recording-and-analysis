function varargout = usvseg08r7(action)

%requires Matlab the signal processing toolbox

% global variable
global gh;

% Initiation ------------------------------------------------------
if nargin==0
    gh.version = 'USV_processing';
    gh.fftsize = 512; % FFT size = window size
    gh.pth = cd;
    % load parameters
    prmfile = [fileparts(mfilename('fullpath')) filesep 'usvseg_prm.mat'];
    if exist(prmfile)==2, l=load(prmfile); gh.prm = l.prm;
    else % default values
        prm.timestep = 0.0005; prm.freqmin = 30000; prm.freqmax = 120000;
        prm.threshval = 4.5;   prm.durmin = 0.005;  prm.durmax = 0.3;
        prm.gapmin = 0.030;    prm.margin = 0.015;  prm.wavfileoutput = 1;
        prm.imageoutput = 1;   prm.imagetype = 0;   prm.traceoutput = 0;
        prm.readsize = 10;     prm.mapL = 1000;     prm.mapH = 6000;
        gh.prm = prm;
    end    % make GUI
    makegui;
    cla(gh.haxes(1)); cla(gh.haxes(2)); cla(gh.haxes(3));
    return
end
% open ----------------------------------------------------------------
if strcmp(action,'open')
    % clear axes
    cla(gh.haxes(1)); cla(gh.haxes(2)); cla(gh.haxes(3));
    % unable
    set([gh.hpush_thrs gh.hpush_dtct gh.hpush_save gh.hpush_sgmt gh.hpush_play gh.hpush_swav],'Enable','Off');
    % fetch parameter
    fetchparams;
    timestep = gh.prm.timestep;
    readsize = gh.prm.readsize;
    fftsize = gh.fftsize;
    % get file
    [fn,pth] = uigetfile([gh.pth filesep '*.wav']);
    % show file name
    set(gh.htext_fnam,'string',['File: ' pth fn]);
    busytoggle(1);
    % read wav
    ai = audioinfo([pth fn]);
    fs = ai.SampleRate;
    wav = audioread([pth fn],[1 min(readsize*fs,ai.TotalSamples)]);
    % show fs
    set(gh.htext_fsvl,'string',num2str(fs/1000));
    % multitaper spectrogram
    step = round(timestep*fs);
    mtsp = multitaperspec(wav,fs,fftsize,timestep);
    mtsimrng = median(mtsp(:)) + [0 40];
    fvec = [0; (1:(fftsize/2))'/fftsize*fs];
    tvec = ((0:(size(mtsp,2)-1))*step+fftsize/2)/fs;
    % default view range : 2 seconds
    xl = [0 2];
    yl = [0 max(fvec)/1000];
    % save in global variable
    gh.fn = fn;
    gh.pth = pth;
    gh.fs = fs;
    gh.wav = wav;
    gh.step = step;
    gh.mtsp = mtsp;
    gh.mtsimrng = mtsimrng;
    gh.fs = fs;
    gh.fvec = fvec;
    gh.tvec = tvec;
    gh.xl = xl;
    gh.yl = yl;
    % draw mtsp
    drawmtsp;
    % enable
    set(gh.hpush_flat,'Enable','On');
    busytoggle(0);
    return;
end
% flattening ------------------------------------------------------------
if strcmp(action,'flat')
    busytoggle(1);
    % clear axes
    cla(gh.haxes(2)); cla(gh.haxes(3));
    % flattening
    fltnd = flattening(gh.mtsp);
    fltndimrng = [0 30];
    % save
    gh.fltnd = fltnd;
    gh.fltndimrng = fltndimrng;
    % draw fltnd
    drawfltnd;
    % enable
    set(gh.hpush_thrs,'Enable','On');
    busytoggle(0);
    return;
end
% thresholding ------------------------------------------------------------
if strcmp(action,'thrs')
    % clear axes
    cla(gh.haxes(3));
    % fetch
    fetchparams;
    freqmin = gh.prm.freqmin; freqmax = gh.prm.freqmax; threshval = gh.prm.threshval; 
    fs = gh.fs; fvec = gh.fvec; tvec = gh.tvec; fltnd = gh.fltnd;
    busytoggle(1);    
    % threshold calculation with threshval*sigma (SD) of background noise 
    thresh = estimatethresh(fltnd,fs,freqmin,freqmax,threshval);
    % thresholding
    thrshd = thresholding(fltnd,fs,freqmin,freqmax,thresh);
    % draw
    set(0,'CurrentFigure',gh.hfig);
    set(gh.hfig,'CurrentAxes',gh.haxes(3));
    imagesc(tvec,fvec/1000,thrshd); axis xy
    colormap(gca,gray);
    set(gca,'tickdir','out'); box off;
    ylabel('Frequency (kHz)');
    set(gh.haxes,'xlim',gh.xl,'ylim',gh.yl);
    drawnow;
    % enable
    set(gh.hpush_dtct,'Enable','On');
    % save
    gh.thrshd = thrshd; gh.thresh = thresh;
    busytoggle(0);
    return;
end
% detection ------------------------------------------------------------
if strcmp(action,'dtct')
    busytoggle(1);
    % clear axes
    cla(gh.haxes(3));   
    % fetch
    fetchparams;
    fs = gh.fs; timestep = gh.prm.timestep;
    freqmin = gh.prm.freqmin; freqmax = gh.prm.freqmax;
    durmin = gh.prm.durmin; durmax = gh.prm.durmax;
    gapmin = gh.prm.gapmin; margin = gh.prm.margin;
    mtsp = gh.mtsp; fltnd = gh.fltnd; thrshd = gh.thrshd;
    % onset/offset detection
    [onoffset,onoffsetm] = detectonoffset(thrshd,fs,timestep,gapmin,durmin,durmax,margin);
    % peak tracking
    [freqtrace,amptrace,maxampval,maxampidx,maxfreq] = specpeaktracking(mtsp,fltnd,fs,timestep,freqmin,freqmax,onoffset,margin);
    % save
    gh.onoffset = onoffset; gh.onoffsetm = onoffsetm;
    gh.freqtrace = freqtrace; gh.amptrace = amptrace;
    gh.maxampval = maxampval; gh.maxfreq = maxfreq; gh.maxampidx = maxampidx;
    % draw trace
    drawtrace;
    % enable
    set(gh.hpush_save,'Enable','On');
    set(gh.hpush_sgmt,'Enable','On');
    set(gh.hpush_play,'Enable','On');
    busytoggle(0);
    return;
end
% save ----------------------------------------------------------------
if strcmp(action,'save')
    busytoggle(0);
    % fetch
    onoffset = gh.onoffset;
    maxampval = gh.maxampval; maxfreq = gh.maxfreq;
    dur = diff(onoffset,[],2);
    [~,prefix,~] = fileparts(gh.fn);
    % save CSV
    [fn,pth] = uiputfile([gh.pth filesep prefix '_dat.csv'],'file name');
    fid = fopen([pth fn],'wt');
    fprintf(fid,'#,start,end,duration,maxfreq,maxamp\n');
    if ~isempty(onoffset)
        for n=1:size(onoffset,1)
            fprintf(fid,'%d,%f,%f,%f,%f,%f\n',n,onoffset(n,1),onoffset(n,2),dur(n),maxfreq(n)/1000,maxampval(n));
        end
    end
    fclose(fid);
    busytoggle(1);
    return;
end
% segment -------------------------------------------------------------
if strcmp(action,'sgmt')
    busytoggle(0);
    % fetch
    fetchparams;
    margin = gh.prm.margin; timestep = gh.prm.timestep;
    wavflg = gh.prm.wavfileoutput; imgflg = gh.prm.imageoutput; 
    fltflg = gh.prm.imagetype; trcflg = gh.prm.traceoutput;
    wav = gh.wav; fs = gh.fs; mtsp = gh.mtsp; fltnd = gh.fltnd;
    tvec = gh.tvec; onoffset = gh.onoffset;
    freqtrace = gh.freqtrace;amptrace = gh.amptrace;
    fltndimrng = gh.fltndimrng; mtsimrng = gh.mtsimrng;
    % file name
    [~,prefix,~] = fileparts(gh.fn);
    outp = uigetdir(gh.pth);
    % output
    startid = 1;
    if fltflg==1, inputimg = fltnd; imrng = fltndimrng;
    else,         inputimg = mtsp;   imrng = mtsimrng; end
    segfun(startid,outp,prefix,inputimg,imrng,wav,fs,timestep,margin,onoffset,tvec,amptrace,freqtrace,wavflg,imgflg,trcflg);
    busytoggle(0);
    return;
end
% play  ----------------------------------------------------------------
if strcmp(action,'play')
    rfs = 44100;
    % fetch
    fetchparams;
    mapL = gh.prm.mapL; mapH = gh.prm.mapH;
    fs = gh.fs; tvec = gh.tvec;
    freqtrace = gh.freqtrace; amptrace = gh.amptrace;
    % audible sound synthesis
    snd = soundsynthesis(freqtrace(:,1),amptrace(:,1),tvec,fs,rfs,[mapL mapH]);
    % get range
    rng = max(1,round(gh.xl(1)*rfs)):min(length(snd),round(gh.xl(2)*rfs));
    % play
    ap = audioplayer(snd(rng),rfs);
    playblocking(ap);
    delete(ap);
    set(gh.hpush_swav,'Enable','On');
    % save
    gh.rfs = rfs;
    gh.snd = snd;
    return;
end
% save wav  ----------------------------------------------------------------
if strcmp(action,'swav')
    newfn = [gh.fn(1:end-4) '_syn.wav'];
    [f,p] = uiputfile([gh.pth newfn]);
    audiowrite([p f],gh.snd,gh.rfs);
    return;
end
% long ----------------------------------------------------------------
if strcmp(action,'long')
    % fetch parameter
    fetchparams;
    params = gh.prm;
    timestep = gh.prm.timestep; margin =  gh.prm.margin; durmax = gh.prm.durmax;
    wavflg = gh.prm.wavfileoutput; imgflg = gh.prm.imageoutput;
    fltflg = gh.prm.imagetype; trcflg = gh.prm.traceoutput;
    readsize = gh.prm.readsize; fftsize = gh.fftsize;
    % clear axes
    cla(gh.haxes(1)); cla(gh.haxes(2)); cla(gh.haxes(3));
    % disable
    set(gh.hpush_open,'enable','on');
    set(gh.hpush_dtct,'enable','off');
    set(gh.hpush_save,'enable','off');
    set(gh.hpush_sgmt,'enable','off');
    set(gh.hpush_long,'enable','on');
    % get long wavfile
    busytoggle(1);
    [fn,pth] = uigetfile([gh.pth filesep '*.wav']);
    fp = [pth fn];
    ainfo = audioinfo(fp);
    wavsize = ainfo.TotalSamples;
    fs = ainfo.SampleRate;
    nreadsize= round(readsize*fs);
    nread = ceil(wavsize/nreadsize);
    fvec = [0; (1:(fftsize/2))'/fftsize*fs];
    step = round(params.timestep*fs);
    yl = [0 max(fvec)/1000];
    gh.fs = fs;
    gh.yl = yl;
    % show fs
    set(gh.htext_fsvl,'string',num2str(fs/1000));
    % CSV setting
    [~,prefix,~] = fileparts(fn);
    [sfn,spth] = uiputfile([pth filesep prefix '_dat.csv'],'save file name');
    savefp = [spth sfn];
    fid = fopen(savefp,'wt');
    fprintf(fid,'#,start,end,duration,maxfreq,maxamp\n');
    fclose(fid);
    % segmentation setting
    outp = uigetdir(spth,'output directory');
    % start
    prevn = 0;
    prevlast = 0;
    med = [];
    thresh = [];
    for r=1:nread
        set(gh.htext_fnam,'string',['File: ' fp ' ...  (' num2str(r) '/' num2str(nread) ' blocks)']);
        busytoggle(1);
        % read
        rng = [prevlast+1 min(r*nreadsize,wavsize)];
        if diff(rng)<fftsize*2, break; end
        [wav,fs] = audioread(fp,rng);
        % process
        [mtsp,fltnd,onoffset,onoffsetm,freqtrace,amptrace,maxampval,maxampidx,maxfreq,thresh,med,contflg] = procfun(wav,fs,fftsize,params,med,thresh);
        busytoggle(0);
        dur = diff(onoffset,[],2);
        ronoffset = onoffset+(prevlast+1)/fs;
        ronoffsetm = onoffsetm+(prevlast+1)/fs;
        nstep = size(mtsp,2);
        tvec = ((0:(nstep-1))*step+fftsize/2)/fs;
        rtvec = tvec+(prevlast+1)/fs;
        mtsimrng = median(mtsp(:)) + [0 40];
        fltndimrng = [0 30];
        % save to global 
        gh.tvec = rtvec; gh.fvec = fvec; gh.mtsp = mtsp; gh.fltnd = fltnd;
        gh.onoffset = ronoffset; gh.onoffsetm = ronoffsetm;
        gh.freqtrace = freqtrace; gh.amptrace = amptrace;
        gh.maxampval = maxampval; gh.maxampidx = maxampidx; gh.maxfreq = maxfreq;
        gh.mtsimrng = mtsimrng; gh.fltndimrng = fltndimrng;
        gh.xl = [min(rtvec) max(rtvec)];
        % draw MTSP
        drawmtsp;
        % draw Flattened
        drawfltnd;
        % draw trace
        drawtrace;
        % save CSV
        fid = fopen(savefp,'at');
        if ~isempty(ronoffset)
            for n=1:size(ronoffset,1)
                fprintf(fid,'%d,%f,%f,%f,%f,%f\n',n+prevn,ronoffset(n,1),ronoffset(n,2),dur(n),maxfreq(n)/1000,maxampval(n));
            end
        end
        fclose(fid);
        % segmentation
        if fltflg==1, inputimg = fltnd; imrng = fltndimrng;
        else,         inputimg = mtsp;   imrng = mtsimrng; end
        segfun(prevn+1,outp,prefix,inputimg,imrng,wav,fs,timestep,margin,onoffset,rtvec,amptrace,freqtrace,wavflg,imgflg,trcflg)
        % calc
        prevn = size(onoffset,1)+prevn;
        if contflg==1
            if isempty(onoffset)
                prevlast = rng(2) - durmax * fs;
            else
                prevlast = round(onoffset(end,2)*fs)+prevlast;
            end
        else
            prevlast = rng(2);
        end        
    end
    msgbox('Done!');
    return;
end
% folder (multiple files) ----------------------------------------------------------------
if strcmp(action,'fldr')
    set(gh.htext_fnam,'string','');
    % clear axes
    cla(gh.haxes(1)); cla(gh.haxes(2)); cla(gh.haxes(3));
    % disable
    set(gh.hpush_open,'enable','on');
    set(gh.hpush_dtct,'enable','off');
    set(gh.hpush_save,'enable','off');
    set(gh.hpush_sgmt,'enable','off');
    set(gh.hpush_long,'enable','on');
    % fetch parameters
    fetchparams;
    params = gh.prm;
    timestep = gh.prm.timestep;
    durmax = gh.prm.durmax;
    margin = gh.prm.margin;
    wavflg = gh.prm.wavfileoutput;
    imgflg = gh.prm.imageoutput;
    fltflg = gh.prm.imagetype;
    trcflg = gh.prm.traceoutput;
    fftsize = gh.fftsize;
    % get reading path and files
    pth = uigetdir(gh.pth);
    gh.pth = pth;
    d = dir([pth filesep '*.wav']);
    fns = {d.name};
    nfile = length(fns);
    for f=1:nfile
        % get wavfile
        fp = [pth filesep fns{f}];
        ainfo = audioinfo(fp);
        wavsize = ainfo.TotalSamples;
        fs = ainfo.SampleRate;
        nreadsize= round(params.readsize*fs);
        nread = ceil(wavsize/nreadsize);
        fvec = [0; (1:(fftsize/2))'/fftsize*fs];
        step = round(params.timestep*fs);
        yl = [0 max(fvec)/1000];
        gh.fs = fs;
        gh.yl = yl;
        % show fs
        set(gh.htext_fsvl,'string',num2str(fs/1000));
        % CSV setting
        [~,prefix,~] = fileparts(fns{f});
        sfn = [prefix '_dat.csv'];
        savefp = [pth filesep sfn];
        fid = fopen(savefp,'wt');
        fprintf(fid,'#,start,end,duration,maxfreq,maxamp\n');
        fclose(fid);
        % segmentation setting
        mkdir(pth,prefix);
        outp = [pth filesep prefix];
        % start
        prevn = 0;
        prevlast = 0;
        med = [];
        thresh = [];
        for r=1:nread
            % read
            set(gh.htext_fnam,'string',sprintf('File: %s  ... (%d/%d blocks) [%d/%d files]',fp,r,nread,f,nfile));
            busytoggle(1);
            rng = [prevlast+1 min(r*nreadsize,wavsize)];
            if diff(rng)<fftsize*2, break; end
            [wav,fs] = audioread(fp,rng);
            % process
            [mtsp,fltnd,onoffset,onoffsetm,freqtrace,amptrace,maxampval,maxampidx,maxfreq,thresh,med,contflg] = procfun(wav,fs,fftsize,params,med,thresh);
            dur = diff(onoffset,[],2);
            ronoffset = onoffset+(prevlast+1)/fs;
            ronoffsetm = onoffsetm+(prevlast+1)/fs;
            nstep = size(mtsp,2);
            tvec = ((0:(nstep-1))*step+fftsize/2)/fs;
            rtvec = tvec+(prevlast+1)/fs;
            mtsimrng = median(mtsp(:)) + [0 40];
            fltndimrng = [0 30];
            busytoggle(0);
            % save to global 
            gh.tvec = rtvec; gh.fvec = fvec; gh.mtsp = mtsp; gh.fltnd = fltnd;
            gh.onoffset = ronoffset; gh.onoffsetm = ronoffsetm;
            gh.freqtrace = freqtrace; gh.amptrace = amptrace;
            gh.maxampval = maxampval; gh.maxampidx = maxampidx; gh.maxfreq = maxfreq;
            gh.mtsimrng = mtsimrng; gh.fltndimrng = fltndimrng;
            gh.xl = [min(rtvec) max(rtvec)];
            % draw MTSP
            drawmtsp;
            % draw Flattened
            drawfltnd;
            % draw trace
            drawtrace;
            % save CSV
            fid = fopen(savefp,'at');
            if ~isempty(ronoffset)
                for n=1:size(ronoffset,1)
                    fprintf(fid,'%d,%f,%f,%f,%f,%f\n',n+prevn,ronoffset(n,1),ronoffset(n,2),dur(n),maxfreq(n)/1000,maxampval(n));
                end
            end
            fclose(fid);
            % segmentation
            if fltflg==1, inputimg = fltnd; imrng = fltndimrng;
            else,         inputimg = mtsp;   imrng = mtsimrng; end
            segfun(prevn+1,outp,prefix,inputimg,imrng,wav,fs,timestep,margin,onoffset,rtvec,amptrace,freqtrace,wavflg,imgflg,trcflg);
            % calc
            prevn = size(onoffset,1)+prevn;
            if contflg==1
                if isempty(onoffset)
                    prevlast = rng(2) - durmax * fs;
                else
                    prevlast = round(onoffset(end,2)*fs)+prevlast;
                end
            else
                prevlast = rng(2);
            end        
        end
    end
    msgbox('Done!');
    return;
end
% zoom  -----------------------------------------------------
if strcmp(action,'z')
    set(gcf,'Pointer','crosshair')
    waitforbuttonpress;
    p1 = get(gca,'currentpoint');
    rbbox;
    p2 = get(gca,'currentpoint');
    newxl = sort([p1(1,1) p2(1,1)]);
    newyl = sort([p1(1,2) p2(1,2)]);
    set(gcf,'Pointer','arrow')
    gh.xl = newxl;
    gh.yl = newyl;
    set([gh.haxes],'xlim',gh.xl,'ylim',gh.yl);
    return;
end
% zoom in -----------------------------------------------------
if strcmp(action,'i')
    xl = get(gca,'XLim');
    center = xl(1)+diff(xl)/2;
    newxl = center + [-diff(xl)/2 diff(xl)/2]*0.75;
    gh.xl = newxl;
    set([gh.haxes],'xlim',gh.xl);
    return;
end
% zoom out -----------------------------------------------------
if strcmp(action,'o')
    xl = get(gca,'XLim');
    center = xl(1)+diff(xl)/2;
    newxl = center + [-diff(xl)/2 diff(xl)/2]*1.25;
    gh.xl = newxl;
    set([gh.haxes],'xlim',gh.xl);
    return;
end
% zoom in veritcal --------------------------------------------------
if strcmp(action,'u')
    yl = get(gca,'ylim');
    center = yl(1)+diff(yl)/2;
    newyl = center + [-diff(yl)/2 diff(yl)/2]*0.75;
    gh.yl = newyl;
    set([gh.haxes],'ylim',gh.yl);
    return;
end
% zoom out vertical ---------------------------------------------------
if strcmp(action,'j')
    yl = get(gca,'ylim');
    center = yl(1)+diff(yl)/2;
    newyl = center + [-diff(yl)/2 diff(yl)/2]*1.25;
    gh.yl = newyl;
    set([gh.haxes],'ylim',gh.yl);
    return;
end
% slide forward ----------------------------------------------
if double(action)==29 % ->
    xl = get(gca,'xlim');
    newxl = xl+diff(xl)*0.05;
    gh.xl = newxl;
    set([gh.haxes],'xlim',gh.xl);
    return;
end
% slide back --------------------------------------------------
if double(action)==28 % <-
    xl = get(gca,'xlim');
    newxl = xl-diff(xl)*0.05;
    gh.xl = newxl;
    set([gh.haxes],'xlim',gh.xl);
    return;
end
% slide up ----------------------------------------------
if double(action)==30 % up
    yl = get(gca,'ylim');
    newyl = yl+diff(yl)*0.1;
    gh.yl = newyl;
    set([gh.haxes],'ylim',gh.yl);
    return;
end
% slide down --------------------------------------------------
if double(action)==31 % down
    xl = get(gca,'ylim');
    newyl = xl-diff(xl)*0.1;
    gh.yl = newyl;
    set([gh.haxes],'ylim',gh.yl);
    return;
end
% save PDF file ------------------------------------------------
if strcmp(action,'p')
    [f,p] = uiputfile([gh.pth filesep '*.pdf']);
    set(gcf,'renderer','painters','PaperOrientation','landscape','PaperUnits','normalized','PaperPosition',[0 0 1 1])
    print([p f],'-dpdf');
    return;
end
% help --------------------------------------------------
if strcmp(action,'h')
    msg = {
    'i ... zoom in horizontal'
    'o ... zoom out horizontal'
    'u ... zoom in vertical'
    'j ... zoom out vertical'
    'right ... slide to right'
    'left ... slide to left'
    'up ... slide to upper'
    'down ... slide to lower'
    'p ... save PDF file'
    'h ... help'
    };
    msgbox(msg,'Key Command List');    
    return;
end
% for other things --------------------------------------
if length(action)==1
    return;
end
% for internal functions --------------------------------------
if length(action)>5
    varargout{1} = eval(['@' action]);
    return;
end
% ////////////////////////////////////////////////////////////////////////
function drawmtsp
global gh;

set(0,'CurrentFigure',gh.hfig);
set(gh.hfig,'CurrentAxes',gh.haxes(1));
imagesc(gh.tvec,gh.fvec/1000,gh.mtsp);  axis xy
set(gca,'tickdir','out'); box off
ylabel('Frequency (kHz)');
colormap(flipud(gray));
caxis(gh.mtsimrng);
set(gh.haxes,'xlim',gh.xl,'ylim',gh.yl);
drawnow;
% ////////////////////////////////////////////////////////////////////////
function drawfltnd
global gh;
% fetch
fetchparams;    
freqmin = gh.prm.freqmin;
freqmax = gh.prm.freqmax;
fvec = gh.fvec;
tvec = gh.tvec;
fltnd = gh.fltnd;
fltndimrng = gh.fltndimrng;
% draw
set(0,'CurrentFigure',gh.hfig);
set(gh.hfig,'CurrentAxes',gh.haxes(2));
imagesc(tvec,fvec/1000,fltnd); axis xy
set(gca,'tickdir','out'); box off;
caxis(fltndimrng);
hls = line([min(tvec) max(tvec); min(tvec) max(tvec)]',[freqmax freqmax;freqmin freqmin]'/1000); set(hls,'color','r','linestyle','--');
ylabel('Frequency (kHz)');
set(gh.haxes,'xlim',gh.xl,'ylim',gh.yl);
drawnow;
% ////////////////////////////////////////////////////////////////////////
function drawtrace
global gh;
% fetch
fvec = gh.fvec;
tvec = gh.tvec;
fltnd = gh.fltnd;
fltndimrng = gh.fltndimrng;
freqtrace = gh.freqtrace;
onoffset = gh.onoffset;
onoffsetm = gh.onoffsetm;
maxampidx = gh.maxampidx;
maxfreq = gh.maxfreq;
margin = gh.prm.margin;
fs = gh.fs;
% draw
set(0,'CurrentFigure',gh.hfig);
set(gh.hfig,'CurrentAxes',gh.haxes(3));
imagesc(tvec,fvec/1000,fltnd); axis xy
if ~isempty(onoffset)
    nonsylzone = [[0 onoffsetm(1,1)]; [onoffsetm(1:end-1,2) onoffsetm(2:end,1)]; [onoffsetm(end,2) max(tvec)]];
    temp = [[onoffset(:,1)-margin onoffset(:,1)]; [onoffset(:,2) onoffset(:,2)+margin]];
    [~,sid] = sort(temp(:,1));
    margzone = temp(sid,:);
    idx = find((margzone(2:end,1)-margzone(1:end-1,2))>0);
    ons = [margzone(1,1);margzone(idx+1)];
    offs = [margzone(idx,2);margzone(end,2)];
    margzone = [ons offs];
else
    nonsylzone = [0 tvec(end)];
    margzone = [0 0];
end
set(0,'CurrentFigure',gh.hfig);
set(gh.hfig,'CurrentAxes',gh.haxes(3));
hp1 = patch(margzone(:,[1 2 2 1])',repmat([0 0 fs fs]/2000,size(margzone,1),1)',1);
hp2 = patch(nonsylzone(:,[1 2 2 1])',repmat([0 0 fs fs]/2000,size(nonsylzone,1),1)',1);
set(hp1,'linestyle','none','facecolor','k','facealpha',0.1);
set(hp2,'linestyle','none','facecolor','k','facealpha',0.3);
hlf = line(tvec,freqtrace(:,1:3)/1000);
set(hlf,'linestyle','none','marker','.','color','b');
hlm = line(tvec(maxampidx),maxfreq/1000);
set(hlm,'linestyle','none','marker','+','color','r');
colormap(gca,flipud(gray));
caxis(fltndimrng);
set(gca,'tickdir','out'); box off;
ylabel('Frequency (kHz)');
xlabel('Time (s)');
set(gh.haxes,'xlim',gh.xl,'ylim',gh.yl);    
drawnow;
% ////////////////////////////////////////////////////////////////////////
function makegui
global gh
% figure
hfig = figure;
set(hfig,'MenuBar','none','NumberTitle','off','ToolBar','none','Name',gh.version);
set(hfig,'KeyPressFcn','usvseg08r7(get(gcf,''CurrentCharacter''));');
set(hfig,'DeleteFcn',@finishfunc);
% axes
haxes(1) = axes;
haxes(2) = axes;
haxes(3) = axes;
set(haxes,'tickdir','out')
% uicontrol
hpush_open = uicontrol(hfig,'style','pushbutton');
hpush_flat = uicontrol(hfig,'style','pushbutton');
hpush_thrs = uicontrol(hfig,'style','pushbutton');
hpush_dtct = uicontrol(hfig,'style','pushbutton');
hpush_save = uicontrol(hfig,'style','pushbutton');
hpush_sgmt = uicontrol(hfig,'style','pushbutton');
hpush_play = uicontrol(hfig,'style','pushbutton');
hpush_swav = uicontrol(hfig,'style','pushbutton');
hpush_long = uicontrol(hfig,'style','pushbutton');
hpush_fldr = uicontrol(hfig,'style','pushbutton');
hedit_step = uicontrol(hfig,'style','edit');
hedit_fmin = uicontrol(hfig,'style','edit');
hedit_fmax = uicontrol(hfig,'style','edit');
hedit_thre = uicontrol(hfig,'style','edit');
hedit_dmin = uicontrol(hfig,'style','edit');
hedit_dmax = uicontrol(hfig,'style','edit');
hedit_gmin = uicontrol(hfig,'style','edit');
hedit_marg = uicontrol(hfig,'style','edit');
htggl_wavo = uicontrol(hfig,'style','toggle');
htggl_imgo = uicontrol(hfig,'style','toggle');
htggl_imgt = uicontrol(hfig,'style','toggle');
htggl_trac = uicontrol(hfig,'style','toggle');
hedit_read = uicontrol(hfig,'style','edit');
hedit_mapL = uicontrol(hfig,'style','edit');
hedit_mapH = uicontrol(hfig,'style','edit');
htext_fnam = uicontrol(hfig,'style','text');
htext_fslb = uicontrol(hfig,'style','text');
htext_fsvl = uicontrol(hfig,'style','text');
htext_step = uicontrol(hfig,'style','text');
htext_fmin = uicontrol(hfig,'style','text');
htext_fmax = uicontrol(hfig,'style','text');
htext_dmin = uicontrol(hfig,'style','text');
htext_dmax = uicontrol(hfig,'style','text');
htext_gmin = uicontrol(hfig,'style','text');
htext_thre = uicontrol(hfig,'style','text');
htext_marg = uicontrol(hfig,'style','text');
htext_wavo = uicontrol(hfig,'style','text');
htext_imgo = uicontrol(hfig,'style','text');
htext_imgt = uicontrol(hfig,'style','text');
htext_trac = uicontrol(hfig,'style','text');
htext_maps = uicontrol(hfig,'style','text');
htext_read = uicontrol(hfig,'style','text');
% busy mark
hpanl_busy = uipanel(hfig);
set(hpanl_busy,'backgroundcolor',[0.6 0.6 0.6]);
% strings
set(hpush_open,'string','open');
set(hpush_flat,'string','flatten');
set(hpush_thrs,'string','threshold');
set(hpush_dtct,'string','detect');
set(hpush_save,'string','save csv');
set(hpush_sgmt,'string','segment');
set(hpush_play,'string','play');
set(hpush_swav,'string','save');
set(hpush_long,'string','long file');
set(hpush_fldr,'string','multiple files');
set(htext_fnam,'string','');
set(hedit_step,'string',num2str(gh.prm.timestep*1000));
set(hedit_fmin,'string',num2str(gh.prm.freqmin/1000));
set(hedit_fmax,'string',num2str(gh.prm.freqmax/1000));
set(hedit_thre,'string',num2str(gh.prm.threshval));
set(hedit_dmin,'string',num2str(gh.prm.durmin*1000));
set(hedit_dmax,'string',num2str(gh.prm.durmax*1000));
set(hedit_gmin,'string',num2str(gh.prm.gapmin*1000));
set(hedit_marg,'string',num2str(gh.prm.margin*1000));
set(hedit_read,'string',num2str(gh.prm.readsize));
set(hedit_mapL,'string',num2str(gh.prm.mapL/1000));
set(hedit_mapH,'string',num2str(gh.prm.mapH/1000));
set(htext_fslb,'string',{'sampling', '(kHz)'});
set(htext_fsvl,'string','-');
set(htext_step,'string',{'time step','(ms)'});
set(htext_fmin,'string',{'freq min','(kHz)'});
set(htext_fmax,'string',{'freq max','(kHz)'});
set(htext_thre,'string',{'threshold','(SD)'});
set(htext_dmin,'string',{'dur min','(ms)'});
set(htext_dmax,'string',{'dur max','(ms)'});
set(htext_gmin,'string',{'gap min','(ms)'});
set(htext_marg,'string',{'margin','(ms)'});
set(htext_wavo,'string',{'wavfile','output'});
set(htext_imgo,'string',{'image','output'});
set(htext_imgt,'string',{'image','type'});
set(htext_trac,'string',{'trace','output'});
set(htext_maps,'string',{'map','(kHz)'});
set(htext_read,'string',{'read size','(s)'});
set(htext_fnam,'string','File: filename');
% toggle
str1 = {'off','on'};
str2 = {'orig','flat'};
set(htggl_wavo,'value',gh.prm.wavfileoutput,'string',str1{gh.prm.wavfileoutput+1});
set(htggl_imgo,'value',gh.prm.imageoutput,'string',str1{gh.prm.imageoutput+1});
set(htggl_imgt,'value',gh.prm.imagetype,'string',str2{gh.prm.imagetype+1});
set(htggl_trac,'value',gh.prm.traceoutput,'string',str1{gh.prm.traceoutput+1});
set(htggl_wavo,'callback','str={''off'',''on''};    set(gco,''string'',str{get(gco,''value'')+1});');
set(htggl_imgo,'callback','str={''off'',''on''};    set(gco,''string'',str{get(gco,''value'')+1});');
set(htggl_imgt,'callback','str={''orig'',''flat''}; set(gco,''string'',str{get(gco,''value'')+1});');
set(htggl_trac,'callback','str={''off'',''on''};    set(gco,''string'',str{get(gco,''value'')+1});');
% callbacks
set(hpush_open,'callback','usvseg08r7(''open'');');
set(hpush_flat,'callback','usvseg08r7(''flat'');');
set(hpush_thrs,'callback','usvseg08r7(''thrs'');');
set(hpush_dtct,'callback','usvseg08r7(''dtct'');');
set(hpush_save,'callback','usvseg08r7(''save'');');
set(hpush_sgmt,'callback','usvseg08r7(''sgmt'');');
set(hpush_play,'callback','usvseg08r7(''play'');');
set(hpush_swav,'callback','usvseg08r7(''swav'');');
set(hpush_long,'callback','usvseg08r7(''long'');');
set(hpush_fldr,'callback','usvseg08r7(''fldr'');');
% text set
set(htext_fnam,'horizontalalign','left');
% set position
set(hfig,'position',[10 80 1000 750]);
hts = [htext_fslb,htext_step,htext_fmin,htext_fmax,htext_thre,htext_dmin,htext_dmax,htext_gmin,htext_marg,htext_wavo,htext_imgo,htext_imgt,htext_trac,htext_read];
hed = [htext_fsvl,hedit_step,hedit_fmin,hedit_fmax,hedit_thre,hedit_dmin,hedit_dmax,hedit_gmin,hedit_marg,htggl_wavo,htggl_imgo,htggl_imgt,htggl_trac,hedit_read];
for n=1:length(hts)
    set(hts(n),'units','normalized','position',[0.02 0.94-(n-1)*0.035 0.06 0.035]);
    set(hed(n),'units','normalized','position',[0.08 0.94-(n-1)*0.035 0.04 0.035]);
end
set(hpush_open,'units','normalized','position',[0.02 0.42 0.10 0.04]);
set(hpush_flat,'units','normalized','position',[0.02 0.38 0.10 0.04]);
set(hpush_thrs,'units','normalized','position',[0.02 0.34 0.10 0.04]);
set(hpush_dtct,'units','normalized','position',[0.02 0.30 0.10 0.04]);
set(hpush_save,'units','normalized','position',[0.02 0.26 0.10 0.04]);
set(hpush_sgmt,'units','normalized','position',[0.02 0.22 0.10 0.04]);
set(htext_maps,'units','normalized','position',[0.02 0.16 0.04 0.04]);
set(hedit_mapL,'units','normalized','position',[0.06 0.16 0.03 0.04]);
set(hedit_mapH,'units','normalized','position',[0.09 0.16 0.03 0.04]);
set(hpush_play,'units','normalized','position',[0.02 0.12 0.05 0.04]);
set(hpush_swav,'units','normalized','position',[0.07 0.12 0.05 0.04]);
set(hpush_long,'units','normalized','position',[0.02 0.06 0.10 0.04]);
set(hpush_fldr,'units','normalized','position',[0.02 0.02 0.10 0.04]);
% filename and axis
set(hpanl_busy,'units','normalized','position',[0.175 0.955 0.02 0.025]);
set(htext_fnam,'units','normalized','position',[0.20 0.96 0.77 0.02]);
set(haxes(1), 'units','normalized','position',[0.20 0.68 0.77 0.25]);
set(haxes(2), 'units','normalized','position',[0.20 0.38 0.77 0.25]);
set(haxes(3), 'units','normalized','position',[0.20 0.08 0.77 0.25]);
% disable
set(hpush_open,'enable','on');
set(hpush_flat,'enable','off');
set(hpush_thrs,'enable','off');
set(hpush_dtct,'enable','off');
set(hpush_save,'enable','off');
set(hpush_sgmt,'enable','off');
set(hpush_play,'enable','off');
set(hpush_swav,'enable','off');
set(hpush_long,'enable','on');
set(hpush_fldr,'enable','on');
% save handles
gh.hfig = hfig;
gh.haxes = haxes;
gh.hpush_open = hpush_open;
gh.hpush_flat = hpush_flat;
gh.hpush_thrs = hpush_thrs;
gh.hpush_dtct = hpush_dtct;
gh.hpush_save = hpush_save;
gh.hpush_sgmt = hpush_sgmt;
gh.hpush_play = hpush_play;
gh.hpush_swav = hpush_swav;
gh.hpush_long = hpush_long;
gh.hpush_fldr = hpush_fldr;
gh.htext_fnam = htext_fnam;
gh.htext_fsvl = htext_fsvl;
gh.hedit_step = hedit_step;
gh.hedit_fmin = hedit_fmin;
gh.hedit_fmax = hedit_fmax;
gh.hedit_thre = hedit_thre;
gh.hedit_dmin = hedit_dmin;
gh.hedit_dmax = hedit_dmax;
gh.hedit_gmin = hedit_gmin;
gh.hedit_marg = hedit_marg;
gh.hedit_read = hedit_read;
gh.hedit_mapL = hedit_mapL;
gh.hedit_mapH = hedit_mapH;
gh.htggl_wavo = htggl_wavo;
gh.htggl_imgo = htggl_imgo;
gh.htggl_imgt = htggl_imgt;
gh.htggl_trac = htggl_trac;
gh.hpanl_busy = hpanl_busy;
% delete function
function finishfunc(~,~)
global gh
prmname = [fileparts(mfilename('fullpath')) filesep 'usvseg_prm.mat'];
fetchparams;
prm = gh.prm;
save(prmname,'prm');
clearvars -global gh;
disp('bye!')
% ////////////////////////////////////////////////////////////////////////
function busytoggle(num)
global gh
if num==1
    set(gh.hpanl_busy,'backgroundcolor',[1 0 0]);
    drawnow;
else
    set(gh.hpanl_busy,'backgroundcolor',[0 1 0]);
    drawnow;
end 

% ////////////////////////////////////////////////////////////////////////
function fetchparams
global gh
prm.timestep = str2num(get(gh.hedit_step,'string'))/1000;
prm.freqmin = str2num(get(gh.hedit_fmin,'string'))*1000;
prm.freqmax = str2num(get(gh.hedit_fmax,'string'))*1000;
prm.threshval = str2num(get(gh.hedit_thre,'string'));
prm.durmin = str2num(get(gh.hedit_dmin,'string'))/1000;
prm.durmax = str2num(get(gh.hedit_dmax,'string'))/1000;
prm.gapmin = str2num(get(gh.hedit_gmin,'string'))/1000;
prm.margin = str2num(get(gh.hedit_marg,'string'))/1000;
prm.wavfileoutput= get(gh.htggl_wavo,'value');
prm.imageoutput = get(gh.htggl_imgo,'value');
prm.imagetype = get(gh.htggl_imgt,'value');
prm.traceoutput = get(gh.htggl_trac,'value');
prm.readsize = str2num(get(gh.hedit_read,'string'));
prm.mapL = str2num(get(gh.hedit_mapL,'string'))*1000;
prm.mapH = str2num(get(gh.hedit_mapH,'string'))*1000;
gh.prm = prm;

% ////////////////////////////////////////////////////////////////////////
function [mtsp,fltnd,onoffset,onoffsetm,freqtrace,amptrace,maxampval,maxampidx,maxfreq,thresh,med,contflg] = procfun(wav,fs,fftsize,params,med,thresh)
timestep = params.timestep;
gapmin = params.gapmin;
durmin = params.durmin;
durmax = params.durmax;
margin = params.margin;
freqmin = params.freqmin;
freqmax = params.freqmax;
% multitaper spec
mtsp = multitaperspec(wav,fs,fftsize,timestep);
% flattening
[fltnd,med] = flattening(mtsp,med);
% threshold calculation with n*sigma (SD) of background noise 
if isempty(thresh)
    thresh = estimatethresh(fltnd,fs,freqmin,freqmax,params.threshval);
end
% thresholding
thrshd = thresholding(fltnd,fs,freqmin,freqmax,thresh);
% onset/offset detection
[onoffset,onoffsetm,~,contflg] = detectonoffset(thrshd,fs,timestep,gapmin,durmin,durmax,margin);
% peak tracking
[freqtrace,amptrace,maxampval,maxampidx,maxfreq] = specpeaktracking(mtsp,fltnd,fs,timestep,freqmin,freqmax,onoffset,margin);

% ////////////////////////////////////////////////////////////////////////
function segfun(startid,outp,prefix,inputimg,imrng,wav,fs,timestep,margin,onoffset,tvec,amptrace,freqtrace,wavflg,imgflg,trcflg)
fftsize = (size(inputimg,1)-1)*2;
step = round(timestep*fs);
onset = onoffset(:,1);
offset = onoffset(:,2);
% just show busy ***************************
im = flipud(uint8((inputimg-imrng(1))/(imrng(2)-imrng(1))*64));
cm = flipud(gray(64));    
for n=1:length(onset)
    rng = [max(round((onset(n)-margin)*fs),1) min(round((offset(n)+margin)*fs),length(wav))];
    rngs = round((rng-fftsize/2)/step);
    rng2 = [max(rngs(1),1) min(rngs(2),size(im,2))];
    % wave write
    if wavflg==1
        fname = sprintf('%s_%04d.wav',prefix,n+startid-1);
        audiowrite([outp filesep fname],wav(rng(1):rng(2)),fs);
    end
    % jpg write
    if imgflg==1
        imseg = im(:,rng2(1):rng2(2));
        fname = sprintf('%s_%04d.jpg',prefix,n+startid-1);
        imwrite(imseg,cm,[outp filesep fname]);
    end
    % trace write
    if trcflg==1
        af = [amptrace(:,1) freqtrace(:,1) amptrace(:,2) freqtrace(:,2) amptrace(:,3) freqtrace(:,3)];
        dat = [tvec(rng2(1):rng2(2))' af(rng2(1):rng2(2),:)];
        % CSV file
        fname = sprintf('%s_%04d.csv',prefix,n+startid-1);
        dlmwrite([outp filesep fname],dat,'precision',6);
    end
end

% ////////////////////////////////////////////////////////////////////////
function [freqtrace,amptrace,maxampval,maxampidx,maxfreq] = specpeaktracking(mtsp,fltnd,fs,timestep,freqmin,freqmax,onoffset,margin)
if isempty(onoffset)
    freqtrace = nan(size(mtsp,2),4);
    amptrace =  nan(size(mtsp,2),4);
    maxampval = [];
    maxampidx = [];
    maxfreq = [];
    return;
end
fftsize = (size(fltnd,1)-1)*2;
step = round(fs*timestep);
% spectral peak saliency
spcsal = spectralsaliency(fltnd);
% get peak and reorganize
bandwidth = 9;
ncandidates = 4;
contmin = 10;
ampthr = 0;
nstep = size(fltnd,2);
fnmin = floor(freqmin/fs*fftsize)+1;
fnmax = ceil(freqmax/fs*fftsize)+1;
onidx = round((onoffset(:,1)-margin)*fs/step); % add margin
offidx = round((onoffset(:,2)+margin)*fs/step); % add margin
onidx(1) = max(1,onidx(1));
offidx(end) = min(nstep,offidx(end));
freqmat = nan(nstep,ncandidates);
ampmat = nan(nstep,ncandidates);
maxampval = nan(length(onoffset(:,1)),1);
maxampidx = ones(length(onoffset(:,1)),1);
maxampfreq = nan(size(onoffset,1),1);
for n=1:size(onoffset,1)
    idx = onidx(n):offidx(n);
    [peakfreq,peakamp] = searchpeak(spcsal(fnmin:fnmax,idx),fltnd(fnmin:fnmax,idx),ncandidates,bandwidth);
    [peakfreqsg,peakampsg] = segregatepeak(peakfreq+fnmin-1,peakamp,contmin,ampthr);
    freqmat(idx,:) = peakfreqsg;
    ampmat(idx,:) = peakampsg;
    if any(~isnan(peakampsg(:)))~=0
        [mvC,miC] = max(peakampsg,[],2);
        [~,miR] = max(mvC);
        maxampidx(n) = miR+idx(1)-1;
        maxampfreq(n) = peakfreqsg(miR,miC(miR));
        if ~isnan(maxampfreq(n))
            maxampval(n) = mtsp(round(maxampfreq(n)),maxampidx(n));
        end
    end
end
freqtrace = (freqmat-1)/fftsize*fs;
amptrace = ampmat;
maxfreq = (maxampfreq-1)/fftsize*fs;

% ////////////////////////////////////////////////////////////////////////
function snd = soundsynthesis(freq,amp,tvec,fs,rfs,freqmap)
% time vector
rt = (1:round(rfs*max(tvec)))'/rfs;
% process frequency
nid = ~isnan(freq);
npf = freq(nid);
npf = [mean(npf); npf; mean(npf)];
nT = [1/fs; tvec(nid)'; max(tvec)];
p = interp1(nT,npf,rt);
pm = p/(fs/2)*(freqmap(2)-freqmap(1))+freqmap(1);
pm(pm<100) = 100; % lower limit
% process amplitude
a2 = amp;
a2(isnan(a2)) = -120;
a3 = interp1(tvec,a2,rt);
a4 = 10.^(a3/20);
afil = filter(ones(128,1)/128,1,[a4;zeros(64,1)]);
afil = afil(65:end);
ampli = 0.2 * afil/max(afil);
ampli(20*log10(afil)<-120) = 0;
% synthesis
omega = 2*pi*pm;
ph = cumsum(omega/rfs);
sig = sin(ph);
snd = sig.*ampli + 2^-16*randn(size(sig));

% ////////////////////////////////////////////////////////////////////////
function [peakfreq,peakamp] = searchpeak(specsaliency,specamp,ncandidates,bandwidth)
num_steps = size(specsaliency,2);
search_range = bandwidth-1;
remove_range = bandwidth*2-1;
peakfreq = nan(num_steps,ncandidates);
peakamp = nan(num_steps,ncandidates);
specsaliency(specsaliency<0) = 0;
for n=1:num_steps
    slice = specsaliency(:,n);
    for c=1:ncandidates
        [~,mi] = max(slice);
        % center of gravity
        rng = max(mi-search_range,1):min(mi+search_range,size(slice,1));
        temp = specsaliency(rng,n);
        peak = sum(temp.*rng')/sum(temp); 
        % store
        peakfreq(n,c) = peak;
        peakamp(n,c) = specamp(mi,n);
        % remove
        idx = max(round(peak)-remove_range,1):min(round(peak)+remove_range,size(slice,1));
        slice(idx) = -Inf;
    end
end
% ////////////////////////////////////////////////////////////////////////
function [peakfreqsg,peakampsg] = segregatepeak(peakfreq,peakamp,conthr,ampthr)
% amplitude thresholding
peakfreq(peakamp<ampthr) = nan;
peakamp(peakamp<ampthr) = nan;
%object segregatin with putting object number
%allow skipping two frames (steps)
distthr = 0.05; % 5 percent: fixed parameter
[nstep,ncand] = size(peakfreq);
objmat = reshape((1:(nstep*ncand)),ncand,nstep)';
objmat(isnan(peakfreq)) = nan;
nskip = 2; % can skip max 2 frames if intermediate framse are NaN
distmat = nan(nstep-3,ncand,nskip+1);
pathmat = nan(nstep-3,ncand,nskip+1);
for n=1:nstep-nskip-1
    for m=1:nskip+1
        temp = abs((1./peakfreq(n,:)'*peakfreq(n+m,:))-1);
        [mv,mid] = min(temp,[],2);
        distmat(n,:,m) = mv;
        pathmat(n,:,m) = mid;
    end
    pm = pathmat;
    pm(distmat>distthr) = nan;
    pm(isnan(distmat)) = nan;
    pp = pm(:,:,1);
    if any(~isnan(pm(n,:,1)))
        pp = pm(:,:,1);
        x = n+1;
    elseif any(~isnan(pm(n,:,2)))
        pp = pm(:,:,2);
        x = n+2;
    elseif any(~isnan(pm(n,:,3)))
        pp = pm(:,:,3);
        x = n+3;
    end
    for m=1:ncand
        if ~isnan(pp(n,m))
            if objmat(x,pp(n,m)) < objmat(n,m)
                val = objmat(x,pp(n,m));
                objmat(objmat==val) = objmat(n,m);
            else
                objmat(x,pp(n,m)) = objmat(n,m);
            end
        end
    end
end
% thresholding
objnum = unique(objmat(:));
objnum = objnum(~isnan(objnum));
peaks2 = peakfreq;
ampmat2 = peakamp;
objmat2 = objmat;
objlen = zeros(length(objnum),1);
objamp = zeros(length(objnum),1);
for n=1:length(objnum)
    idx = find(objmat==objnum(n));
    objlen(n) = length(idx);
    objamp(n) = mean(ampmat2(objmat==objnum(n)));
end
for n=1:length(objlen)
    if objlen(n)<conthr
        objlen(n) = nan;
        peaks2(objmat==objnum(n)) = nan;
        objmat2(objmat==objnum(n)) = nan;
        ampmat2(objmat==objnum(n)) = nan;
    end
end
objnum = objnum(~isnan(objlen));
objamp = objamp(~isnan(objlen));
objlen = objlen(~isnan(objlen));
% sorting
peakfreqsg = nan(size(peaks2));
peakampsg = nan(size(peakamp));
for n=1:nstep
    on = objmat2(n,:);
    oa = nan(length(on),1);
    for m=1:length(on)
        if ~isempty(find(objnum==on(m)))
            oa(m) = objamp(objnum==on(m));
        end
    end
    oa2 = oa;
    oa2(isnan(oa)) = -Inf;
    [~,sid] = sort(oa2,'descend');
    peakfreqsg(n,:) = peaks2(n,sid);
    peakampsg(n,:) = ampmat2(n,sid);
end

% ////////////////////////////////////////////////////////////////////////
function spcsal = spectralsaliency(fltnd)
fftsize = (size(fltnd,1)-1)*2;
tfsp = fftshift(sum(abs(fft(dpsstapers,fftsize)),2));
dtfsp = -diff(tfsp,2); % second-order differential
rng = fftsize/2+(-6:6);
rdtfsp = dtfsp(rng);
salfilt = (rdtfsp-mean(rdtfsp))/std(rdtfsp);
fil = filter(salfilt,1,[fltnd;zeros(6,size(fltnd,2))]);
spcsal = fil(7:end,:);

% ////////////////////////////////////////////////////////////////////////
function [onoffset,onoffsetm,onoffsig,contflg] = detectonoffset(thrshd,fs,timestep,gapmin,durmin,durmax,margin,onoffthresh)
if nargin<8
    onoffthresh = 5; % optimized for multitaper spectrogram
end
fftsize = (size(thrshd,1)-1)*2;
step = round(timestep*fs);
% onset/offset detection
onoff = max(filter(ones(onoffthresh,1),1,thrshd))'>=onoffthresh;
% merge fragmented pieces
ndurmin = round(durmin*fs/step);
f = filter(ones(ndurmin,1)/ndurmin,1,[onoff;zeros(round(ndurmin/2),1)]);
monoff = f(round(ndurmin/2)+1:end)>0.5;
monoff(1) = 0; 
monoff(end) = onoff(end);
onidx = find(diff(monoff)>0)+1;
offidx = find(diff(monoff)<0)+1;
if isempty(onidx)||isempty(offidx)
    onoffset = zeros(0,2);
    onoffsetm = zeros(0,2);
    onoffsig = zeros(size(thrshd,2),1);
    contflg = 0;
    return;
end
offidx(end) = min(offidx(end),length(monoff));
if ~isempty(onidx) && monoff(end)==1
    onidx = onidx(1:end-1);
end
% continuity flag: check if the end of read file is "ON"
contflg = monoff(end);
% gap thresholding
gap = (onidx(2:end)-offidx(1:end-1))*timestep;
gap(end+1) = 0;
gid = find(gap>=gapmin);
if ~isempty(gid)
    onidx = [onidx(1); onidx(gid+1)];
    offidx = [offidx(gid); offidx(end)];
else
    onidx = onidx(1);
    offidx = offidx(end);
end
% syllable duration threholding
dur = (offidx-onidx)/fs*step;
did = find(durmin<=dur & dur<=durmax);
onidx = onidx(did);
offidx = offidx(did);
tvec = ((0:(size(thrshd,2)-1))'*step+fftsize/2)/fs;
onset = tvec(onidx);
offset = tvec(offidx);
if isempty(onset)||isempty(offset)
    onoffset = zeros(0,2);
    onoffsetm = zeros(0,2);
    onoffsig = zeros(size(thrshd,2),1);
    contflg = 0;
    return;
end
% margin addition
onsetm = onset-margin;
offsetm = offset+margin;
% syllables whose margins are overlapped are integrated in onoffsetm but not in onoffset
idx = find((onsetm(2:end)-offsetm(1:end-1))>0);
onsetI = [onset(1);onset(idx+1)];
offsetI = [offset(idx);offset(end)];
onsetm = onsetI-margin;
onsetm(1) = max(1/fs*step,onsetm(1));
offsetm = offsetI+margin;
offsetm(end) = min(max(size(thrshd,2)*step/fs),offsetm(end));    
% output 
onoffset = [onset offset];
onoffsetm = [onsetm offsetm];
% on/off signal
temp = zeros(size(onoff));
onidx2 = round((onset*fs-fftsize/2)/step+1);
offidx2 = round((offset*fs-fftsize/2)/step+1);
temp(onidx2) = 1;
temp(offidx2+1) = -1;
onoffsig = cumsum(temp);

% ////////////////////////////////////////////////////////////////////////
function thrshd = thresholding(fltnd,fs,freqmin,freqmax,thresh)
fftsize = (size(fltnd,1)-1)*2;
fnmin = floor(freqmin/fs*fftsize)+1;
fnmax = ceil(freqmax/fs*fftsize)+1;
mask = zeros(size(fltnd));
mask(fnmin:fnmax,:) = 1;
thrshd = (fltnd>thresh).*mask;

% ////////////////////////////////////////////////////////////////////////
function thresh = estimatethresh(fltnd,fs,freqmin,freqmax,threshval)
fftsize = (size(fltnd,1)-1)*2;
fnmin = floor(freqmin/fs*fftsize)+1;
fnmax = ceil(freqmax/fs*fftsize)+1;
cut = fltnd(fnmin:fnmax,:);
bin = -0.05:0.1:10;
bc = bin(1:end-1)+diff(bin)/2;
h = histcounts(cut(:),bin);
fwhm = bc(find(h<h(1)/2,1)-1)*2;
sigma = fwhm/2.35;
thresh = sigma*threshval;

% ////////////////////////////////////////////////////////////////////////
function [fltnd,med] = flattening(mtsp,med)
% generate flattned spectrogram
if nargin<2
    med = [];
end
% liftering
liftercutoff = 3; % fixed parameter
fftsize = (size(mtsp,1)-1)*2;
cep = fft([mtsp;flipud(mtsp(2:end-1,:))]);
lifter = ones(size(cep));
lifter(1:liftercutoff,:) = 0;
lifter((fftsize-liftercutoff+1):fftsize,:) = 0;
temp = real(ifft(cep.*lifter));
liftered = temp(1:(fftsize/2+1),:);
% median centering
if isempty(med)
    med = median(liftered,2);
end
liftmed = liftered-repmat(med,1,size(liftered,2));
% 5-point median filter on frequency axis
if exist('movmedian')==2
    fltnd = movmedian(liftmed,5);
else % for R2015b or earlier
    pad1 = repmat(liftmed(1,:),2,1);
    pad2 = repmat(liftmed(end,:),2,1);
    padded = [pad1;liftmed;pad2];
    fltnd = zeros(size(liftmed));
    for n=1:size(liftmed,1)
        fltnd(n,:) = median(padded(n+(0:4),:));
    end
end

% ////////////////////////////////////////////////////////////////////////
function [whtnd,med] = whitening(mtps,med)
% generate flattned spectrogram
if nargin<2
    med = [];
end
% median centering
if isempty(med)
    med = median(mtps,2);
end
submed = mtps-repmat(med,1,size(mtps,2));
% 5-point median filter on frequency axis
whtnd = medfilt1(submed,5);

% ////////////////////////////////////////////////////////////////////////
function mtsp = multitaperspec(wav,fs,fftsize,timestep)
% generate multitaper spectrogram
step = round(timestep*fs);
wavlen = length(wav);
tapers = dpsstapers;
ntapers = size(tapers,2);
nsteps = floor((wavlen-fftsize+step)/step);
spgsize = fftsize/2+1;
idx = repmat((1:fftsize)',1,nsteps)+repmat((0:nsteps-1)*step,fftsize,1);
wavslice = wav(idx);
spmat = zeros(spgsize,nsteps,ntapers);
parfor n=1:ntapers % use "parfor" if Parallel Computing Toolbox installed
    ft = fft(wavslice.*repmat(tapers(:,n),1,nsteps),fftsize);
    spmat(:,:,n) = abs(ft(1:(fftsize/2+1),:));
end
mtsp = 20*log10(mean(spmat,3)*sqrt(1/(2*pi*fftsize)));

% ////////////////////////////////////////////////////////////////////////
function hnsp = hanntaperspec(wav,fs,fftsize,timestep)
% generate multitaper spectrogram
step = round(timestep*fs);
wavlen = length(wav);
hanwin = hanning(fftsize);
nsteps = floor((wavlen-fftsize+step)/step);
idx = repmat((1:fftsize)',1,nsteps)+repmat((0:nsteps-1)*step,fftsize,1);
wavslice = wav(idx);
ft = fft(wavslice.*repmat(hanwin,1,nsteps),fftsize);
spmat = abs(ft(1:(fftsize/2+1),:));
hnsp = 20*log10(spmat);

% ////////////////////////////////////////////////////////////////////////
function D = dpsstapers
% DPSS windows for multitaper method
D = dpss(512,3,6);
