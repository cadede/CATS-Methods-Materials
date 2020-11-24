vid = mmread([movieloc movies{n}], [],[starttime endtime],false,true);

% these are catches for some unusual situations, should not be relevant now
readingweird = false;
try TD = vid.totalDuration; catch; TD = 160*60*60; end
if length(vid.frames)<150
    disp(['Video ' movies{n} ' is reading weirdly, checking all frames for information']);
    readingweird = true;
    kk = 0;
    while length(vid.frames)<150 && endtime<TD+30*60 && endtime < 160*60*60;
        clear vid; endtime = endtime + 30;
        vid = mmread([movieloc movies{n}], [],[starttime endtime],false,true);
        if mod(round((endtime-starttime)*2/60)/2,300) == 0; kk = kk+5; disp([num2str(kk) ' "hours" read']);
        end
    end
end
% if starttime>406;
%     disp(3)
% end
if length(vid.frames) == 0
    disp(['Video read ending at total duration ' datestr(frameTimes{movN(n)}(end)-frameTimes{movN(n)}(1),'HH:MM:SS.fff')]);
    return
end


if starttime ~= 0 %double check you didn't get the same start frame
    newf = ftime;
    while newf<=ftime || ftime+(23/24+59/60/24+59.5/60/60/24)<newf
        [newf, flag] = gettimestamp(vid.frames(1).cdata,false);
        fr = median(diff(oframeTimes{movN(n)}));
        if ~flag; oi = datevec(newf); oi = oi(6);
            if round((oi-floor(oi))*1000)/1000 < 0.101 && round((oi-floor(oi))*1000)/1000 > 0.099
                [oi2, flag2] = gettimestamp(vid.frames(2).cdata,false); if ~flag2; oi2 = datevec(oi2); oi2 = oi2(6); end
                [oi3, flag3] = gettimestamp(vid.frames(3).cdata,false); if ~flag3; oi3 = datevec(oi3); oi3 = oi3(6); end
                if ~flag2 && ~flag3 && oi2>0.9+oi && oi2<oi+1 && oi3>oi2+fr*0.5&& oi3<oi2+fr*2.5;
                    newf = newf + 0.9/24/60/60;
                end
            end
        end
        if ~flag %as long as you can read the timestamps okay
            if newf<=ftime || ftime+(23/24+59/60/24+59.5/60/60/24)<newf; vid.frames(1) = [];
                % for vid times, it's complicated because you can't
                % tell where mmread is missing files here.  Assumes
                % that the time stamp was right but it missed the
                % wrong frame (only important for synching the
                % audio)
                if vid.times(1)>oframeTimes{movN(n)}(end); vid.times(end) = []; else vid.times(1) = []; end
                if length(vid.frames) == 0;
                    disp(['Video read ending at total duration ' datestr(frameTimes{movN(n)}(end)-frameTimes{movN(n)}(1),'HH:MM:SS.fff')]);
                    return
                end
            end
        else
            newf = ftime+(vid.times(1)-oframeTimes{movN(n)}(end))/24/60/60;
            if oframeTimes{movN(n)}(end)+fr/100>vid.times(1);
                vid.frames(1) = []; vid.times(1) = [];
                if length(vid.frames) == 0;
                    disp(['Video read ending at total duration ' datestr(frameTimes{movN(n)}(end)-frameTimes{movN(n)}(1),'HH:MM:SS.fff')]);
                    return
                end
            end
        end
    end
end
ovid = [0 vid.times]; if ~isempty(oframeTimes{movN(n)}); ovid(1) = oframeTimes{movN(n)}(end); end
oframeTimes{movN(n)} = [oframeTimes{movN(n)} vid.times];
if starttime == 0; stillfirst = true; end
badframes = []; ctime = nan; sketchy = [];
fr = median(diff(vid.times))/24/60/60;
if DAY == 0; changeday = false; else changeday = true; end
totalflags = 0;
for k = 1:length(vid.frames)
    [ftime, flag] = gettimestamp(vid.frames(k).cdata,false);
    trygray = true;
    if k == 1 && starttime >3; lasttime = frameTimes{movN(n)}(end); kk = 0; %check if there is some big gap more than 0.5 seconds
    elseif k>1; kk = find(~isnan(vid.times(1:k-1)),1,'last'); lasttime = vid.times(kk); if isempty(lasttime); lasttime = frameTimes{movN(n)}(end); kk = 0; end
    else kk = 1; lasttime = ftime; trygray = false;
    end
    kk = k-kk;
    grayThresh = 150:20:230; iv = 1;
    while trygray && ~stillfirst
        if (ftime-lasttime)*24*60*60<0 || (ftime-lasttime)*24*60*60>kk/30+0.5 % if there is a big gap, rerun gettime stamp adding gray
            [ftime, flag] = gettimestamp(vid.frames(k).cdata,false,grayThresh(iv));
        else trygray = false; continue;
        end
        iv = iv+1; if iv>length(grayThresh); [ftime, flag] = gettimestamp(vid.frames(k).cdata,false); trygray = false; end
    end
    if ~flag; oi = datevec(ftime); oi = oi(6);
        if round((oi-floor(oi))*1000)/1000 < 0.101 && round((oi-floor(oi))*1000)/1000 > 0.099
            try
                [oi2, flag2] = gettimestamp(vid.frames(k+1).cdata,false); if ~flag2; oi2 = datevec(oi2); oi2 = oi2(6); end
                [oi3, flag3] = gettimestamp(vid.frames(k+2).cdata,false); if ~flag3; oi3 = datevec(oi3); oi3 = oi3(6); end
                if ~flag2 && ~flag3 && oi2>0.9+oi && oi2<oi+1 && oi3>oi2+1/150 && oi3<oi2+1;
                    newf = newf + 0.9/24/60/60;
                end
            catch
                [oi2, flag2] = gettimestamp(vid.frames(k-1).cdata,false); if ~flag2; oi2 = datevec(oi2); oi2 = oi2(6); end
                if ~flag2 && oi2>oi+0.9-2.5*fr*24*60*60 && oi2<oi+0.9
                    ftime = ftime + 0.9/24/60/60;
                end
            end
        end
    end
    if (flag && ~stillfirst)
        badframes = [badframes; k];
        vid.times(k) = nan;
        if k == length(vid.frames)
            bb = find(~isnan(vid.times),1,'last');
            ftime = vid.times(bb)+(k-bb)*fr;
        end
        continue;
    end
    if flag || (k<15 && starttime == 0 && length(vid.times)>16); totalflags = totalflags+1; continue; end %sometimes first few frames are dark and get read as wrong time stamps, but if the video is really really short, use what you've got.
    if (~stillfirst && ftime + DAY<ctime) || (~stillfirst && abs(ftime+DAY-ctime)>2.2*fr) %mark as sketchy
        if ftime+0.9999+DAY<ctime && ~changeday; DAY = DAY+1; changeday = true;
        else
            sketchy = [sketchy; k];
        end
        vid.times(k) = ftime + DAY;
        ctime = ftime + DAY;
        continue;
    end
    if ~stillfirst && (abs(ftime + DAY -ctime)>1/24/60/60||ftime + DAY <ctime) % if something skips or is out of order, just stop
        error(['bad time stamp read in video ' num2str(movN(n)) ' at frame ' num2str(k) ' in currently loaded portion of video']);
    end
    vid.times(k) = ftime + DAY;
    ctime = ftime + DAY;
    if stillfirst
        if readingweird
            vid.times(1:k-1) = ftime + DAY - fr*(k-1):fr:ftime+DAY-fr;
            readingweird = false;
            stillfirst = false;
        else
            vid.times(1:k-1) = ftime + DAY -(ovid(k+1)-ovid(2:k))/24/60/60;
            stillfirst = false;
        end
    end
end
if totalflags == length(vid.times); disp(['Video ' num2str(movN(n)) ' has no readable timestamps.  Moving to bad movies folder']); vidDN(movN(n)) = nan; badmovie = true; return; else badmovie = false; end
if ~isempty(sketchy)
    [s,e] = consec(sketchy);
    s = s-1; if s(1) == 0; s(1) = 1; end; sketchy = unique([sketchy; s']);
    for k = 1:length(s)
        oi = (s(k):e(k))';
        if e(k)+1>length(vid.times)
            badframes = [badframes; oi]; vid.times(oi) = nan;
        else
            nextgood = e(k)+find(~ismember(e(k)+1:e(k)+1000,[badframes;sketchy]),1);
            if nextgood>length(vid.frames); badframes = [badframes; oi]; vid.times(oi) = nan; continue; end
            for kk = 1:length(oi)
                DIF = diff([vid.times(oi(kk)) vid.times(nextgood)]) - fr*(nextgood-oi(kk));
                if  DIF>1.2*fr || DIF <-1.2*fr
                    badframes = [badframes; oi(kk)]; vid.times(oi(kk)) = nan;
                end
                %                         oi = oi(abs(diff(vid.times(s(k):e(k)+1)))>2.2*fr | isnan(vid.times(oi)));
            end
        end
    end
    badframes = unique(badframes);
end
x = (1:length(vid.times))'; y = vid.times'; x(isnan(y)) = []; y (isnan(y)) = [];

% this now fixes all bad gaps after the whole video is read
% if length(x)>1 && (any(diff(y)<0) || any(diff(y)>3/24/60/60)) % last two conditions only runs this if there is a bad frame (the negative diff y) or if there is a huge jump bigger than 3 seconds (else the jump is probably legit)
%     [obj,gof] = fit(x,y,'poly1');
%     newy = x*obj.p1+obj.p2;
%     if gof.rsquare<0.999
%         [obj2,gof2] =  fit(x(y>newy),y(y>newy),'poly1');
%         newy2 = x*obj2.p1+obj2.p2;
%         if gof2.rsquare>0.999; newy = newy2;
%         else
%             [obj3,gof3] =  fit(x(y<newy),y(y<newy),'poly1');
%             newy3 = x*obj3.p1+obj3.p2;
%             if gof3.rsquare>0.999; newy = newy3;
%             elseif length(vid.times)<50 && (max(vid.times)-min(vid.times))*24*60*60<2; disp(['Video ' num2str(movN(n)) ' seems very short with a jump in frame.  See plot for more info but saving as is.']);
%                 resid = diff([newy y],[],2); bads = find(abs(resid)>fr*2.2); figure; plot(x,y,x,newy,x(bads),y(bads),'rx'); title(['Video ' num2str(movN(n)) ]); xlabel('Frame #'); set(gca,'yticklabel',datestr(get(gca,'ytick'),'HH:MM:SS.fff')); legend('orig','regression')
%             elseif gettagnum(datafile)>=40
%                 resid = diff([newy y],[],2);
%                 bads = find(abs(resid)>fr*2.2);
%                 figure; plot(x,y,x,newy,x(bads),y(bads),'rx');
%                 set(gca,'yticklabel',datestr(get(gca,'ytick'),'HH:MM:SS.fff'));
%                 title (['Video ' movies{n}]);
%                 II = vid.times(~isnan(vid.times));
%                 if ~any(diff(II)<0); disp('All increasing'); newy = y;
%                 else
%                     II = find(~isnan(vid.times));
%                     disp(['Bad Points: ' num2str(II(find(diff(vid.times(II))<0))) ', or:']);
%                     disp(['Bad Points: ' num2str(II(find(diff(vid.times(II))<0)+1))]);
%                     disp('Check Video');
%                     error('Check the comments under this line in readwirelessvideo for instructions on how to correct for this error.');
%                 end
%                 % Tag 40 reads funny sometimes.  This is set up to deal with slightly different frame rates, but if there are errors both higher
%                 % and lower than a median line, you have to remove those manually.  Comment out line 221 (with the error), put a breakpoint at line 220, whenever it breaks, run:
%                 % I = [92:96 196:226]; vid.times(I) = nan; badframes = unique([badframes; I']);
%                 % but with "92:96 196:226" replaced with the bad frames from the x-axis of the associated graph.
%                 % Then run:
%                 % x = (1:length(vid.times))'; y = vid.times'; x(isnan(y)) = []; y (isnan(y)) = []; [obj,gof] = fit(x,y,'poly1');   newy = x*obj.p1+obj.p2; gof.rsquare<0.999
%                 % if the result is false (0), continue.  If the result is true, run:
%                 % [obj2,gof2] =  fit(x(y>newy),y(y>newy),'poly1'); newy2 = x*obj2.p1+obj2.p2; [obj3,gof3] = fit(x(y<newy),y(y<newy),'poly1'); newy3 = x*obj3.p1+obj3.p2; if gof2.rsquare>0.999; newy = newy2; elseif gof3.rsquare>0.999; newy = newy3; else error('You''ll have to figure out a new solution to this error...'); end
%             else error('Some frames appear to be read with smaller or larger time stamps and I can''t tell which one automatically, suggest running the plot on line 237, and perhaps adding this tag number to the else in line 208 if this looks like a one-off error');
%             end
%         end
%     end
%     resid = diff([newy y],[],2);
%     bads = find(abs(resid)>fr*2.2);
%     %     figure; plot(x,y,x,newy,x(bads),y(bads),'rx');
%     vid.times(x(bads)) = nan;
%     badframes = unique([badframes; x(bads)]);
% end

if length(badframes) == length(vid.times) && length(vid.times)<50 %if you don't have enough timestamps
    error(['All frames bad in video ' num2str(movN(n)) ', could use data.Camera variable?']);
    %         load([dataloc datafile(1:end-3) 'mat']);
    %         if ~exist('UTC','var'); UTC = input('UTC offset? '); end
    %         aa = find(data.Time+data.Date<=frameTimes{movN(n)-1}(end)/24/60/60+vidDN(movN(n)-1)+UTC/24+data.Date(1),1,'last'); % find the last camon time from the last video
    %         fs = round(1/(median(diff(data.Time(100:150)))*24*60*60));
    %         a1 = data.Time(find(data.Camera(aa+round(fs/2):end),1)+aa+round(fs/2)-1)-UTC/24;
    %         if a1<frameTimes{movN(n)-1}(end)/24/60/60+vidDN(movN(n)-1); a1 = a1+1; end
    %         vid.times = oframeTimes{movN(n)}/24/60/60+a1;
    %         badframes = [];
end

% fix this at the end
% for k = 1:length(badframes)
%     nextgood = find(~isnan(vid.times(badframes(k)+1:end)),1);
%     if isempty(nextgood); vid.times(badframes(k)) = vid.times(badframes(k)-1)+mean(diff(vid.times(badframes(k)-4:badframes(k)-1)));
%     else try DIF = diff([vid.times(badframes(k)-1) vid.times(badframes(k)+nextgood)])/(nextgood+1);
%             vid.times(badframes(k)) = vid.times(badframes(k)-1) + DIF;
%         catch
%             vid.times(badframes(k)) = vid.times(badframes(k)+nextgood)-fr*nextgood;
%         end
%     end
%     disp(['Bad frame read at frame ' num2str(badframes(k)+length(frameTimes{movN(n)})) ', replaced with ' datestr(vid.times(badframes(k)),'HH:MM:SS.fff')]);
% end
