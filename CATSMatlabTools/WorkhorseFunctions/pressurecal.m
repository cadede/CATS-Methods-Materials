function [Depth,CAL] = pressurecal(data,DN,CAL,nopress,ofs,df,tagon,pHz)
    
try pconst = CAL.pconst; pcal = CAL.pcal;
catch
    CAL.pconst = 0; CAL.pcal = 1;
    pconst = 0; pcal = 1;
end
try pressTemp = data.TempDepthInternal; catch; try pressTemp = data.Temp1; catch; pressTemp = data.Temp; end; end

if ~nopress %&& tagnum>5

    [~,pc] = fix_pressure((data.Pressure(tagon)-pconst)*pcal,pressTemp(tagon),ofs); %cal procedure from animaltags.org
    %     Depth = nan(size(data.Pressure));
    %     Depth(tagon) = pT;
    %     Depth = decdc(calDepth((data.Pressure-pconst)*pcal,data.Temp,tagon,ofs),df);
    %     [~,peakmag] = peakfinder(Depth,3,min(Depth)+6,-1); Depth = Depth-max(min(peakmag),max(peakmag)-2);
    % else %Depth = decdc(calDepth(data.Pressure,data.Temp,tagon,ofs),df);
    %     [~,pc] = decdc(fix_pressure(data.Pressure(tagon),data.Temp1q(tagon),ofs),df); %Mark Johnson's method
    % end
    Depth = (data.Pressure-pconst)*pcal+polyval([pc.tcomp,pc.poly(2)],pressTemp-pc.tref);
    
    Depth = decimateM(Depth,ofs,pHz,df,length(DN),'pHz');
%     if length(Depth) == length(DN)+1; Depth(end) = []; end
%     if abs(round(ofs/pHz)-ofs/pHz)>.01; error('pHz does not divide evenly into original sample rate'); end
%     Depth = Depth(1:ofs/pHz:end);
%     if abs(round(df/(ofs/pHz)) - df/(ofs/pHz)) > .01; error('decimation factor does not divide evenly into pHz'); end
%     Depth = decdc(Depth,df/(ofs/pHz));
    
    figure(4); clf;
    oDepth = interp2length(decdc((data.Pressure-pconst)*pcal,df),ofs/df,ofs/df,length(DN));
    plot(DN,oDepth,'b',DN,Depth,'g--');
    set(gca,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'),'ydir','rev');
    legend('Bench cal','In situ cal (animaltags.org method)');
    title('Examine data, then look at main screen to choose calibration method');
    % to add an additional offset, run this script where offset = what you are trying to add or subtract off the observed values:
    %     pc.poly(2) = pc.poly(2)-offset; Depth = decdc((data.Pressure-pconst)*pcal+polyval([pc.tcomp,pc.poly(2)],tTempI-pc.tref),df);
    %      figure(4); clf;  plot(data.Date+data.Time+timedif/24,data.Pressure,DN,Depth); set(gca,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS')); legend('Orig','Calibrated Pressure (Johnson Method)');
    z1 = zoom(4); %z2 = zoom(sp2); z3 = zoom(sp3);% z4 = zoom(s4);
    set(z1,'enable','on','Motion','both');
    Pchoice = input('bench cal (blue) = 1, in situ cal (green dashes) = 2, add an offset = 3. Which do you want to use? (Enter 1 or 2) ');
    switch Pchoice
        case 1
            Depth = oDepth;
            pc = [];
            disp('bench cal chosen');
        case 2
            disp('in situ cal chosen');
        case 3
            P2 = input('Choose your baseline: bench cal (blue) = 1, in situ cal (green) = 2: ');
            offset = input('Linear offset (what you want to add to (+) or subtract from (-) the observed values: ');
            switch P2
                case 2
                    pc.poly(2) = pc.poly(2)+offset; Depth = interp2length(decdc((data.Pressure-pconst)*pcal+polyval([pc.tcomp,pc.poly(2)],pressTemp-pc.tref),df),ofs/df,ofs/df,length(DN));
                    disp('in situ cal chosen');
                case 1
                    pconst = pconst - offset/pcal;
                    Depth = decdc((data.Pressure-pconst)*pcal,df);
                    pc = [];
                    disp('bench cal chosen');
            end
            oDepth = interp2length(decdc((data.Pressure-CAL.pconst)*CAL.pcal,df),ofs/df,ofs/df,length(DN));
            figure(4); clf;  plot(DN,oDepth,'b',DN,Depth,'g--'); set(gca,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'),'ydir','rev'); legend('original data','recalibrated pressure');
        otherwise
            disp('No calibration chosen, so in situ cal was used by default');
    end


else
    Depth = interp2length(decdc((data.Pressure-pconst)*pcal,df),ofs/df,ofs/df,length(DN)); 
    pc = [];
end
CAL.pconst = pconst;
CAL.pcal = pcal;
CAL.pc = pc;
