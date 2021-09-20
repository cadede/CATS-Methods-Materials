% function list = mmwrite(filename,...options...)
% 
% mmwrite is able to write AVI,WMV,WMA,ASF files.  For AVI files you can 
% choose from the available codecs to compress the audio and video streams.
% For WMV,WMA and ASF the encoding defaults to Windows Media 9 44100Hz 
% 16bit stereo 98% quality for the audio and Windows Media 9 Video with 98
% quality.  The quality can be specified for both audio and video.  
% Surround sound only seems to work with AVI and multi-pass encoding is not
% supported.  Writing any other file type is not supported.  This uses the
% Windows DirectX infrastructure to generate the movie, so other OSs are
% out of luck.
%
% INPUT:
%   filename:  This must be the first parameter and specifies the filename
%     to write.
%
%   video structure:  The video structure matches the output of mmread.  At
%     a minimum it must have 4 fields "frames", "times", "height" and 
%     "width".  The "frames" field must be a struct array with a field 
%     "cdata" that contains the raw frame data encoded as height by width 
%     by color(3) as UINT8s.  The "times" field contains the time stamps of
%     the data stored in frames.  "times" and "frames.cdata" must be the 
%     same length.  The reason we use times, is to support variable frame 
%     rates as several codecs generate.
%
%   audio structure:  The audio structure matches the output of mmread.  At
%     a minimum it must have 3 fields "data", "rate" and "times".  The
%     "data" field is a matrix nrSamples by nrChannels (the same format was
%     wavread/wavplay).  The field "rate" is the sampling rate of the data
%     in Hz, eg. 44100.  The field "frames" is used to specify the time
%     that audio should start, the rest of the time is extrapolated based
%     upon the "rate" and the nrSamples.
%
%   AVI config structure:  
%     videoCompressor   Specify which video compressor/codec to use.  Use
%       'ListAviVideoEncoders' to determine what are valid codecs on your
%       machine.
%     audioCompressor   Specify which audio compressor/codec to use.  Use
%       'ListAviAudioEncoders' to determine what are valid codecs on your
%       machine.
%
%   WMV/WMA/ASF config structure:  can have any of the following fields
%     videoQuality  the quality of the video, between 0 and 100 default 98.
%     audioQuality  the quality of the audio, between 0 and 100 default 98.
%     outputHeight  the height of the video to be generated.
%     outputWidth   the width of the video to be generated.
%     outputFrameRate   the frame rate of the video to be generated.
%     prxFile   Specify a custom encoding file.  The settings here overwide
%       the all of the other config options.  To create a custom file, use
%       Windows Media Encoder and use the Export feature of the Compression
%       tab.
%
%   'ListAviVideoEncoders': Use this option to list availble video encoders
%     Eg. list = mmwrite('','ListAviVideoEncoders');
%
%   'ListAviAudioEncoders': Use this option to list availble audio encoders
%     Eg. list = mmwrite('','ListAviAudioEncoders');
%
%   'Continue':  Keep the encoding going so that more data can be added by
%     a later call to mmwrite.  Defaults to false.  To have a usable output
%     file, you must later call mmwrite with only the 'Initialized' option.
%
%   'Initialized':  Indicates that mmwrite has already been initialized (by
%     a call with 'Continue') and to just add the data specified.  Warning,
%     the order of audio and video structures must be the same as the
%     first 'Continue' command otherwise the streams will get mixed.
%
%   OUTPUT
%       list:   Only 'ListAviVideoEncoders' and 'ListAviAudioEncoders' have
%         an output, which is the list of encoders installed on the system.
%
%   EXAMPLES:
%
%   % write a simple WMV file with audio and video
%   mmwrite('blah.wmv',audio,video);
%
%   % make a WMA file from the audio taken from another video
%   [video, audio] = mmread('your movie');
%   mmwrite('blah.wma',audio);
%
%   %specify the encoding quality
%   conf.videoQuality = 50;
%   conf.audioQuality = 50;
%   mmwrite('blah.wmv',audio,video,conf);
%
%   % make a video progressively
%   mmwrite('blah.wmv',audio,video,'Continue'); %initialize the movie
%   ...
%   % the "times" fields for both audio and video must start after the last
%   "times" in the previous call to mmwrite.
%   mmwrite('blah.wmv',audio2,video2,'Continue','Initialized'); %don't initialize or stop
%   ...
%   mmwrite('blah.wmv',audio3,video3,'Initialized'); %don't initialize and STOP
%
%   % make an AVI with custom compressors.
%   audioList = mmwrite('','ListAviAudioEncoders');
%   videoList = mmwrite('','ListAviVideoEncoders');
%   if ~any(ismember(videoList,'ffdshow video encoder'))
%     conf.videoCompressor = 'Cinepak Codec by Radius'; % default to this if ffdshow isn't installed...
%   else
%     conf.videoCompressor = 'ffdshow video encoder';
%   end
%   mmwrite('balh.avi',audio,video,conf);
%
%   % write just a subsample (the time range 10 to 20s) of a movie using mmread
%   [video, audio] = mmread('your file',[],[10 20]);
%   % subtract 10 seconds off the time stamps so that the audio and video 
%   will start at the beginning of the movie.
%   video.times = video.times - 10; 
%   audio.times = audio.times - 10; 
%   mmwrite('blah.wmv',audio,video); % make the movie...
%
%   % multiple audio and video streams are supported
%   % I'm not sure why you would ever want to do this though...
%   mmwrite('blah.wmv',video1,video2,audio1,audio2);
%
%   % make a movie ourselves...
%   f=figure;
%   x = -pi:.1:pi;
%   set(gca,'xlim',[-length(x) length(x)],'ylim',[-length(x) length(x)]);
%   for i=1:length(x)
%       patch(sin(x)*i,cos(x)*i,[abs(cos(x(i))) 0 0]);
%       v.frames(i) = getframe(gca);
%       v.times(i) = i/30; % display at 30fps
%   end
%   v.width=size(v.frames(1).cdata,2);
%   v.height=size(v.frames(1).cdata,1);
%   close(f);
%   mmwrite('blah.wmv',v);
