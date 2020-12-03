function data = datcomimport( files, varargin )
%  DATCOMIMPORT Bring DATCOM file into MATLAB.
%   AERO = DATCOMIMPORT( FILE ) takes a filename as a string or cell array
%   of filenames as strings, FILE, and imports aerodynamic data from FILE
%   into a cell array of structures, AERO. Prior to reading DATCOM file
%   values are initialized to 99999, in order to show when there is not a
%   full set of data for the DATCOM case.
%
%   AERO = DATCOMIMPORT( FILE, USENAN ) is an alternate method allowing the
%   specification of using NaN or zero to replace data points where no
%   DATCOM methods exist or where the method is not applicable.  The
%   default value for USENAN is true.
%
%   AERO = DATCOMIMPORT( FILE, USENAN, VERBOSE ) is an alternate method
%   allowing the an additional specification of how the status of the
%   DATCOM file read is displayed. The default value for VERBOSE is 2 
%   which displays a waitbar.  Other options are 0 which displays no
%   information, and 1 which displays text to the MATLAB command window.
%
%   AERO = DATCOMIMPORT( FILE, USENAN, VERBOSE, FILETYPE ) is an alternate
%   method allowing the an additional specification of which type of DATCOM
%   file is read. The default value for FILETYPE is 6 which reads the
%   for006.dat file output by DATCOM.  Other options are 21 and 42 which
%   read the for021.dat and for042.csv file output by DATCOM 2007 or DATCOM
%   2008.  Note: For FILETYPE = 21, the last entry in AERO is table and
%   breakpoints collected from all of the cases read. 
%
%   The fields of AERO are dependent on the data within the DATCOM file.
%   See documentation for list of possible fields and their definitions.
%
%   Limitation: 
%
%   The operational limitations of Digital DATCOM apply to the data
%   contained in AERO.  For more information on Digital DATCOM limitations,
%   see section 2.4.5 of reference 1 listed in the documentation.
%
%   USAF Digital DATCOM data for wing section, horizontal tail section,  
%   vertical tail section and ventral fin section are not read.
%
%   Examples:
%
%   Read the USAF Digital DATCOM output file astdatcom.out:
%      aero = datcomimport('astdatcom.out')
%
%   Read the USAF Digital DATCOM output file astdatcom.out using zeros to
%   replace data points where no DATCOM methods exist and displaying no
%   status information: 
%      usenan = false;
%      aero = datcomimport('astdatcom.out', usenan, 0 )
%

%   Copyright 2000-2017 The MathWorks, Inc.

%   References: 
%
%   [1] AFFDL-TR-79-3032: "The USAF Stability and Control DATCOM
%   Volume 1, Users Manual."
%    
%   [2] AFRL-VA-WP-TR-1998-3009: "MISSILE DATCOM, Users Manual - 1997 FORTRAN 90 Revision"
%
%   [3] AFRL-RB-WP-TR-2009-3015: "MISSILE DATCOM, Users Manual - 2008 Revision"
%
%   [4] AFRL-RB-WP-TR-2011-3071: "MISSILE DATCOM, Users Manual - 2011 Revision"
%
%   [5] AFRL-RQ-WP-TR-2014-3999: "MISSILE DATCOM, Users Manual - 2014 Revision"

narginchk(1, 4);

if ~ischar( files ) && ~iscell( files ) && ~isstring( files )
    error (message('aero:datcomimport:notCharOrCell'));
end

if iscell( files )
    if ~all(cellfun('isclass',files,'char')) && ~all(cellfun('isclass',files,'string'))
        error (message('aero:datcomimport:notChar'));
    end
end

if ischar( files ) || isstring( files )
    files = { files };
end
[usenan,filetype,verbose] = checkarguments(nargin, varargin{:});


switch filetype
    case 6
        % read datcom 006 files
        [casenoout,casedata] = usafdatcom(files,usenan,verbose);
        
        % post-process data
        data = postprocess6();
    case 21
        % read datcom 021 files
        [caseno,data] = read021datcom(files,usenan,verbose);

        % post-process data collect into tables and breakpoints and put into 
        % last structure entry in data
        data = postprocess21(data);
        
    case 42
        % read datcom 042 files
        [caseno,data] = read042datcom(files,usenan,verbose);
        
    otherwise
        error (message('aero:datcomimport:unrecognizedFileType'));
end


%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    function [usenan,filetype,verbose] = checkarguments(nargin, varargin)
        if nargin == 1
            usenan = true;
            verbose = 2;
            filetype = 6;
        elseif nargin == 2
            if ~islogical(varargin{1})
                error(message('aero:datcomimport:notLogical'));
            end
            usenan = varargin{1};
            verbose = 2;
            filetype = 6;
        else
            if ~islogical(varargin{1})
                error(message('aero:datcomimport:notLogical'));
            end
            if ~(isnumeric(varargin{2}) && any(varargin{2} == [0 1 2]))
                error(message('aero:datcomimport:notNumeric'));
            end
            usenan = varargin{1};
            verbose = varargin{2};   
            if nargin == 3
                filetype = 6;
            else
                if ~isnumeric(varargin{3})
                    error(message('aero:datcomimport:notNumericFileType'));
                end
                filetype = varargin{3};
            end
        end
    end
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    function data = postprocess6()
        for k = casenoout:-1:1
            data{k} = casedata{k}; 
        end
    end
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    function data = postprocess21(data)
        % if no data was read return
        if (caseno == 0)
            return;
        end
        
        % collect, sort and find unique arrays of inputs
        mach      = [];
        alt       = [];
        alpha     = [];
        beta      = [];
        total_col = [];
        deriv_col = [];
        delta1    = [];
        delta2    = [];
        delta3    = [];
        delta4    = [];
        for m = 1:caseno
            mach      = [mach data{m}.mach]; %#ok<AGROW>
            alt       = [alt data{m}.alt]; %#ok<AGROW>
            alpha     = [alpha data{m}.alpha]; %#ok<AGROW>
            beta      = [beta data{m}.beta]; %#ok<AGROW>
            total_col = [total_col data{m}.total_col]; %#ok<AGROW>
            deriv_col = [deriv_col data{m}.deriv_col]; %#ok<AGROW>
            delta1    = [delta1; data{m}.config.fin1.delta]; %#ok<AGROW>
            delta2    = [delta2; data{m}.config.fin2.delta]; %#ok<AGROW>
            delta3    = [delta3; data{m}.config.fin3.delta]; %#ok<AGROW>
            delta4    = [delta4; data{m}.config.fin4.delta]; %#ok<AGROW>
        end
        endcaseno = caseno + 1;
        data{endcaseno}.mach      = unique(sort(mach));
        data{endcaseno}.alt       = unique(sort(alt));
        data{endcaseno}.alpha     = unique(sort(alpha));
        data{endcaseno}.nalpha    = length(data{endcaseno}.alpha);
        data{endcaseno}.beta      = unique(sort(beta));
        data{endcaseno}.total_col = unique(sort(total_col));
        data{endcaseno}.deriv_col = unique(sort(deriv_col));
        data{endcaseno}.config.fin1.delta = unique(sortrows(delta1),'rows');
        data{endcaseno}.config.fin2.delta = unique(sortrows(delta2),'rows');
        data{endcaseno}.config.fin3.delta = unique(sortrows(delta3),'rows');
        data{endcaseno}.config.fin4.delta = unique(sortrows(delta4),'rows');
        data{endcaseno}.version = data{1}.version;
        
        nmach   = length(data{endcaseno}.mach);
        nalt    = length(data{endcaseno}.alt);
        nbeta   = length(data{endcaseno}.beta);
        ndelta1 = size(data{endcaseno}.config.fin1.delta,1);
        ndelta2 = size(data{endcaseno}.config.fin2.delta,1);
        ndelta3 = size(data{endcaseno}.config.fin3.delta,1);
        ndelta4 = size(data{endcaseno}.config.fin4.delta,1);
        
        % initialize static coefficients
        data{endcaseno}.cn    = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4);
        data{endcaseno}.cm    = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4);
        data{endcaseno}.ca    = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4);
        data{endcaseno}.cy    = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4);
        data{endcaseno}.cln   = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4);
        data{endcaseno}.cll   = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4);
        
        for m = 1:caseno
            if isfield(data{m},'cmad')
                % initialize dynamic derivatives
                data{endcaseno}.cnad = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4); 
                data{endcaseno}.cmad = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4); 
                data{endcaseno}.cnq  = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4); 
                data{endcaseno}.cmq  = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4); 
                data{endcaseno}.caq  = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4); 
                data{endcaseno}.cyq  = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4); 
                data{endcaseno}.clnq = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4); 
                data{endcaseno}.cllq = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4); 
                data{endcaseno}.cnp  = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4); 
                data{endcaseno}.cmp  = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4); 
                data{endcaseno}.cap  = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4); 
                data{endcaseno}.cyp  = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4); 
                data{endcaseno}.clnp = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4);  
                data{endcaseno}.cllp = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4); 
                data{endcaseno}.cnr  = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4); 
                data{endcaseno}.cmr  = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4); 
                data{endcaseno}.car  = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4); 
                data{endcaseno}.cyr  = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4); 
                data{endcaseno}.clnr = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4); 
                data{endcaseno}.cllr = 99999*ones(data{endcaseno}.nalpha,nmach,nalt,nbeta,ndelta1,ndelta2,ndelta3,ndelta4); 
                break;
            end
        end
        
        for m = 1:caseno
            % find indices in combined N-D matrix
            idxmach = find( data{m}.mach == data{endcaseno}.mach );
            idxalt  = find( data{m}.alt == data{endcaseno}.alt );
            idxbeta = find( data{m}.beta == data{endcaseno}.beta );
            for p = 1:ndelta1
                idxdel1 = all( data{m}.config.fin1.delta == data{endcaseno}.config.fin1.delta(p,:) )*p;
                if idxdel1 > 0
                    break;
                end
            end
            for p = 1:ndelta2
                idxdel2 = all( data{m}.config.fin2.delta == data{endcaseno}.config.fin2.delta(p,:) )*p;
                if idxdel2 > 0
                    break;
                end
            end
            for p = 1:ndelta3
                idxdel3 = all( data{m}.config.fin3.delta == data{endcaseno}.config.fin3.delta(p,:) )*p;
                if idxdel3 > 0
                    break;
                end
            end
            for p = 1:ndelta4
                idxdel4 = all( data{m}.config.fin4.delta == data{endcaseno}.config.fin4.delta(p,:) )*p;
                if idxdel4 > 0
                    break;
                end
            end
            
            for p = 1:data{m}.nalpha
                idxaph  = find( data{m}.alpha(p) == data{endcaseno}.alpha );
                
                % set static coefficients
                data{endcaseno}.cn(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4)  = data{m}.cn(p); 
                data{endcaseno}.cm(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4)  = data{m}.cm(p);
                data{endcaseno}.ca(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4)  = data{m}.ca(p); 
                data{endcaseno}.cy(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4)  = data{m}.cy(p); 
                data{endcaseno}.cln(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4) = data{m}.cln(p); 
                data{endcaseno}.cll(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4) = data{m}.cll(p); 
                
                if isfield(data{m},'cmad')
                    % set dynamic derivatives
                    data{endcaseno}.cnad(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4) = data{m}.cnad(p);
                    data{endcaseno}.cmad(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4) = data{m}.cmad(p);
                    data{endcaseno}.cnq(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4)  = data{m}.cnq(p);
                    data{endcaseno}.cmq(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4)  = data{m}.cmq(p);
                    data{endcaseno}.caq(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4)  = data{m}.caq(p);
                    data{endcaseno}.cyq(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4)  = data{m}.cyq(p);
                    data{endcaseno}.clnq(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4) = data{m}.clnq(p);
                    data{endcaseno}.cllq(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4) = data{m}.cllq(p);
                    data{endcaseno}.cnp(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4)  = data{m}.cnp(p);
                    data{endcaseno}.cmp(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4)  = data{m}.cmp(p);
                    data{endcaseno}.cap(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4)  = data{m}.cap(p);
                    data{endcaseno}.cyp(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4)  = data{m}.cyp(p);
                    data{endcaseno}.clnp(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4) = data{m}.clnp(p);
                    data{endcaseno}.cllp(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4) = data{m}.cllp(p);
                    data{endcaseno}.cnr(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4)  = data{m}.cnr(p);
                    data{endcaseno}.cmr(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4)  = data{m}.cmr(p);
                    data{endcaseno}.car(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4)  = data{m}.car(p);
                    data{endcaseno}.cyr(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4)  = data{m}.cyr(p);
                    data{endcaseno}.clnr(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4) = data{m}.clnr(p);
                    data{endcaseno}.cllr(idxaph,idxmach,idxalt,idxbeta,idxdel1,idxdel2,idxdel3,idxdel4) = data{m}.cllr(p);
                end
            end
        end
        
    end
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
end

%   The fields of AERO are dependent on the data within the DATCOM file.
%   Common fields are the following:
%   case             :a string containing the CASEID. The default value is [].
%   mach             :an array of Mach numbers. The default value is [].
%   alt              :an array of altitudes. The default value is [].
%   alpha            :an array of angle of attacks. The default value is [].
%   nmach            :the number of Mach numbers. The default value is 0.
%   nalt             :the number of altitudes. The default value is 0.
%   nalpha           :the number of angle of attacks. The default value is 0.
%   rnnub            :an array of Reynolds numbers. The default value is [].
%   hypers           :a logical denoting when true that mach number above
%                     tsmach are hypersonic. The default value is false and
%                     those value are supersonic.
%   loop             :a scalar denoting the type of looping done to
%                     generate the DATCOM file. When loop is 1, mach and
%                     altitude are varied together. When loop is 2, mach
%                     varies while altitude is fixed. Altitude is then
%                     updated and mach numbers are cycled through again.
%                     When loop is 3, mach is fixed while altitude varies.
%                     Mach is then updated and altitudes are cycled through
%                     again. The default value is 1.
%   sref             :a scalar denoting the reference area for the case.
%                     The default value is [].
%   cbar             :a scalar denoting the longitudinal reference length.
%                     The default value is [].
%   blref            :a scalar denoting the lateral reference length. The
%                     default value is [].
%   dim              :a string denoting the specified system of units for
%                     the case. The default value is 'ft'.
%   deriv            :a string denoting the specified angle units for the
%                     case. The default value is 'deg'.
%   stmach           :a scalar value setting the upper limit of subsonic
%                     mach numbers. The default value is 0.6.
%   tsmach           :a scalar value setting the lower limit of supersonic
%                     mach numbers. The default value is 1.4.
%   save             :a logical denoting the input values for this case are
%                     used in the next case. The default value is false.
%   stype            :a scalar denoting the type of asymmetric flap for the
%                     case. The default value is [].
%   trim             :a logical denoting the reading of trim data for the
%                     case.  When trim runs are read this value is set to
%                     true. The default value is false.
%   damp             :a logical denoting the reading of dynamic derivative
%                     data for the case.  When dynamic derivative runs are
%                     read this value is set to true. The default value is
%                     false.
%   build            :a scalar denoting the reading of build data for the
%                     case.  When build runs are read this value is set to
%                     10. The default value is 1.
%   part             :a logical denoting the reading of partial data for the
%                     case.  When partial runs were written for each mach
%                     number this value is set to true. The default value
%                     is false.
%   highsym          :a logical denoting the reading of symmetric flap high
%                     lift data for the case.  When symmetric flap runs are
%                     read this value is set to true. The default value is
%                     false.
%   highasy          :a logical denoting the reading of asymmetric flap high
%                     lift data for the case.  When asymmetric flap runs are
%                     read this value is set to true. The default value is
%                     false.
%   highcon          :a logical denoting the reading of control/trim tab high
%                     lift data for the case.  When control/trim tab runs are
%                     read this value is set to true. The default value is
%                     false.
%   tjet             :a logical denoting the reading of transverse-jet
%                     control data for the case.   When transverse-jet
%                     control runs are read this value is set to true. The
%                     default value is false.
%   hypeff           :a logical denoting the reading of hypersonic flap
%                     effectiveness data for the case. When hypersonic flap
%                     effectiveness runs are read this value is set to
%                     true. The default value is false.
%   lb               :a logical denoting the reading of low aspect ratio
%                     wing or lifting body data for the case. When low
%                     aspect ratio wing or lifting body runs are read this
%                     value is set to true. The default value is false.
%   pwr              :a logical denoting the reading of power effects
%                     data for the case. When power effects runs are read
%                     this value is set to true. The default value is false.
%   grnd             :a logical denoting the reading of ground effects
%                     data for the case. When ground effects runs are read
%                     this value is set to true. The default value is false.
%   wsspn            :a scalar denoting the semi-span theoretical panel for
%                     wing.  This value is used to determine if the
%                     configuration contains a canard. The default value is 1.
%   hsspn            :a scalar denoting the semi-span theoretical panel for
%                     horizontal tail.  This value is used to determine if
%                     the configuration contains a canard. The default
%                     value is 1.
%   ndelta           :the number of control surface deflections: delta,
%                     deltal, or deltar. The default value is 0.
%   delta            :an array of control-surface streamwise deflection
%                     angles. The default value is [].
%   deltal           :an array of left lifting surface streamwise control
%                     deflection angles. The default value is [] and is
%                     defined positive for trailing-edge down.
%   deltar           :an array of right lifting surface streamwise control
%                     deflection angles. The default value is [] and is
%                     defined positive for trailing-edge down.
%   ngh              :a scalar denoting the number of ground altitudes. The
%                     default value is 0.
%   grndht           :a array denoting of ground altitudes. The default
%                     value is [].
%   config           :a logical denoting the case contains horizontal
%                     tails. The default value is false.
%
%   Static longitude and lateral stability fields available are:
%   cd               :a matrix of drag coefficients.  These coefficients
%                     are a function of alpha, mach, alt, build, ground
%                     height, and delta and defined positive for an aft
%                     acting load.
%   cl               :a matrix of lift coefficients. These coefficients are
%                     a function of alpha, mach, alt, build, ground
%                     height, and delta and defined positive for an up
%                     acting load. 
%   cm               :a matrix of pitching-moment coefficients.  These
%                     coefficients are a function of alpha, mach, alt,
%                     build, ground height, and delta and defined positive
%                     for a nose-up rotation.
%   cn               :a matrix of normal-force coefficients. These
%                     coefficients are a function of alpha, mach, alt,
%                     build, ground height, and delta and defined positive
%                     for a normal force in the +Z direction.
%   ca               :a matrix of axial-force coefficients. These
%                     coefficients are a function of alpha, mach, alt,
%                     build, ground height, and delta and defined positive
%                     for a normal force in the +X direction.
%   xcp              :a matrix of distances between moment reference center
%                     and the center of pressure divided by the
%                     longitudinal reference length. These distances are
%                     a function of alpha, mach, alt, build, ground height,
%                     and delta and defined positive for a location forward
%                     of the center of gravity.
%   cla              :a matrix of derivatives of lift coefficients w.r.t
%                     alpha. These derivatives are a function of alpha,
%                     mach, alt, build, ground height, and delta.
%   cma              :a matrix of derivatives of pitching-moment
%                     coefficients w.r.t alpha. These derivatives are a
%                     function of alpha, mach, alt, build, ground height,
%                     and delta. 
%   cyb              :a matrix of derivatives of side-force coefficients
%                     w.r.t sideslip angle. These derivatives are a
%                     function of alpha, mach, alt, build, ground height,
%                     and delta. 
%   cnb              :a matrix of derivatives of yawing-moment coefficients
%                     w.r.t sideslip angle. These derivatives are a
%                     function of alpha, mach, alt, build, ground height,
%                     and delta. 
%   clb              :a matrix of derivatives of rolling-moment
%                     coefficients w.r.t sideslip angle. These derivatives
%                     are a function of alpha, mach, alt, build, ground
%                     height, and delta. 
%   qqinf            :a matrix of ratios of dynamic pressure at the
%                     horizontal tail to the freestream value. These
%                     ratios are a function of alpha, mach, alt, build,
%                     ground height, and delta.
%   eps              :a matrix of downwash angle at horizontal tail in
%                     degrees. These angles are a function of alpha,
%                     mach, alt, build, ground height, and delta.
%   depsdalp         :a matrix of downwash angle w.r.t angle of attack.
%                     These angles are a function of alpha, mach,
%                     alt, build, ground height, and delta.
%
%   Dynamic derivative fields are:
%   clq              :a matrix of lift force derivatives due to pitch
%                     rate. These derivatives are a function of alpha,
%                     mach, alt, and build.
%   cmq              :a matrix of pitching moment derivatives due to pitch
%                     rate. These derivatives are a function of alpha,
%                     mach, alt, and build.
%   clad             :a matrix of lift force derivatives due to rate of
%                     angle of attack. These derivatives are a function
%                     of alpha, mach, alt, and build.
%   cmad             :a matrix of pitching moment derivatives due to rate
%                     of angle of attack. These derivatives are a
%                     function of alpha, mach, alt, and build.
%   clp              :a matrix of rolling moment derivatives due to roll
%                     rate. These derivatives are a function of alpha,
%                     mach, alt, and build.
%   cyp              :a matrix of lateral force derivatives due to roll
%                     rate. These derivatives are a function of alpha,
%                     mach, alt, and build.
%   cnp              :a matrix of yawing moment derivatives due to roll
%                     rate. These derivatives are a function of alpha,
%                     mach, alt, and build.
%   cnr              :a matrix of yawing moment derivatives due to yaw
%                     rate. These derivatives are a function of alpha,
%                     mach, alt, and build.
%   clr              :a matrix of rolling moment derivatives due to yaw
%                     rate. These derivatives are a function of alpha,
%                     mach, alt, and build.
%
%   High lift and control fields for symmetric flaps are:
%   dcl_sym          :a matrix of incremental lift coefficients in the
%                     linear-lift angle of attack range due to deflection
%                     of control surface. These coefficients are a function
%                     of delta, mach, and altitude.
%   dcm_sym          :a matrix of incremental pitching-moment coefficients
%                     due to deflection of control surface valid in the
%                     linear-lift angle of attack range. These coefficients
%                     are a function of delta, mach, and altitude.
%   dclmax_sym       :a matrix of incremental maximum lift coefficients.
%                     These coefficients are a function of delta, mach, and
%                     altitude.
%   dcdmin_sym       :a matrix of incremental minimum drag coefficients due
%                     to control or flap deflection. These coefficients
%                     are a function of delta, mach, and altitude.
%   clad_sym         :a matrix of the lift-curve slope of the deflected,
%                     translated surface. These coefficients are a function
%                     of delta, mach, and altitude.
%   cha_sym          :a matrix of control-surface hinge-moment derivatives
%                     due to angle of attack. These derivatives are a
%                     function of delta, mach, and altitude and when
%                     defined positive will tend to rotate the flap
%                     trailing edge down.
%   chd_sym          :a matrix of control-surface hinge-moment derivatives
%                     due to control deflection. These derivatives are a
%                     function of delta, mach, and altitude and when
%                     defined positive will tend to rotate the flap
%                     trailing edge down.
%   dcdi_sym         :a matrix of incremental induced drag coefficients due
%                     to flap deflection. These coefficients are a function
%                     of alpha, delta, mach, and altitude.
%
%   High lift and control fields available for asymmetric flaps are:
%   xsc              :a matrix of streamwise distances from wing leading
%                     edge to spoiler tip. These distances are a function
%                     of delta, mach, and altitude.
%   hsc              :a matrix of projected height of spoiler measured from
%                     normal to airfoil meanline. These distances are a
%                     function of delta, mach, and altitude.
%   ddc              :a matrix of projected height of deflector for
%                     spoiler-slot-deflector control. These distances are a
%                     function of delta, mach, and altitude.
%   dsc              :a matrix of projected height of spoiler control.
%                     These distances are a function of delta, mach, and
%                     altitude.
%   clroll           :a matrix of incremental rolling moment coefficients
%                     due to asymmetrical deflection of control surface.
%                     These coefficients are a function of delta, mach, and
%                     altitude or a function of alpha, delta, mach and
%                     altitude for differential horizontal stabilizer and
%                     defined positive when right wing is down.
%   cn_asy           :a matrix of incremental yawing moment coefficients
%                     due to asymmetrical deflection of control surface.
%                     These coefficients are a function of delta, mach, and
%                     altitude or a function of alpha, delta, mach and
%                     altitude for plain flaps and defined positive when
%                     nose is right.  
%
%   High lift and control fields available for control/trim tabs are:
%   fc_con           :a matrix of stick forces or stick force coefficients.
%                     These forces or coefficients are a function of alpha,
%                     delta, mach, and altitude.
%   fhmcoeff_free    :a matrix of flap hinge moment coefficients tab free.
%                     These coefficients are a function of alpha, delta,
%                     mach, and altitude.
%   fhmcoeff_lock    :a matrix of flap hinge moment coefficients tab
%                     locked. These coefficients are a function of alpha,
%                     delta, mach, and altitude.
%   fhmcoeff_gear    :a matrix of flap hinge moment coefficients due to
%                     gearing. These coefficients are a function of alpha,
%                     delta, mach, and altitude.
%   ttab_def         :a matrix of trim tab deflections for zero stick
%                     force. These deflections are a function of alpha,
%                     delta, mach, and altitude.
%
%   High lift and control fields available for trim are:
%   cl_utrim         :a matrix of untrimmed lift coefficients.  These
%                     coefficients are a function of alpha, mach, and
%                     altitude and defined positive for an up acting load.
%   cd_utrim         :a matrix of untrimmed drag coefficients.  These
%                     coefficients are a function of alpha, mach, and
%                     altitude and defined positive for an aft acting load.
%   cm_utrim         :a matrix of untrimmed pitching moment coefficients.
%                     These coefficients are a function of alpha, mach, and
%                     altitude and defined positive for a nose-up rotation.
%   delt_trim        :a matrix of trimmed control-surface streamwise
%                     deflection angles. These coefficients are a function
%                     of alpha, mach, and altitude.
%   dcl_trim         :a matrix of trimmed incremental lift coefficients in
%                     the linear-lift angle of attack range due to deflection
%                     of control surface. These coefficients are a function
%                     of alpha, mach, and altitude.
%   dclmax_trim      :a matrix of trimmed incremental maximum lift
%                     coefficients. These coefficients are a function of
%                     alpha, mach, and altitude.
%   dcdi_trim        :a matrix of trimmed incremental induced drag
%                     coefficients due to flap deflection. These
%                     coefficients are a function of alpha, mach, and
%                     altitude.
%   dcdmin_trim      :a matrix of trimmed incremental minimum drag
%                     coefficients due to control or flap deflection. These
%                     coefficients are a function of alpha, mach, and
%                     altitude.
%   cha_trim         :a matrix of trimmed control-surface hinge-moment
%                     derivatives due to angle of attack. These derivatives
%                     are a function of alpha, mach, and altitude.
%   chd_trim         :a matrix of trimmed control-surface hinge-moment
%                     derivatives due to control deflection. These
%                     derivatives are a function of alpha, mach, and
%                     altitude.
%   cl_tailutrim     :a matrix of untrimmed stabilizer lift coefficients.
%                     These coefficients are a function of alpha, mach, and
%                     altitude and defined positive for an up acting load.
%   cd_tailutrim     :a matrix of untrimmed stabilizer drag coefficients.
%                     These coefficients are a function of alpha, mach, and
%                     altitude and defined positive for an aft acting load.
%   cm_tailutrim     :a matrix of untrimmed stabilizer pitching moment
%                     coefficients. These coefficients are a function of
%                     alpha, mach, and altitude and defined positive for a
%                     nose-up rotation.
%   hm_tailutrim     :a matrix of untrimmed stabilizer hinge moment
%                     coefficients. These coefficients are a function of
%                     alpha, mach, and altitude and defined positive for a
%                     stabilizer rotation with leading edge up and trailing
%                     edge down.
%   aliht_tailtrim   :a matrix of stabilizer incidence required to trim.
%                     These coefficients are a function of alpha, mach, and
%                     alt.
%   cl_tailtrim      :a matrix of trimmed stabilizer lift coefficients.
%                     These coefficients are a function of alpha, mach, and
%                     altitude and defined positive for an up acting load.
%   cd_tailtrim      :a matrix of trimmed stabilizer drag coefficients.
%                     These coefficients are a function of alpha, mach, and
%                     altitude and defined positive for an aft acting load.
%   cm_tailtrim      :a matrix of trimmed stabilizer pitching moment
%                     coefficients. These coefficients are a function of
%                     alpha, mach, and altitude and defined positive for a
%                     nose-up rotation.
%   hm_tailtrim      :a matrix of trimmed stabilizer hinge moment
%                     coefficients. These coefficients are a function of
%                     alpha, mach, and altitude and defined positive for a
%                     stabilizer rotation with leading edge up and trailing
%                     edge down.
%   cl_trimi         :a matrix of lift coefficients at trim incidence.
%                     These coefficients are a function of alpha, mach, and
%                     altitude and defined positive for an up acting load.
%   cd_trimi         :a matrix of drag coefficients at trim incidence.
%                     These coefficients are a function of alpha, mach, and
%                     altitude and defined positive for an aft acting load.
%
%   Transverse jet control fields are:
%   time             :a matrix of times. These times are stored with
%                     indices of mach, alt, and alpha.
%   ctrlfrc          :a matrix of control forces. These forces are stored
%                     with indices of mach, alt, and alpha.
%   locmach          :a matrix of local Mach numbers. These Mach numbers
%                     are stored with indices of mach, alt, and alpha.
%   reynum           :a matrix of Reynolds numbers. These Reynolds numbers
%                     are stored with indices of mach, alt, and alpha.
%   locpres          :a matrix of local pressures. These pressures are
%                     stored with indices of mach, alt, and alpha.
%   dynpres          :a matrix of dynamic pressures. These pressures are
%                     stored with indices of mach, alt, and alpha.
%   blayer           :a cell array of strings containing the state of the
%                     boundary layer. These states are stored with indices
%                     of mach, alt, and alpha.
%   ctrlcoeff        :a matrix of control force coefficients. These
%                     coefficients are stored with indices of mach, alt,
%                     and alpha.
%   corrcoeff        :a matrix of corrected force coefficients. These
%                     coefficients are stored with indices of mach, alt,
%                     and alpha.
%   sonicamp         :a matrix of sonic amplification factors. These
%                     factors are stored with indices of mach, alt, and
%                     alpha.
%   ampfact          :a matrix of amplification factors. These factors are
%                     stored with indices of mach, alt, and alpha.
%   vacthr           :a matrix of vacuum thrusts. These thrusts are stored
%                     with indices of mach, alt, and alpha.
%   minpres          :a matrix of minimum pressure ratios. These ratios are
%                     stored with indices of mach, alt, and alpha.
%   minjet           :a matrix of minimum jet pressures. These pressures are
%                     stored with indices of mach, alt, and alpha.
%   jetpres          :a matrix of jet pressures. These pressures are stored
%                     with indices of mach, alt, and alpha.
%   massflow         :a matrix of mass flow rates. These rates are stored
%                     with indices of mach, alt, and alpha.
%   propelwt         :a matrix of propellant weights. These weights are
%                     stored with indices of mach, alt, and alpha.
%
%   Hypersonic fields are:
%   df_normal        :a matrix of increments in normal force per spanwise
%                     foot of control. These increments are stored with
%                     indices of alpha, delta, and mach.
%   df_axial         :a matrix of increments in axial force per spanwise
%                     foot of control. These increments are stored with
%                     indices of alpha, delta, and mach.
%   cm_normal        :a matrix of increments pitching moment in due to
%                     normal force per spanwise foot of control. These
%                     increments are stored with indices of alpha, delta,
%                     and mach.
%   cm_axial         :a matrix of increments pitching moment in due to
%                     axial force per spanwise foot of control. These
%                     increments are stored with indices of alpha, delta,
%                     and mach.
%   cp_normal        :a matrix of center of pressure locations of normal
%                     force. These locations are stored with indices of
%                     alpha, delta, and mach.
%   cp_axial         :a matrix of center of pressure locations of axial
%                     force. These locations are stored with indices of
%                     alpha, delta, and mach.
%
%   Auxiliary and partial fields available are:
%   wetarea_b        :a matrix of body wetted area. These areas are stored
%                     with indices of mach, altitude, and number of runs.
%   xcg_b            :a matrix of longitudinal locations of cg. These
%                     locations are stored with indices of mach, altitude,
%                     and number of runs (normally 1, 2 for hypers = true). 
%   zcg_b            :a matrix of vertical locations of cg. These
%                     locations are stored with indices of mach, altitude,
%                     and number of runs (normally 1, 2 for hypers = true). 
%   basearea_b       :a matrix of body base area. These areas are stored
%                     with indices of mach, altitude, and number of runs
%                     (normally 1, 2 for hypers = true). 
%   cd0_b            :a matrix of body zero lift drags. These drags are
%                     stored with indices of mach and altitude, and number
%                     of runs (normally 1, 2 for hypers = true). 
%   basedrag_b       :a matrix of body base drags. These drags are stored
%                     with indices of mach, altitude, and number of runs
%                     (normally 1, 2 for hypers = true). 
%   fricdrag_b       :a matrix of body friction drags. These drags are
%                     stored with indices of mach, altitude, and number of
%                     runs (normally 1, 2 for hypers = true). 
%   presdrag_b       :a matrix of body pressure drags. These drags are
%                     stored with indices of mach, altitude, and number of
%                     runs (normally 1, 2 for hypers = true). 
%   lemac            :a matrix of leading edge mean aerodynamic chords.
%                     These chords are stored with indices of mach and
%                     altitude.
%   sidewash         :a matrix of SIDEWASH. These values are stored with
%                     indices of mach and altitude.
%   hiv_b_w          :a matrix of IV-B(W). These values are stored with
%                     indices of alpha, mach and altitude.
%   hiv_w_h          :a matrix of IV-W(H). These values are stored with
%                     indices of alpha, mach and altitude.
%   hiv_b_h          :a matrix of IV-B(H). These values are stored with
%                     indices of alpha, mach and altitude.
%   gamma            :a matrix of GAMMA*2*PI*ALPHA*V*R. These values are
%                     stored with indices of alpha, mach and altitude.
%   gamma2pialpvr    :a matrix of GAMMA*(2*PI*ALPHA*V*R)T. These values are
%                     stored with indices of alpha, mach and altitude.
%   clpgammacl0      :a matrix of CLP(GAMMA=CL=0). These values are stored
%                     with indices of mach and altitude.
%   clpgammaclp      :a matrix of CLP(GAMMA)/CL (GAMMA=0). These values are
%                     stored with indices of mach and altitude.
%   cnptheta         :a matrix of CNP/THETA. These values are stored with
%                     indices of mach and altitude.
%   cypgamma         :a matrix of CYP/GAMMA. These values are stored with
%                     indices of mach and altitude.
%   cypcl            :a matrix of CYP/CL (CL=0). These values are stored
%                     with indices of mach and altitude.
%   clbgamma         :a matrix of CLB/GAMMA. These values are stored with
%                     indices of mach and altitude.
%   cmothetaw        :a matrix of (CMO/THETA)W. These values are stored
%                     with indices of mach and altitude.
%   cmothetah        :a matrix of (CMO/THETA)W. These values are stored
%                     with indices of mach and altitude.
%   espeff           :a matrix of (EPSOLN)EFF. These values are stored with
%                     indices of alpha, mach and altitude.
%   despdalpeff      :a matrix of D(EPSOLN)/D(ALPHA) EFF. These values are
%                     stored with indices of alpha, mach and altitude.
%   dragdiv          :a matrix of DRAG DIVERGENCE MACH NUMBER. These values
%                     are stored with indices of mach and altitude.
%   cd0mach          :a matrix of 4 mach numbers for the zero lift drag.
%                     These mach numbers are stored with indices of index,
%                     mach, and altitude.
%   cd0              :a matrix of 4 zero lift drags. These mach numbers are
%                     stored with indices of index, mach, and altitude.
%   clbclmfb_****    :a matrix of (CLB/CL)MFB where **** is either wb
%                     (wing-body) or bht (body-horz. tail). These values
%                     are stored with indices of mach and altitude.
%   cnam14_****      :a matrix of (CNA)M=1.4 where **** is either wb
%                     (wing-body) or bht (body-horz. tail). These values
%                     are stored with indices of mach and altitude.
%   area_*_**        :a matrix of areas where * is either w (wing), ht
%                     (horz. tail), vt (vert. tail), or vf (vent. fin) and
%                     ** is either tt (tot. theoretical), ti (theor.
%                     inboard), te (tot. exposed), ei (exp. inboard), or o
%                     (outboard). These areas are stored with indices of
%                     mach, altitude, and number of runs (normally 1, 2 for
%                     hypers = true). 
%   taperratio_*_**  :a matrix of taper ratios where * is either w (wing),
%                     ht (horz. tail), vt (vert. tail), or vf (vent. fin)
%                     and ** is either tt (tot. theoretical), ti (theor.
%                     inboard), te (tot. exposed), ei (exp. inboard), or o
%                     (outboard). These ratios are stored with indices of
%                     mach, altitude, and number of runs (normally 1, 2 for
%                     hypers = true).
%   aspectratio_*_** :a matrix of aspect ratios where * is either w (wing),
%                     ht (horz. tail), vt (vert. tail), or vf (vent. fin)
%                     and ** is either tt (tot. theoretical), ti (theor.
%                     inboard), te (tot. exposed), ei (exp. inboard), or o
%                     (outboard). These ratios are stored with indices of
%                     mach, altitude, and number of runs (normally 1, 2 for
%                     hypers = true).
%   qcsweep_*_**     :a matrix of quarter chord sweeps where * is either w
%                     (wing), ht (horz. tail), vt (vert. tail), or vf
%                     (vent. fin) and ** is either tt (tot. theoretical),
%                     ti (theor. inboard), te (tot. exposed), ei (exp.
%                     inboard), or o (outboard). These sweeps are stored
%                     with indices of mach, altitude, and number of runs
%                     (normally 1, 2 for hypers = true).
%   mac_*_**         :a matrix of mean aerodynamic chords where * is either
%                     w (wing), ht (horz. tail), vt (vert. tail), or vf
%                     (vent. fin) and ** is either tt (tot. theoretical),
%                     ti (theor. inboard), te (tot. exposed), ei (exp.
%                     inboard), or o (outboard). These chords are stored
%                     with indices of mach, altitude, and number of runs
%                     (normally 1, 2 for hypers = true).
%   qcmac_*_**       :a matrix of quarter chord X(MAC) where * is either w
%                     (wing), ht (horz. tail), vt (vert. tail), or vf
%                     (vent. fin) and ** is either tt (tot. theoretical),
%                     ti (theor. inboard), te (tot. exposed), ei (exp.
%                     inboard), or o (outboard). These values are stored
%                     with indices of mach, altitude, and number of runs
%                     (normally 1, 2 for hypers = true).
%   ymac_*_**        :a matrix of Y(MAC) where * is either w (wing), ht
%                     (horz. tail), vt (vert. tail), or vf (vent. fin) and
%                     ** is either tt (tot. theoretical), ti (theor.
%                     inboard), te (tot. exposed), ei (exp. inboard), or o
%                     (outboard). These values are stored with indices of
%                     mach, altitude, and number of runs (normally 1, 2 for
%                     hypers = true).
%   cd0_*_**         :a matrix of zero lift drags where * is either w
%                     (wing), ht (horz. tail), vt (vert. tail), or vf
%                     (vent. fin) and ** is either tt (tot. theoretical),
%                     ti (theor. inboard), te (tot. exposed), ei (exp.
%                     inboard), or o (outboard). These drags are stored
%                     with indices of mach, altitude, and number of runs
%                     (normally 1, 2 for hypers = true).
%   friccoeff_*_**   :a matrix of friction coefficients where * is either w
%                     (wing), ht (horz. tail), vt (vert. tail), or vf
%                     (vent. fin) and ** is either tt (tot. theoretical),
%                     ti (theor. inboard), te (tot. exposed), ei (exp.
%                     inboard), or o (outboard). These values are stored
%                     with indices of mach, altitude, and number of runs
%                     (normally 1, 2 for hypers = true). 
%   cla_b_***        :a matrix of CLA-B(***) where *** is either w (wing)
%                     or ht (stabilizer). These values are stored with
%                     indices of mach, altitude, and number of runs
%                     (normally 1, 2 for hypers = true).
%   cla_***_b        :a matrix of CLA-***(B) where *** is either w (wing)
%                     or ht (stabilizer). These values are stored with
%                     indices of mach, altitude, and number of runs
%                     (normally 1, 2 for hypers = true).
%   k_b_***          :a matrix of K-B(***) where *** is either w (wing)
%                     or ht (stabilizer). These values are stored with
%                     indices of mach, altitude, and number of runs
%                     (normally 1, 2 for hypers = true).
%   k_***_b          :a matrix of K-***(B) where *** is either w (wing)
%                     or ht (stabilizer). These values are stored with
%                     indices of mach, altitude, and number of runs
%                     (normally 1, 2 for hypers = true). 
%   xacc_b_***       :a matrix of XAC/C-B(***) where *** is either w (wing)
%                     or ht (stabilizer). These values are stored with
%                     indices of mach, altitude, and number of runs
%                     (normally 1, 2 for hypers = true).
%   cdlcl2_***       :a matrix of CDL/CL^2 where *** is either w (wing)
%                     or ht (stabilizer). These values are stored with
%                     indices of mach and altitude.
%   clbcl_***        :a matrix of CLB/CL where *** is either w (wing), ht
%                     (stabilizer), wb (wing-body) or bht (body-horz. tail).
%                     These values are stored with indices of mach and
%                     altitude.
%   fmach0_***       :a matrix of force break mach number with zero sweep
%                     where *** is either w (wing) or ht (stabilizer).
%                     These values are stored with indices of mach and
%                     altitude.
%   fmach_***        :a matrix of force break mach number with sweep where
%                     *** is either w (wing) or ht (stabilizer). These
%                     values are stored with indices of mach and altitude.
%   macha_***        :a matrix of MACH(A) where *** is either w (wing)
%                     or ht (stabilizer). These values are stored with
%                     indices of mach and altitude.
%   machb_***        :a matrix of MACH(B) where *** is either w (wing)
%                     or ht (stabilizer). These values are stored with
%                     indices of mach and altitude.
%   claa_***         :a matrix of CLA(A) where *** is either w (wing)
%                     or ht (stabilizer). These values are stored with
%                     indices of mach and altitude.
%   clab_***         :a matrix of CLA(B) where *** is either w (wing)
%                     or ht (stabilizer). These values are stored with
%                     indices of mach and altitude.
%   clbm06_***       :a matrix of (CLB/CL)M=0.6 where *** is either w (wing)
%                     or ht (stabilizer). These values are stored with
%                     indices of mach and altitude.
%   clbm14_***       :a matrix of (CLB/CL)M=1.4 where *** is either w
%                     (wing), ht (stabilizer), wb (wing-body) or bht
%                     (body-horz. tail). These values are stored with
%                     indices of mach and altitude.
%   clalpmach_***    :a matrix of 5 mach numbers for the lift curve slope
%                     where *** is either w (wing) or ht (stabilizer).
%                     These mach numbers are stored with indices of index,
%                     mach, and altitude.
%   clalp_***        :a matrix of 5 lift curve slope values where *** is
%                     either w (wing) or ht (stabilizer). These values are
%                     stored with indices of index, mach, and altitude.
