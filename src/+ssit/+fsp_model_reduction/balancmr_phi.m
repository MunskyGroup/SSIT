function [Gred,info,SLBIG,SRBIG] = balancmr_phi(G,varargin)  %BEM - Record Transformation Matrices

%BALANCMR  balanced truncation model reduction for continuous/discrete &
% stable/unstable plant.
%
%[GRED,REDINFO] = BALANCMR(G,ORDER,KEY1,VALUE1,KEY2,VALUE2,...) performs 
%    balanced truncation model reduction on LTI system G such that the 
%    sigma plot of the error (GRED - G) < 2*(the sum of the discarded 
%    states' Hankel singular values). For unstable G the algorithm works 
%    by first splitting G into a sum of stable and antistable part, 
%    reduces the stable part, then add the unstable part back for 
%    the final output G. The unstable part is kept as default.
%
% Inputs:
%    G  - LTI system to be reduced
%
% Optional inputs:
%    ORDER - an integer array with desired order of reduced model. 
%		 A batch run of [m:n] can be specified for a pack of reduced 
% 		 order models to be generated. 
%
%   KEY       |VALUE        | MEANING
%  -------------------------------------------------------------------------
%  'MaxError' |real no. or  | Reduce to achieve H-infinity error
%             |vector       | If present, 'MaxError' overrides ORDER input.
%  ------------------------------------------------------------------------- 
%  'Weights'  |{Wout,Win}   | Optional 1x2 cell array  of LTI weights Wout
%             | cell        | (output)and Win (input); default is both
%             |             | identity. Weights have to be invertible.
%  -------------------------------------------------------------------------
%  'Display'  |'off' or 'on'| Display HANKELSV plots (default 'off')
%
%  Outputs:
%    GRED    - LTI reduced order system
%    REDINFO - a struct of REDINFO.ErrorBound, REDINFO.StabSV and
%              REDINFO.UnstabSV
%
% See also	REDUCE, SCHURMR, BSTMR, NCFMR, HANKELMR. 

% R. Y. Chiang & M. G. Safonov 4/15/02
%   Copyright 1988-2008 The MathWorks, Inc.
% All Rights Reserved. 

if nargin == 1
   Type = 3;
end

if hasInternalDelay(G)
    error(message('Robust:analysis:balancmr1'));
end

try 
    G = ss(G);
catch
    error(message('Robust:analysis:balancmr2'));
end

if ndims(G) > 2
    error(message('Robust:analysis:balancmr3'));
end   

D = get(G,'d');

% Handle discrete case:

paaz = 2^-4*sqrt(eps);
Ts0 = getTs(G);
Ts = abs(Ts0); % take care of the inherited case (Ts = -1)

if Ts,
    a = get(G,'a');
    nx = size(a,1);
    if (max(abs(eig(a)+1))>paaz)  % if no poles near z=-1
        paaz = 0;
    end
    G = bilin(G,-1,'S_Tust',[Ts, 1-paaz]);
end

Gorig = G;
[a,b,c,d] = ssdata(G);

% Find no. of stable poles + jw-axis poles
lamda = eig(a);
ra = length(lamda);
mjws = length(find(real(lamda) < sqrt(eps)));
no_unstable = ra-mjws;

% Find no. of jw-axis poles
no_jw_pole = 0;
Gjw = ss([],[],[],zeros(size(d)));
Gu =  ss([],[],[],zeros(size(d)));
no_jw_pole = length((find(abs(real(pole(G))) < sqrt(eps))));

if no_jw_pole ~= 0
    gsjw = G;
    if no_unstable ~= 0
        [Gu,gsjw] = modreal(G, no_unstable ,'real');
        [aa,bb,cc,dd] = ssdata(Gu);
        Gu = ss(aa,bb,cc,zeros(size(dd)));
    end
    [Gjw,Gs] = modreal(gsjw, no_jw_pole, 'real');
    G = Gu+Gs;
    [aa,bb,cc,dd] = ssdata(Gjw);
    Gjw = ss(aa,bb,cc,zeros(size(dd)));
end

% Check stability of G
 
if no_unstable ~=0
    [G,P,m] = stabproj(G);
else
    P = ss([],[],[],zeros(size(get(G,'d'))));
    m = size(G,'o');
end

% Processing input arguments

vvaagg = varargin;

% Defaults:

Type = 3;   % reduced style default (when order or error is not specified)
Typeoriginal = 1;
redval = {[],[]};
dispval = 'off';
wtflag = 0;
displayflag = 0;

% Checking the "order" case (2nd argument is not char string)

if ~isempty(vvaagg)
   if isnumeric(vvaagg{1,1}) 
      redval = vvaagg{1,1}-no_unstable-no_jw_pole; % Reserve unstable part in reduced model
      Type = 1;
      vvaagg = varargin(2:end);
   end 
end

if ~isempty(vvaagg)
   if mod(length(vvaagg),2) 
     error(message('Robust:analysis:balancmr4'));
   end
end

% Loop through varargin:

if ~isempty(vvaagg)
for i = 1:2:length(vvaagg)-1
    curarg0 = vvaagg{1,i};
    curarg = lower(curarg0(1:3));
    if ~isstr(curarg)
        error(message('Robust:analysis:balancmr5'));
    end
    
    switch curarg
      
    % Checking 'order'
    case 'ord'
         Type = 1;
         Typeoriginal = Type;
         redval = vvaagg{1,i+1}-no_unstable-no_jw_pole; % Reserve unstable part in reduced model
    
    % Checking 'MaxError'   
    case 'max' 
        
        Type = 2;
        Typeoriginal = Type;
        redval = vvaagg{1,i+1};   % overwrites the "order" value

    % Checking 'Weight'
    case 'wei'
        
        wtflag = 1;
        wtval  = vvaagg{1,i+1};
        sys1 = ss(wtval{1,1});
        sys2 = ss(wtval{1,2});
        wtTs1 = getTs(sys1);
        wtTs2 = getTs(sys2);
        if wtTs1 == 0 || wtTs2 == 0   % take care Tf = 1; for cont./discrete cases
           sys1 = setTs(sys1,wtTs1+wtTs2);
           sys2 = setTs(sys2,wtTs1+wtTs2);
        end
        
        wtval = {sys1,sys2};
        
        if (wtTs1 || wtTs2) 
           wta1 = get(wtval{1,1},'a');
           wta2 = get(wtval{1,2},'a');
           if (max(abs(eig(wta1)+1))>paaz) && (max(abs(eig(wta2)+1))>paaz) % if no poles near z=-1
               paaz = 0;
           end
           if isempty(wta1)
              wt1 = sys1;
           else
              wt1 = bilin(wtval{1,1},-1,'S_Tust',[wtTs1, 1-paaz]); 
           end
           if isempty(wta2)
              wt2 = sys2; 
           else
              wt2 = bilin(wtval{1,2},-1,'S_Tust',[wtTs2, 1-paaz]);
           end
           wtval = {wt1,wt2};
        end
        if length(find(real(pole(wtval{1,1}))> 0)) || length(find(real(pole(wtval{1,2}))>0)) ...
         || length(find(real(tzero(wtval{1,1}))>0)) || length(find(real(tzero(wtval{1,2}))>0))
            error(message('Robust:analysis:balancmr6'));
        end        
        try
           G = wtval{1,1}'\G/wtval{1,2}';
        catch
           error(message('Robust:analysis:balancmr7'));
        end
        [G,Gus] = stabproj(G);     
        
    % Checking 'Display'
    case 'dis' 
        dispval = vvaagg{1,i+1};
        if strcmp(dispval(1:2),'on')
           Type = 3;
           displayflag = 1;
    %       redval = [];
        end
    otherwise        
        error(message('Robust:analysis:balancmr8',curarg0));
    end
end % for loop
end % if vvaagg is not empty

% ----- Start execution:

if ~iscell(redval)
   redval = {redval};
end

if displayflag == 1 || Type == 3
   hankelsv(Gorig)
   if isempty(redval{1,1})
      ord = input(getString(message('Robust:analysis:balancmr9',num2str(no_unstable+no_jw_pole))));        
      if isempty(ord)
           ord = length(ssdata(Gorig));
           disp(getString(message('Robust:analysis:balancmr10')));           
      end
      close(gcf)
      redval{1,1} = ord-no_unstable-no_jw_pole;
   end
   Type = Typeoriginal;
end

stabexe = length(redval{1,1});
Ebound = zeros(stabexe,1);

Gred = [];
for i = 1:stabexe   
   if redval{1,1}(i) < 0
      no_ns = no_unstable+no_jw_pole;
      warning(message('Robust:analysis:balancmr11',redval{1,1}(i)+no_ns,no_ns,no_ns));
      redval{1,1}(i) = 0;
   end
   [gr,aug,hhsv,SLBIG,SRBIG] = balsq(G,Type,redval{1,1}(i));  %BEM - Record Transformation Matrices
    if wtflag % exists weighting
        ggrr = wtval{1,1}'*(gr+Gus)*wtval{1,2}';
        [gr,gtemp] = stabproj(ggrr);
        gr = gr + get(gtemp,'d');
    end 
     gr = gr+P+Gjw; % add back the unstable & integral part of G
     if Ts, % discrete case
        gr = bilin(gr,1,'S_Tust',[Ts,1-paaz]);
        gr = setTs(gr,Ts0);
     end
     set(gr,'d',D)    % put back original D matrix
     Gred = subsasgn(Gred,substruct('()',{':' ':' i}),gr);  % Gred(:,:,i) = gr; 
     Ebound(i) = aug(2);
     info.ErrorBound = Ebound;
     info.StabSV = [inf*ones(no_jw_pole,1);hhsv];
%     keyboard;
%     info.UnstabSV = hksv(P);
end
  
% ---------- End of BALANCMR.M % 10/26/02
            
    
