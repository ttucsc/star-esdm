function J = cmap_CSC(m,top_or_bottom)
%cmap_CSC(m,top_or_bottom)  
%   returns a variation of JET colormap which ranges from blue (low) to red (high), but is white in the middle.
%   Inputs:
%       m               integer.  # of rows in colormap
%       top_or_bottom   optional.  logical. Leave off to return full colormap. 
%                                           true:  returns top half only, starting with white, going to red
%                                           false: returns bottom half only, starting with blue, going to white.
%                           note:  m is size of returned map, even if top_or_bottom is present.
%                           note 2:  top_or_bottom may also be specified as "top" or "bottom"
%
%   Here's the description of JET:
%JET    Variant of HSV
%   JET(M) returns an M-by-3 matrix containing the jet colormap, a variant
%   of HSV(M). The colors begin with dark blue, range through shades of
%   blue, cyan, green, yellow and red, and end with dark red. JET, by
%   itself, is the same length as the current figure's colormap. If no
%   figure exists, MATLAB uses the length of the default colormap.
%
%   See also PARULA, HSV, HOT, PINK, FLAG, COLORMAP, RGBPLOT.

%   Copyright 1984-2015 The MathWorks, Inc.

    if nargin < 1
       f = get(groot,'CurrentFigure');
       if isempty(f)
          m = size(get(groot,'DefaultFigureColormap'),1);
       else
          m = size(f.Colormap,1);
       end
    end
    
    if (exist('top_or_bottom','var')), m = 2*m+1; end

    J = zeros(m,3);    
    
    n = floor((m-1)/4);         % n is # of 1's in each color band.
    if (mod(m,2))
        n = n + ~mod(n,2);      % n odd if m is odd
    else
        n = n + mod(n,2);       % m even if n is even
    end
    
    nz = (m-n)/2;               % # of non-1's at each end.
    
    u = [(0:(nz-1))/nz, ones(1,n), ((nz-1):-1:0)/nz];   % first 1 is at (nz+1)
    
    mid = ceil((m+1)/2);    % where we want the first red 1.
    
    rix1 = mid - nz;
    rix = rix1:m;
    J(rix,1) = u(1:length(rix));
    if (J(end,1)==0), J(:,1) = circshift(J(:,1),1,1); end
    if (mod(m,2))
        J(:,3) = flipud(J(:,1));
    else
        shft = ceil(n/2);
        bix = shft:m;
        J(1:length(bix),3) = u(bix);
    end
    
%   gix1 = max(floor(n/10),floor(n/3));
    if (m < 32)
        gix1 = ceil(n/2);
    elseif (m < 64)
        gix1 = ceil(n/2.5);
    else
        gix1 = ceil(n/3);
    end
    gix2 = mid - gix1+1;
    if(2*gix2-1<=m)
        endix = 2;
    else
        endix = 3;
    end
    
    gix = [1:gix2,gix2:-1:endix];
    J(gix1:gix1+length(gix)-1,2) = u(gix);    

    J(1:mid,3) = J(1:mid,3).^.2;
    J(mid:end,1) = J(mid:end,1).^.3;
    J(mid+1,2) = (J(mid+1,2)+J(mid+1,3))/2;
    J(mid+1,3) = (J(mid+1,3)+J(mid+2,3))/2;
    
    if (m < 32)
%       J(:,2) = J(:,2).^1.5;
        J(:,2:3) = J(:,2:3).^.75;
    elseif (m >= 48)
        J(:,2:3) = J(:,2:3).^.75;
        J(:,2) = circshift(J(:,2),1,1);
    end
%     J(:,1) = J(:,1).^.8;
%     J(:,2) = J(:,2).^.6;
%     J(:,3) = J(:,3).^.7;

    if (exist('top_or_bottom','var'))
        if (ischar_s(top_or_bottom))
            if (strncmpi(top_or_bottom,'top',3))
                top_or_bottom = true;
            else
                top_or_bottom = false;
            end
        end
        ix = ceil(m/2);
        if (top_or_bottom)
            J = J(ix:end,:);
        else
            J = J(1:ix,:);
        end
    end            
end

