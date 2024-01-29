function [Coords]=brlinexya(Sx,Sy,Ex,Ey)
% function [Coords]=brlinexya(Sx,Sy,Ex,Ey)
% Bresenham line algorithm.
% Sx, Sy, Ex, Ey - desired endpoints
% Coords - nx2 ordered list of x,y coords.
% Author: Andrew Diamond;
%
%	if(length(M) == 0)
%		M = zeros(max([Sx,Sy]),max([Ex,Ey]));
%	end
	Dx = Ex - Sx;
	Dy = Ey - Sy;
%	Coords = [];
	CoordsX = zeros(2 .* ceil(abs(Dx)+abs(Dy)),1);
	CoordsY = zeros(2 .* ceil(abs(Dx)+abs(Dy)),1);
    iCoords=0;
	if(abs(Dy) <= abs(Dx))
		if(Ey >= Sy)
			if(Ex >= Sx)
				D = 2*Dy - Dx;
				IncH = 2*Dy;
				IncD = 2*(Dy - Dx);
				X = Sx;
				Y = Sy;
%				M(Y,X) = Value;
				% Coords = [Sx,Sy];
                iCoords = iCoords + 1;
                CoordsX(iCoords) = Sx;
                CoordsY(iCoords) = Sy;
				while(X < Ex)
					if(D <= 0)
						D = D + IncH;
						X = X + 1;
					else
						D = D + IncD;
						X = X + 1;
						Y = Y + 1;
					end
%					M(Y,X) = Value;
                    iCoords = iCoords + 1;
                    CoordsX(iCoords) = X;
                    CoordsY(iCoords) = Y;
					% Coords = [Coords; [X,Y]];
				end
			else % Ex < Sx
				D = -2*Dy - Dx;
				IncH = -2*Dy;
				IncD = 2*(-Dy - Dx);
				X = Sx;
				Y = Sy;
%				M(Y,X) = Value;
				% Coords = [Sx,Sy];
                iCoords = iCoords + 1;
                CoordsX(iCoords) = Sx;
                CoordsY(iCoords) = Sy;
				while(X > Ex)
					if(D >= 0)
						D = D + IncH;
						X = X - 1;
					else
						D = D + IncD;
						X = X - 1;
						Y = Y + 1;
					end
%					M(Y,X) = Value;
                    iCoords = iCoords + 1;
                    CoordsX(iCoords) = X;
                    CoordsY(iCoords) = Y;
%					Coords = [Coords; [X,Y]];
				end
			end
		else % Ey < Sy
			if(Ex >= Sx)
				D = 2*Dy + Dx;
				IncH = 2*Dy;
				IncD = 2*(Dy + Dx);
				X = Sx;
				Y = Sy;
%				M(Y,X) = Value;
				% Coords = [Sx,Sy];
                iCoords = iCoords + 1;
                CoordsX(iCoords) = Sx;
                CoordsY(iCoords) = Sy;
				while(X < Ex)
					if(D >= 0)
						D = D + IncH;
						X = X + 1;
					else
						D = D + IncD;
						X = X + 1;
						Y = Y - 1;
					end
%					M(Y,X) = Value;
                    iCoords = iCoords + 1;
                    CoordsX(iCoords) = X;
                    CoordsY(iCoords) = Y;
					% Coords = [Coords; [X,Y]];
				end
			else % Ex < Sx
				D = -2*Dy + Dx;
				IncH = -2*Dy;
				IncD = 2*(-Dy + Dx);
				X = Sx;
				Y = Sy;
%				M(Y,X) = Value;
				% Coords = [Sx,Sy];
                iCoords = iCoords + 1;
                CoordsX(iCoords) = Sx;
                CoordsY(iCoords) = Sy;
				while(X > Ex)
					if(D <= 0)
						D = D + IncH;
						X = X - 1;
					else
						D = D + IncD;
						X = X - 1;
						Y = Y - 1;
					end
%					M(Y,X) = Value;
                    iCoords = iCoords + 1;
                    CoordsX(iCoords) = X;
                    CoordsY(iCoords) = Y;
%					Coords = [Coords; [X,Y]];
				end
			end
		end
	else % abs(Dy) > abs(Dx) 
		Tmp = Ex;
		Ex = Ey;
		Ey = Tmp;
		Tmp = Sx;
		Sx = Sy;
		Sy = Tmp;
		Dx = Ex - Sx;
		Dy = Ey - Sy;
		if(Ey >= Sy)
			if(Ex >= Sx)
				D = 2*Dy - Dx;
				IncH = 2*Dy;
				IncD = 2*(Dy - Dx);
				X = Sx;
				Y = Sy;
%				M(X,Y) = Value;
				% Coords = [Sx,Sy];
                iCoords = iCoords + 1;
                CoordsX(iCoords) = Sy;
                CoordsY(iCoords) = Sx;
				while(X < Ex)
					if(D <= 0)
						D = D + IncH;
						X = X + 1;
					else
						D = D + IncD;
						X = X + 1;
						Y = Y + 1;
					end
%					M(X,Y) = Value;
                    iCoords = iCoords + 1;
                    CoordsX(iCoords) = Y;
                    CoordsY(iCoords) = X;
%					Coords = [Coords; [Y,X]];
				end
			else % Ex < Sx
				D = -2*Dy - Dx;
				IncH = -2*Dy;
				IncD = 2*(-Dy - Dx);
				X = Sx;
				Y = Sy;
%				M(X,Y) = Value;
				% Coords = [Sx,Sy];
                iCoords = iCoords + 1;
                CoordsX(iCoords) = Sy;
                CoordsY(iCoords) = Sx;
				while(X > Ex)
					if(D >= 0)
						D = D + IncH;
						X = X - 1;
					else
						D = D + IncD;
						X = X - 1;
						Y = Y + 1;
					end
%					M(X,Y) = Value;
                    iCoords = iCoords + 1;
                    CoordsX(iCoords) = Y;
                    CoordsY(iCoords) = X;
%					Coords = [Coords; [Y,X]];
				end
			end
		else % Ey < Sy
			if(Ex >= Sx)
				D = 2*Dy + Dx;
				IncH = 2*Dy;
				IncD = 2*(Dy + Dx);
				X = Sx;
				Y = Sy;
%				M(X,Y) = Value;
				% Coords = [Sx,Sy];
                iCoords = iCoords + 1;
                CoordsX(iCoords) = Sy;
                CoordsY(iCoords) = Sx;
				while(X < Ex)
					if(D >= 0)
						D = D + IncH;
						X = X + 1;
					else
						D = D + IncD;
						X = X + 1;
						Y = Y - 1;
					end
%					M(X,Y) = Value;
                    iCoords = iCoords + 1;
                    CoordsX(iCoords) = Y;
                    CoordsY(iCoords) = X;
%					Coords = [Coords; [Y,X]];
				end
			else % Ex < Sx
				D = -2*Dy + Dx;
				IncH = -2*Dy;
				IncD = 2*(-Dy + Dx);
				X = Sx;
				Y = Sy;
%				M(X,Y) = Value;
				% Coords = [Sx,Sy];
                iCoords = iCoords + 1;
                CoordsX(iCoords) = Sy;
                CoordsY(iCoords) = Sx;
				while(X > Ex)
					if(D <= 0)
						D = D + IncH;
						X = X - 1;
					else
						D = D + IncD;
						X = X - 1;
						Y = Y - 1;
					end
%					M(X,Y) = Value;
                    iCoords = iCoords + 1;
                    CoordsX(iCoords) = Y;
                    CoordsY(iCoords) = X;
%					Coords = [Coords; [Y,X]];
				end
			end
		end
	end
Coords = [CoordsX(1:iCoords),CoordsY(1:iCoords)];
