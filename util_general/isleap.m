function yn = isleap(yr)
    %   returns true if yr is a leap year, false if not.
    yn=zeros(size(yr));
    for i=1:length(yr)
        yr(i) = floor(yr(i));
        yn(i) = mod(yr(i),4)==0 && (mod(yr(i),400)==0 || mod(yr(i),100)~=0);
    end
end
    

