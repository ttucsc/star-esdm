function [ipaddr, ipdigits] = ic_get_ip_address()
%   returns the IP address of machine as a char vector, and the digits of the ip address as a numerical vector.
%   ignores IP addresses beginning with 192, 10 and 127.
%   returns empty string and vector of 4 nan's if it cannot get the IP address.

    try       
        if (strcmpi(computer(),'MACI64'))
            [~,ipaddr]=system('ifconfig | grep -w inet | awk ''{print $2}'' | grep -v "127." | grep -v "192." | grep -v " 10."');
            [ipaddr, ipdigits] = get_digits(ipaddr);


        elseif (strcmpi(computer(),'GLNXA64'))
            [~,ipaddr]=system('ifconfig | grep "inet" | awk ''{print $2}'' | grep -v ":192." | grep -v ":10." | grep -v ":127."  | awk -F ":" ''{print $2}''');
            [ipaddr, ipdigits] = get_digits(ipaddr);

        elseif (strcmpi(computer(), 'PCWIN64'))
            ipaddr=get_win_ip(ipcout);
            [ipaddr, ipdigits] = get_digits(ipaddr);
        else
            ipaddr='';
            ipdigits=nan(1,4);
        end
    catch
       ipaddr='';
       ipdigits=nan(1,4);
    end
end

function [ipaddr,ipdigs] = get_digits(ipaddr)
    
    if (isempty(ipaddr) || strcmp(ipaddr,newline) || strcmp(ipaddr, char(13)) || strcmp(ipaddr,char([13,10])))
        ipaddr='';
        ipdigs=nan(1,4);
    else
        ipaddr=strtrim(ipaddr);
        z=split(ipaddr,'.');
        ipdigs=nan(1,length(z));
        for i=1:length(z)
            ipdigs(i)=str2double(z{i});
        end
    end
end

function ipaddr = get_win_ip()

    [~,ipcout]=system('ipconfig all');
    z=split(ipcout,newline);
    z=z(contains(z,'IP address'));

    z=z(~contains(' 192.'));
    z=z(~contains(' 10.'));
    z=z(~contains(' 127.'));
    ipaddr=strtrim(z{1});
end