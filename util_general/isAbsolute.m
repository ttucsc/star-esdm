function tf = isAbsolute(filename)
%   returns true if filename is an absolute path.
  file=java.io.File(filename);
  tf = file.isAbsolute();
end