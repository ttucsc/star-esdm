function fullpath = resolvePath(filename, verify_exists)
%   returns full path of filename
  file=java.io.File(filename);
  if file.isAbsolute()
      fullpath = filename;
  else
      fullpath = fullfile(pwd(),filename);
  end
  if (exist('verify_exists','var') && verify_exists)
      if (isfile(fullpath))
          return
      else
          error('resolvePath:CannotResolve', 'Does not exist or failed to resolve absolute path for %s.', filename);
      end
  end
end
