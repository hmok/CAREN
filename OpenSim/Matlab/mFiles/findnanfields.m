function fnames = findnanfields(s)
   %s: a scalar structure
   %fnames: column cell array of field names that are nan. recurses through structures
   fnames = {};
   for fn = fieldnames(s)'
      fieldcontent = s.(fn{1});
      if isstruct(fieldcontent)
           subfields = findnanfields(fieldcontent);
           if ~isempty(subfields)
              fnames = [fnames; strcat(fn{1}, '.', subfields)];
           end
      elseif isnan(fieldcontent)
           fnames = [fnames; fn{1}];
      end
   end
end