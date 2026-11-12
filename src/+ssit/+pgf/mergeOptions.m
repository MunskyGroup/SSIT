function opts = mergeOptions(src, opts)
% opts = ssit.pgf.mergeOptions(src, opts) completes PGF options: empty opts
% are taken from src.pgfOptions (if present), and missing fields get the
% values of ssit.pgf.defaultOptions().
defaults = ssit.pgf.defaultOptions();
if isempty(opts)
    if isPropOrField(src, 'pgfOptions') && ~isempty(src.pgfOptions)
        opts = src.pgfOptions;
    else
        opts = defaults;
    end
end
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i})
        opts.(names{i}) = defaults.(names{i});
    end
end
end

function tf = isPropOrField(obj, name)
tf = (isstruct(obj) && isfield(obj, name)) || (isobject(obj) && isprop(obj, name));
end
