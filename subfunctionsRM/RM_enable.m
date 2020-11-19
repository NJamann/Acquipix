function RM_enable(handles)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	%Enable,'off'
	
	cellNames = fieldnames(handles);
	for intPtr=1:numel(cellNames)
		if ~isempty(strfind(cellNames{intPtr},'ptrButton')) ||...
				~isempty(strfind(cellNames{intPtr},'ptrList')) ||...
				~isempty(strfind(cellNames{intPtr},'ptrEdit'))
			if strcmpi(get(handles.(cellNames{intPtr}),'UserData'),'lock')
				set(handles.(cellNames{intPtr}),'UserData','unlock');
			end
		end
	end
end

