function counts = import_cellcounter_xmls(xml_path)

xdoc = xmlread([xml_path]);
st = xml2struct(xdoc);
    
% keyboard
% Count instances of each marker class
for k = 1:numel(st.CellCounter_Marker_File.Marker_Data.Marker_Type)
    if (k<= numel(st.CellCounter_Marker_File.Marker_Data.Marker_Type)) && ...
            isfield(st.CellCounter_Marker_File.Marker_Data.Marker_Type{k},'Marker');
        counts(k) = numel(st.CellCounter_Marker_File.Marker_Data.Marker_Type{k}.Marker);
    end
end
