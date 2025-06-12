function save_param(ifquantization, path_save, mesh_name, X, T, UV, TUV, sing, E2V_hardedge)

if ~exist(path_save, 'dir')
    mkdir(path_save);
end

id_sing_p = sing > 1/8;
id_sing_m = sing <-1/8;

writeObj([path_save, mesh_name, '_pos.obj'], X(id_sing_p,:), []);
writeObj([path_save, mesh_name, '_neg.obj'], X(id_sing_m,:), []);
writeObj([path_save, mesh_name, '_param.obj'], X, T, UV, TUV, [], [], E2V_hardedge);

if ifquantization
    if  exist('./QuantizationYoann/build/Quantization', 'file') ~= 0
        status = system(['./QuantizationYoann/build/Quantization -s a -sa ', num2str(1), ' -r -o ', path_save, mesh_name, '_quantiz.obj ', path_save, mesh_name, '_param.obj']);
        if status ~= 0
            warning('Quantization: Yoann failed me :(');
        end
    else
        error('Must compile the Quantization program. Go to folder QuantizationYoann/');
    end
end