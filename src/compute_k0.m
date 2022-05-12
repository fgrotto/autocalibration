function [K0] = compute_k0(S, directory)
    image = imread([directory S{1,1}.name_view_i]);
    info = imfinfo([directory S{1,1}.name_view_i]);
    [H, W, C]= size(image);
    SensorW=35;
    Fmm=info.DigitalCamera.FocalLengthIn35mmFilm;
    fp=(Fmm*W)/SensorW;
    u_0=W/2;
    v_0=H/2; 
    K0=[fp 0 u_0; 0 fp v_0; 0 0 1];
end