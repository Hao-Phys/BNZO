parameters();
exchange_tensors();
bond_info();
return_value = local_to_global();

if return_value == 0
    return
end

qx = 0.0;
qy = 0.0;
qz = 0.0;

omega = 0:0.01:1;

sc_int = intensity_v2(omega, qx, qy, qz);
