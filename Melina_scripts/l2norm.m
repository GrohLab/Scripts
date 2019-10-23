
function norm_acc = l2norm(aux_input_data)

norm_acc=sqrt(sum(aux_input_data.*aux_input_data));
norm_acc=norm_acc-mean(norm_acc);

