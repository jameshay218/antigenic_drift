function plot_SIR()
        Y_sl = get_SIR('../outputs/single_fixed_low/.002/');
        Y_sh = get_SIR('../outputs/single_fixed_high/.001/');
        Y_ada = get_SIR('../outputs/single_adaptive/.006/');
        
        Y_slsort = sort(Y_sl','descend');
        total = round(length(Y_slsort(:,1))*0.95);
        leng = length(Y_slsort(1,:));
        Y_sllb = Y_slsort(total,:)';
        Y_slsort = sort(Y_sl');
        Y_slub = Y_slsort(total,:)';
        figure;
        hold;
        plot(mean(Y_sl,2));
        h = fill([1:leng leng:-1:1], [Y_sllb' fliplr(Y_slub')],'r');  
        set(h,'facealpha',.5)
        %hold;
        plot(Y_sh);
        plot(Y_ada);

end

