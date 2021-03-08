for model = 1:2
    if model == 1
        Struct = SalISIstruct;
        mName = 'Saline';
    else
        Struct = CfaISIstruct;
        mName = 'CFA';
    end
    Baseline = Struct(1).Stacks(1).ConditionMap;
    for a = 2: length(Struct)
        Baseline = [Baseline;Struct(a).Stacks(1).ConditionMap];
    end
    Ind =[];
    for e = 1:length(Struct)
        if contains(Struct(e).name, 'Mech_Control', 'IgnoreCase', true) == true
            Ind = [Ind; e];
        end
    end
   
    Evoked = Struct(Ind(1)).Stacks(2).ConditionMap;
    for a  = 2:length(Ind)
        Evoked = [Evoked; Struct(Ind(a)).Stacks(2).ConditionMap];
    end
    Baseline = log10(Baseline);
    Evoked = log10(Evoked);
    for figs  = 1:2
        if figs == 1
            responseName = 'Baseline';
            response = Baseline;
        else
            responseName = 'Evoked';
            response = Evoked;
        end
        figure('Name', [mName, ' ', responseName, ' ISI Return Map'], 'Color', 'White');
        Hist = histogram2(response(:,1),response(:,2), 125, 'Normalization', 'probability');
        title([mName, ' ', responseName, ' ISI Return Map'])
        xlabel('ISI_n _(_m_s_e_c_s_)');
        ylabel('ISI_n_+_1 _(_m_s_e_c_s_)')
        xlim([-3, 0]);
        ylim([-3, 0]);
        zlim([0, 0.008])
        zlabel('Probability')
        fig = gcf;
        ax = gca;
        ax.FontSize = 20;
        xticks(-3:0)
        ax.XTickLabel = (10.^cellfun(@str2double,ax.XTickLabel)) * 1e3;
        ax.XTickLabelRotation = -45;
        yticks(-3:0);
        ax.YTickLabel = (10.^cellfun(@str2double,ax.YTickLabel)) * 1e3;
        ax.YTickLabelRotation = 45;
    end
end
