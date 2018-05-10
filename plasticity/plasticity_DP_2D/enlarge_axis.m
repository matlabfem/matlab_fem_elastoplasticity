function enlarge_axis(alpha1,alpha2)
    ax=axis;
    if strcmp(get(gca,'XScale'),'log')     
        ax(1)=log10(ax(1));
        ax(2)=log10(ax(2));    
    end
    if strcmp(get(gca,'YScale'),'log')     
        ax(3)=log10(ax(3));
        ax(4)=log10(ax(4));    
    end
    
    ax_new=ax*[[1+alpha2 -alpha2; -alpha2 1+alpha2], zeros(2); zeros(2), [1+alpha1 -alpha1; -alpha1 1+alpha1]];
    
    if strcmp(get(gca,'XScale'),'log')
        ax_new(1)=10^ax_new(1);
        ax_new(2)=10^ax_new(2);       
    end
    if strcmp(get(gca,'YScale'),'log')
        ax_new(3)=10^ax_new(3);
        ax_new(4)=10^ax_new(4);       
    end
      
    axis(ax_new);
end

