function DF(problem,NT)
    if problem == "FDA2"
       
        frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1.75*1.75)*10+0.75));
        hold on
        plot(frontx,fronty,'Linewidth',1);
        frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*1/NT)))*(1+(0.75+0.7*sin(0.5*pi*1/NT)))*10+(0.75+0.7*sin(0.5*pi*1/NT)) ));
        plot(frontx,fronty,'Linewidth',1);
        frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*2/NT)))*(1+(0.75+0.7*sin(0.5*pi*2/NT)))*10+(0.75+0.7*sin(0.5*pi*2/NT)) ));
        plot(frontx,fronty,'Linewidth',1);
        frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*3/NT)))*(1+(0.75+0.7*sin(0.5*pi*3/NT)))*10+(0.75+0.7*sin(0.5*pi*3/NT)) ));
        plot(frontx,fronty,'Linewidth',1);
        frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*4/NT)))*(1+(0.75+0.7*sin(0.5*pi*4/NT)))*10+(0.75+0.7*sin(0.5*pi*4/NT)) ));
        plot(frontx,fronty,'Linewidth',1);
        frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*5/NT)))*(1+(0.75+0.7*sin(0.5*pi*5/NT)))*10+(0.75+0.7*sin(0.5*pi*5/NT)) ));
        plot(frontx,fronty,'Linewidth',1);
        frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*11/NT)))*(1+(0.75+0.7*sin(0.5*pi*11/NT)))*10+(0.75+0.7*sin(0.5*pi*11/NT)) ));
        plot(frontx,fronty,'Linewidth',1);
        frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*12/NT)))*(1+(0.75+0.7*sin(0.5*pi*12/NT)))*10+(0.75+0.7*sin(0.5*pi*12/NT)) ));
        plot(frontx,fronty,'Linewidth',1);
        frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*13/NT)))*(1+(0.75+0.7*sin(0.5*pi*13/NT)))*10+(0.75+0.7*sin(0.5*pi*13/NT)) ));
        plot(frontx,fronty,'Linewidth',1);
        frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*14/NT)))*(1+(0.75+0.7*sin(0.5*pi*14/NT)))*10+(0.75+0.7*sin(0.5*pi*14/NT)) ));
        plot(frontx,fronty,'Linewidth',1);
        frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*15/NT)))*(1+(0.75+0.7*sin(0.5*pi*15/NT)))*10+(0.75+0.7*sin(0.5*pi*15/NT)) ));
        plot(frontx,fronty,'Linewidth',1);
        frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*21/NT)))*(1+(0.75+0.7*sin(0.5*pi*21/NT)))*10+(0.75+0.7*sin(0.5*pi*21/NT)) ));
        plot(frontx,fronty,'Linewidth',1);
        frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*22/NT)))*(1+(0.75+0.7*sin(0.5*pi*22/NT)))*10+(0.75+0.7*sin(0.5*pi*22/NT)) ));
        plot(frontx,fronty,'Linewidth',1);
        frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*23/NT)))*(1+(0.75+0.7*sin(0.5*pi*23/NT)))*10+(0.75+0.7*sin(0.5*pi*23/NT)) ));
        plot(frontx,fronty,'Linewidth',1);
        frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*24/NT)))*(1+(0.75+0.7*sin(0.5*pi*24/NT)))*10+(0.75+0.7*sin(0.5*pi*24/NT)) ));
        plot(frontx,fronty,'Linewidth',1);
        frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*25/NT)))*(1+(0.75+0.7*sin(0.5*pi*25/NT)))*10+(0.75+0.7*sin(0.5*pi*25/NT)) ));
        plot(frontx,fronty,'Linewidth',1);
        frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*31/NT)))*(1+(0.75+0.7*sin(0.5*pi*31/NT)))*10+(0.75+0.7*sin(0.5*pi*31/NT)) ));
        plot(frontx,fronty,'Linewidth',1);
        frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*32/NT)))*(1+(0.75+0.7*sin(0.5*pi*32/NT)))*10+(0.75+0.7*sin(0.5*pi*32/NT)) ));
        plot(frontx,fronty,'Linewidth',1);
        frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*33/NT)))*(1+(0.75+0.7*sin(0.5*pi*33/NT)))*10+(0.75+0.7*sin(0.5*pi*33/NT)) ));
        plot(frontx,fronty,'Linewidth',1);
        frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*34/NT)))*(1+(0.75+0.7*sin(0.5*pi*34/NT)))*10+(0.75+0.7*sin(0.5*pi*34/NT)) ));
        plot(frontx,fronty,'Linewidth',1);
        frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*35/NT)))*(1+(0.75+0.7*sin(0.5*pi*35/NT)))*10+(0.75+0.7*sin(0.5*pi*35/NT)) ));
        plot(frontx,fronty,'Linewidth',1);
        set(gca,'YLim',[0 0.35]);%Y轴的数据显示范围
       
        
    end
    if problem == "F5"|| problem == "F6"||problem == "F7"||problem == "F9"
                m=0:0.01:1;frontx=m.^1.25 ;fronty=(1-m).^1.25;
        hold on;
        plot(frontx,fronty,'Linewidth',2);
                m=0:0.01:1;frontx=m.^(1.25+0.75*sin(pi/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi/NT));
        plot(frontx,fronty,'Linewidth',2);
                m=0:0.01:1;frontx=m.^(1.25+0.75*sin(pi*2/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*2/NT));
        plot(frontx,fronty,'Linewidth',2);
                m=0:0.01:1;frontx=m.^(1.25+0.75*sin(pi*3/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*3/NT));
        plot(frontx,fronty,'Linewidth',2);
                m=0:0.01:1;frontx=m.^(1.25+0.75*sin(pi*4/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*4/NT));
        plot(frontx,fronty,'Linewidth',2);
                m=0:0.01:1;frontx=m.^(1.25+0.75*sin(pi*5/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*5/NT));
        plot(frontx,fronty,'Linewidth',2);
                m=0:0.01:1;frontx=m.^(1.25+0.75*sin(pi*11/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*11/NT));
        plot(frontx,fronty,'Linewidth',2);
                m=0:0.01:1;frontx=m.^(1.25+0.75*sin(pi*12/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*12/NT));
        plot(frontx,fronty,'Linewidth',2);
                m=0:0.01:1;frontx=m.^(1.25+0.75*sin(pi*13/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*13/NT));
        plot(frontx,fronty,'Linewidth',2);
                m=0:0.01:1;frontx=m.^(1.25+0.75*sin(pi*14/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*14/NT));
        plot(frontx,fronty,'Linewidth',2);
                m=0:0.01:1;frontx=m.^(1.25+0.75*sin(pi*15/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*15/NT));
        plot(frontx,fronty,'Linewidth',2);
        set(gca,'YLim',[0 1]);%Y轴的数据显示范围
         
    end

    if problem == "dMOP2"|| problem == "dMOP1"
        frontx=0:0.01:1 ;fronty=1-frontx.^1.25;
        hold on
        plot(frontx,fronty,'Linewidth',2,'color',[0.5 0.5 0.5])
        frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*1/NT));
        plot(frontx,fronty,'Linewidth',2,'color',[0.5 0.5 0.5]) 
        frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*2/NT));
        plot(frontx,fronty,'Linewidth',2,'color',[0.5 0.5 0.5])
         frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*3/NT));
        plot(frontx,fronty,'Linewidth',2,'color',[0.5 0.5 0.5])
        frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*4/NT));
        plot(frontx,fronty,'Linewidth',2,'color',[0.5 0.5 0.5])      
        frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*5/NT));
        plot(frontx,fronty,'Linewidth',2,'color',[0.5 0.5 0.5])      
        frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*6/NT));
        plot(frontx,fronty,'Linewidth',2,'color',[0.5 0.5 0.5])       
        frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*7/NT));
        plot(frontx,fronty,'Linewidth',2,'color',[0.5 0.5 0.5])     
        frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*8/NT));
        plot(frontx,fronty,'Linewidth',2,'color',[0.5 0.5 0.5])
        frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*9/NT));
        plot(frontx,fronty,'Linewidth',2,'color',[0.5 0.5 0.5])
        frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*10/NT));
        plot(frontx,fronty,'Linewidth',2,'color',[0.5 0.5 0.5])
        
        frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*21/NT));
        plot(frontx,fronty,'Linewidth',2,'color',[0.5 0.5 0.5])
        frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*22/NT));
        plot(frontx,fronty,'Linewidth',2,'color',[0.5 0.5 0.5])
        frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*23/NT));
        plot(frontx,fronty,'Linewidth',2,'color',[0.5 0.5 0.5])
        frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*24/NT));
        plot(frontx,fronty,'Linewidth',2,'color',[0.5 0.5 0.5])
        frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*25/NT));
        plot(frontx,fronty,'Linewidth',2,'color',[0.5 0.5 0.5])
        frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*26/NT));
        plot(frontx,fronty,'Linewidth',2,'color',[0.5 0.5 0.5])
        frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*27/NT));
        plot(frontx,fronty,'Linewidth',2,'color',[0.5 0.5 0.5])
        frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*28/NT));
        plot(frontx,fronty,'Linewidth',2,'color',[0.5 0.5 0.5])
        frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*29/NT));
        plot(frontx,fronty,'Linewidth',2,'color',[0.5 0.5 0.5])
        frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*30/NT));
        plot(frontx,fronty,'Linewidth',2,'color',[0.5 0.5 0.5])
       set(gca,'YLim',[0 1]);%Y轴的数据显示范围

        
    end
    
    if problem == "FDA3"
        m=0:0.01:1;frontx=m.^(10^(2*sin(0.5*pi*0/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*0/NT)))))*(1+abs(sin(0.5*pi*0/NT)));
        hold on
        plot(frontx,fronty,'b.')
        m=0:0.01:1;frontx=m.^(10^(2*sin(0.5*pi*1/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*1/NT)))))*(1+abs(sin(0.5*pi*1/NT)));
        plot(frontx,fronty,'b.')  
        m=0:0.01:1;frontx=m.^(10^(2*sin(0.5*pi*2/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*2/NT)))))*(1+abs(sin(0.5*pi*2/NT)));
        plot(frontx,fronty,'b.')
        m=0:0.01:1;frontx=m.^(10^(2*sin(0.5*pi*3/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*3/NT)))))*(1+abs(sin(0.5*pi*3/NT)));
        plot(frontx,fronty,'b.')
        m=0:0.01:1;frontx=m.^(10^(2*sin(0.5*pi*4/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*4/NT)))))*(1+abs(sin(0.5*pi*4/NT)));
        plot(frontx,fronty,'b.')
        m=0:0.01:1;frontx=m.^(10^(2*sin(0.5*pi*5/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*5/NT)))))*(1+abs(sin(0.5*pi*5/NT)));
        plot(frontx,fronty,'b.')
        m=0:0.01:1;frontx=m.^(10^(2*sin(0.5*pi*11/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*11/NT)))))*(1+abs(sin(0.5*pi*11/NT)));
        plot(frontx,fronty,'b.')
        m=0:0.01:1;frontx=m.^(10^(2*sin(0.5*pi*12/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*12/NT)))))*(1+abs(sin(0.5*pi*12/NT)));
        plot(frontx,fronty,'b.')
        m=0:0.01:1;frontx=m.^(10^(2*sin(0.5*pi*13/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*13/NT)))))*(1+abs(sin(0.5*pi*13/NT)));
        plot(frontx,fronty,'b.')
        m=0:0.01:1;frontx=m.^(10^(2*sin(0.5*pi*14/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*14/NT)))))*(1+abs(sin(0.5*pi*14/NT)));
        plot(frontx,fronty,'b.')
        m=0:0.01:1;frontx=m.^(10^(2*sin(0.5*pi*15/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*15/NT)))))*(1+abs(sin(0.5*pi*15/NT)));
        plot(frontx,fronty,'b.')
        m=0:0.01:1;frontx=m.^(10^(2*sin(0.5*pi*21/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*21/NT)))))*(1+abs(sin(0.5*pi*21/NT)));
        plot(frontx,fronty,'b.')
        m=0:0.01:1;frontx=m.^(10^(2*sin(0.5*pi*22/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*22/NT)))))*(1+abs(sin(0.5*pi*22/NT)));
        plot(frontx,fronty,'b.')
        m=0:0.01:1;frontx=m.^(10^(2*sin(0.5*pi*23/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*23/NT)))))*(1+abs(sin(0.5*pi*23/NT)));
        plot(frontx,fronty,'b.')
        m=0:0.01:1;frontx=m.^(10^(2*sin(0.5*pi*24/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*24/NT)))))*(1+abs(sin(0.5*pi*24/NT)));
        plot(frontx,fronty,'b.')
        m=0:0.01:1;frontx=m.^(10^(2*sin(0.5*pi*25/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*25/NT)))))*(1+abs(sin(0.5*pi*25/NT)));
        plot(frontx,fronty,'b.')
        m=0:0.01:1;frontx=m.^(10^(2*sin(0.5*pi*31/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*31/NT)))))*(1+abs(sin(0.5*pi*31/NT)));
        plot(frontx,fronty,'b.')
        m=0:0.01:1;frontx=m.^(10^(2*sin(0.5*pi*32/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*32/NT)))))*(1+abs(sin(0.5*pi*32/NT)));
        plot(frontx,fronty,'b.')
        m=0:0.01:1;frontx=m.^(10^(2*sin(0.5*pi*33/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*33/NT)))))*(1+abs(sin(0.5*pi*33/NT)));
        plot(frontx,fronty,'b.')
        m=0:0.01:1;frontx=m.^(10^(2*sin(0.5*pi*34/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*34/NT)))))*(1+abs(sin(0.5*pi*34/NT)));
        plot(frontx,fronty,'b.')
        m=0:0.01:1;frontx=m.^(10^(2*sin(0.5*pi*35/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*35/NT)))))*(1+abs(sin(0.5*pi*35/NT)));
        plot(frontx,fronty,'b.')
        set(gca,'XLim',[0 1]);
    end
    
    if problem == "FDA5"
        P = UniformPoint(2000,3);
        P = P./repmat((1/(1+(sin(0.5*pi*1/NT))^4))*sqrt(sum(P.^2,2)),1,3);
        Draw(P,'ok','Markeredgecolor',[.0 .9 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',1);    
        
        P = UniformPoint(2000,3);
        P = P./repmat((1/(1+(sin(0.5*pi*2/NT))^4))*sqrt(sum(P.^2,2)),1,3);
        Draw(P,'ok','Markeredgecolor',[.0 .9 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',1); 
        
        P = UniformPoint(2000,3);
        P = P./repmat((1/(1+(sin(0.5*pi*3/NT))^4))*sqrt(sum(P.^2,2)),1,3);
        Draw(P,'ok','Markeredgecolor',[.0 .9 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',1); 
        
        P = UniformPoint(2000,3);
        P = P./repmat((1/(1+(sin(0.5*pi*4/NT))^4))*sqrt(sum(P.^2,2)),1,3);
        Draw(P,'ok','Markeredgecolor',[.0 .9 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',1); 
        
        P = UniformPoint(2000,3);
        P = P./repmat((1/(1+(sin(0.5*pi*5/NT))^4))*sqrt(sum(P.^2,2)),1,3);
        Draw(P,'ok','Markeredgecolor',[.0 .9 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',1);
    end

end