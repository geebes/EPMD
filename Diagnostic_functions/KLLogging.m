function stop = KLLogging(optimValues,state,species)
persistent kllog iters stopnow
switch state
    case 'init'
        stopnow = false;
        kllog = [];
        iters = [];
%         h = figure;
%         c = uicontrol('Style','pushbutton','String','Stop','Position', ...
%             [10 10 50 20],'Callback',@stopme);
    case 'iter'
        kllog = [kllog; optimValues.fval,log(norm(optimValues.grad))];
        assignin('base','history',kllog)
        iters = [iters; optimValues.iteration];
        if length(iters) > 1
%             figure(h)
%             subplot(2,1,2)
%             plot(iters,kllog);
%             xlabel('Iterations')
%             ylabel('Loss and Gradient')
%             legend('Divergence','log(norm(gradient))')
%             subplot(2,1,1)
%             clr = optimValues.Y;
%             clr = normalize(clr,'range');
%             scatter(optimValues.Y(:,1),optimValues.Y(:,2),50,clr,'filled')
            cmap=turbo(8);
            cla
            scatter(optimValues.Y(:,1),optimValues.Y(:,2),75,cmap(species,:),'filled')
            hold on
            scatter(optimValues.Y(:,1),optimValues.Y(:,2),75,'k')
            
            axis tight;
            ax=axis;
            axx=abs(optimValues.Y(:,1:2));
            axx=prctile(axx(:),95).*1.5;
            axlim=max([axx,500]);
            axis([-axlim axlim -axlim axlim])
%             title('Embedding')
            legend off
            title('t-distributed Stochastic Neighbour Embedding')
            set(gca,'XTick',[],'YTick',[])
            box  on
            drawnow
        end
    case 'done'
        % Nothing here
end
stop = stopnow;

function stopme(~,~)
stopnow = true;
end
end