
clf

Haxes = axes('position',[.1 .1 .8 .3]);
set(Haxes,'fontsize',20)

for Frame=[0 1 2 4]
  for mq=1:2
    plotframe1
    if mq==1
       axis([-5 5 0.0 2.0])
       end
    if mq==2
       axis([-5 5 -1.0 1.0])
       end

    eval(['print tworaref_t' num2str(Frame) 'p' num2str(mq) ' -deps'])
    query
    end
  end

