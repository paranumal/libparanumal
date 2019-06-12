function myprint(fname)

  fig = gcf;
fig.PaperPositionMode = 'auto'
  fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

print(fig,fname,'-dpdf', '-painters')
