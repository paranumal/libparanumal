
files = ls;

cnt = 1;
for n=1:size(files,1)

  if(strfind(files(n,:), '_occa_'))
    name= files(n,:);
    m = length(name);
    while(m>1)
      if(isspace(name(m)))
	name = name(1:end-1);
      else
	break;
      end
      m=m-1;
    end

    figure(cnt);
    rooflinePlot(name)
    cnt = cnt+1;
  end
end

