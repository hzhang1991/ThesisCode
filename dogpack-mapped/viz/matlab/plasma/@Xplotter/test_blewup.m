function blewup = test_blewup(xplotter, outputdir)
  [blewup,garbage] = system(['tail ' outputdir '/[^q]*.dat | grep -i nan']);
  blewup = ~blewup;
end
