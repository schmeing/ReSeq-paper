#!/bin/bash
cd $1
# Get maximum size from all plots
maxsize=$(cat <(pdfinfo figure1.pdf) <(pdfinfo figure2.pdf) | grep 'Page size:' | awk '{x=x<$3?$3:x; y=y<$5?$5:y}END{print x,y}')
xf1ratio=$(cat <(pdfinfo figure1.pdf) <(pdfinfo figure2.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%f", $3/($3+l[3])}')

# Extend all plots to maximum size (with whitespace, no scaling: All plots are expected to be at the same scalling and the size difference caused by trimming)
cropvals=$(cat <(pdfinfo figure1.pdf) | grep 'Page size:' | awk -v maxval="$maxsize" '{split(maxval,mv," ");printf "%i %i %i %i", (mv[1]-$3), (mv[2]-$5)/2, 0, (mv[2]-$5+1)/2}');  pdfcrop --margins "$cropvals" figure1.pdf cfigure1.pdf 1>/dev/null
cropvals=$(cat <(pdfinfo figure2.pdf) | grep 'Page size:' | awk -v maxval="$maxsize" '{split(maxval,mv," ");printf "%i %i %i %i", 0, (mv[2]-$5)/2, (mv[1]-$3), (mv[2]-$5+1)/2}');  pdfcrop --margins "$cropvals" figure2.pdf cfigure2.pdf 1>/dev/null

# Merge plots to have final width of Landscape A4: 841.89 x 595.276 pts
pdfjam cfigure[1-2].pdf --delta '10px 0px' --landscape --nup 2x1 --outfile final.pdf 2>/dev/null
pdfcrop --margins '0 0 0 1' final.pdf cfinal.pdf 1>/dev/null
pdfjam cfinal.pdf --landscape --outfile jfinal.pdf 2>/dev/null
pdfcrop --margins '0 20 0 1' jfinal.pdf cfinal.pdf 1>/dev/null

# Add letters
echo "a" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > a.pdf 
pdfcrop a.pdf ca.pdf 1>/dev/null #a(11 x 12 pts)
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo ca.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%i %i %i %i", 0, 4, $3-l[3], $5-l[5]-4 }'); pdfcrop --margins "$cropvals" a.pdf ca.pdf 1>/dev/null
pdftk ca.pdf background cfinal.pdf output scfinal.pdf

echo "b" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > b.pdf
pdfcrop b.pdf cb.pdf 1>/dev/null #b(11 x 16 pts)
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo cb.pdf) | grep 'Page size:' | awk -v xf1ratio="$xf1ratio" '{getline line; split(line,l," "); printf "%i %i %i %i", ($3-10)*xf1ratio+10, 0, $3-(($3-10)*xf1ratio+10)-l[3], $5-l[5] }'); pdfcrop --margins "$cropvals" b.pdf cb.pdf 1>/dev/null
pdftk cb.pdf background scfinal.pdf output scfinal2.pdf && mv -f scfinal2.pdf scfinal.pdf
