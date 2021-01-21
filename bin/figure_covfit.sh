#!/bin/bash
cd $1
pdfjam gcbias.pdf surbias.pdf --delta '10px 0' --nup 2x1 --landscape --outfile top.pdf 2>/dev/null
pdfjam samplemean.pdf dispa.pdf dispb.pdf --delta '10px 0' --nup 3x1 --landscape --outfile bottom.pdf 2>/dev/null
pdfcrop top.pdf ctop.pdf 1>/dev/null
pdfcrop bottom.pdf cbottom.pdf 1>/dev/null
pdfjam ctop.pdf cbottom.pdf --nup 1x2 --landscape --delta "0 -$(cat <(pdfinfo ctop.pdf) <(pdfinfo cbottom.pdf) | grep 'Page size:' | awk 'BEGIN{y=596}{y-=$5}END{printf "%i", (y+1)/2-20}')px" --outfile final.pdf 2>/dev/null
pdfcrop --margins '0 20 0 0' final.pdf cfinal.pdf 1>/dev/null

echo "a" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > a.pdf
pdfcrop a.pdf ca.pdf 1>/dev/null
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo ca.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%i %i %i %i", 0, 4, $3-l[3], $5-l[5]-4 }'); pdfcrop --margins "$cropvals" a.pdf ca.pdf 1>/dev/null
pdftk ca.pdf background cfinal.pdf output scfinal.pdf

echo "b" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > b.pdf
pdfcrop b.pdf cb.pdf 1>/dev/null
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo cb.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%i %i %i %i", $3/2, 0, $3-$3/2-l[3], $5-l[5] }'); pdfcrop --margins "$cropvals" b.pdf cb.pdf 1>/dev/null
pdftk cb.pdf background scfinal.pdf output scfinal2.pdf && mv -f scfinal2.pdf scfinal.pdf

echo "c" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > c.pdf
pdfcrop c.pdf cc.pdf 1>/dev/null
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo cc.pdf) <(pdfinfo ctop.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); getline line; split(line,l2," "); printf "%i %i %i %i", 0, 20+l2[5]+4, $3-l[3], $5-20-l2[5]-l[5]-4 }'); pdfcrop --margins "$cropvals" c.pdf cc.pdf 1>/dev/null
pdftk cc.pdf background scfinal.pdf output scfinal2.pdf && mv -f scfinal2.pdf scfinal.pdf

echo "d" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > d.pdf
pdfcrop d.pdf cd.pdf 1>/dev/null
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo cd.pdf) <(pdfinfo ctop.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); getline line; split(line,l2," "); printf "%i %i %i %i", $3/3, 20+l2[5], $3-$3/3-l[3], $5-20-l2[5]-l[5] }'); pdfcrop --margins "$cropvals" d.pdf cd.pdf 1>/dev/null
pdftk cd.pdf background scfinal.pdf output scfinal2.pdf && mv -f scfinal2.pdf scfinal.pdf

echo "e" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > e.pdf
pdfcrop e.pdf ce.pdf 1>/dev/null
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo ce.pdf) <(pdfinfo ctop.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); getline line; split(line,l2," "); printf "%i %i %i %i", $3*2/3, 20+l2[5]+4, $3-$3*2/3-l[3], $5-20-l2[5]-l[5]-4 }'); pdfcrop --margins "$cropvals" e.pdf ce.pdf 1>/dev/null
pdftk ce.pdf background scfinal.pdf output scfinal2.pdf && mv -f scfinal2.pdf scfinal.pdf
