#!/bin/bash
cd $1
# Bring screenshots to height of plots through scaling and adding white space (Considering that they should be twice as long)
img2pdf SRR490124-screen.png | pdfjam /dev/stdin --papersize "{$(pdfinfo SRR490124-errors.pdf | grep 'Page size:' | awk '{printf "%ipx,%ipx", $3*1.5, $5}')}" --outfile SRR490124-screen.pdf 2>/dev/null
pdfcrop SRR490124-screen.pdf cSRR490124-screen.pdf 1>/dev/null
cropvals=$(cat <(pdfinfo cSRR490124-screen.pdf) <(pdfinfo SRR490124-errors.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%i %i %i %i", 0, (l[5]-$5)/2, 0, (l[5]-$5+1)/2}');  pdfcrop --margins "$cropvals" cSRR490124-screen.pdf c2SRR490124-screen.pdf 1>/dev/null && mv -f c2SRR490124-screen.pdf cSRR490124-screen.pdf
img2pdf S5L001-screen.png | pdfjam /dev/stdin --papersize "{$(pdfinfo S5L001-errors.pdf | grep 'Page size:' | awk '{printf "%ipx,%ipx", $3*1.5, $5}')}" --outfile S5L001-screen.pdf 2>/dev/null
pdfcrop S5L001-screen.pdf cS5L001-screen.pdf 1>/dev/null
cropvals=$(cat <(pdfinfo cS5L001-screen.pdf) <(pdfinfo cSRR490124-screen.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%i %i %i %i", (l[3]-$3+1)/2, (l[5]-$5)/2, (l[3]-$3)/2, (l[5]-$5+1)/2}');  pdfcrop --margins "$cropvals" cS5L001-screen.pdf c2S5L001-screen.pdf 1>/dev/null && mv -f c2S5L001-screen.pdf cS5L001-screen.pdf

# Add whitespace to right of plots to prevent the white space from being equally distributed between left and right
cropvals=$(cat <(pdfinfo SRR490124-errors.pdf) <(pdfinfo SRR490124-screen.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%i %i %i %i", 0, 0, l[3]-$3, 0}'); pdfcrop --margins "$cropvals" SRR490124-errors.pdf cSRR490124-errors.pdf 1>/dev/null
cropvals=$(cat <(pdfinfo S5L001-errors.pdf) <(pdfinfo S5L001-screen.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%i %i %i %i", 0, 0, l[3]-$3, 0}'); pdfcrop --margins "$cropvals" S5L001-errors.pdf cS5L001-errors.pdf 1>/dev/null

# Merge columns
pdfjam cSRR490124-screen.pdf cSRR490124-errors.pdf --delta '10px 0' --nup 2x2 --landscape --outfile top.pdf 2>/dev/null
pdfjam cS5L001-screen.pdf cS5L001-errors.pdf --delta '10px 0' --nup 2x2 -landscape --outfile middle.pdf 2>/dev/null

# Crop new rows and make them the same width (Landscape A4: 841.89 x 595.276 pts)
pdfcrop top.pdf ctop.pdf 1>/dev/null
pdfjam ctop.pdf --landscape --outfile jtop.pdf 2>/dev/null
pdfcrop jtop.pdf ctop.pdf 1>/dev/null
pdfcrop middle.pdf cmiddle.pdf 1>/dev/null
pdfjam cmiddle.pdf --papersize "{$(pdfinfo ctop.pdf | grep 'Page size:' | awk '{printf "%ipx,%ipx", $3, $5}')}" --outfile jmiddle.pdf 2>/dev/null
pdfcrop jmiddle.pdf cmiddle.pdf 1>/dev/null
cropvals=$(cat <(pdfinfo cmiddle.pdf) <(pdfinfo ctop.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%i %i %i %i", l[3]-$3, 0, 0, 0}');  pdfcrop --margins "$cropvals" cmiddle.pdf cmiddle2.pdf 1>/dev/null && mv -f cmiddle2.pdf cmiddle.pdf

# Prepare bottom row
img2pdf rep-screen.png | pdfjam /dev/stdin --papersize "{$(pdfinfo ctop.pdf | grep 'Page size:' | awk '{printf "%ipx,%ipx", $3, $5}')}" --outfile rep-screen.pdf 2>/dev/null
pdfcrop rep-screen.pdf crep-screen.pdf 1>/dev/null
cropvals=$(cat <(pdfinfo crep-screen.pdf) <(pdfinfo ctop.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%i %i %i %i", 0, 0, 0, l[5]-$5}');  pdfcrop --margins "$cropvals" rep-screen.pdf crep-screen.pdf 1>/dev/null

# Merge rows  
pdfjam ctop.pdf cmiddle.pdf crep-screen.pdf --delta '0 20px' --nup 1x3 --papersize "{$(pdfinfo ctop.pdf | grep 'Page size:' | awk '{printf "%ipx,%ipx", $3, $5*3+2*20}')}" --outfile final.pdf 2>/dev/null
pdfcrop --margins '0 20 0 0' final.pdf cfinal.pdf 1>/dev/null

# Add letters
echo "a" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > a.pdf
pdfcrop a.pdf ca.pdf 1>/dev/null
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo ca.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%i %i %i %i", 0, 4, $3-l[3], $5-l[5]-4 }'); pdfcrop --margins "$cropvals" a.pdf ca.pdf 1>/dev/null
pdftk ca.pdf background cfinal.pdf output scfinal.pdf

echo "c" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > c.pdf
pdfcrop c.pdf cc.pdf 1>/dev/null
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo cc.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%i %i %i %i", $3*3/5, 4, $3-$3*3/5-l[3], $5-l[5]-4 }'); pdfcrop --margins "$cropvals" c.pdf cc.pdf 1>/dev/null
pdftk cc.pdf background scfinal.pdf output scfinal2.pdf && mv -f scfinal2.pdf scfinal.pdf

echo "b" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > b.pdf
pdfcrop b.pdf cb.pdf 1>/dev/null
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo cb.pdf) <(pdfinfo ctop.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); getline line; split(line,l2," "); printf "%i %i %i %i", 0, 20+l2[5], $3-l[3], $5-20-l2[5]-l[5] }'); pdfcrop --margins "$cropvals" b.pdf cb.pdf 1>/dev/null
pdftk cb.pdf background scfinal.pdf output scfinal2.pdf && mv -f scfinal2.pdf scfinal.pdf

echo "d" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > d.pdf
pdfcrop d.pdf cd.pdf 1>/dev/null
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo cd.pdf) <(pdfinfo ctop.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); getline line; split(line,l2," "); printf "%i %i %i %i", $3*3/5, 20+l2[5], $3-$3*3/5-l[3], $5-20-l2[5]-l[5] }'); pdfcrop --margins "$cropvals" d.pdf cd.pdf 1>/dev/null
pdftk cd.pdf background scfinal.pdf output scfinal2.pdf && mv -f scfinal2.pdf scfinal.pdf

echo "e" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > e.pdf
pdfcrop e.pdf ce.pdf 1>/dev/null
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo ce.pdf) <(pdfinfo ctop.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); getline line; split(line,l2," "); printf "%i %i %i %i", 0, 2*(20+l2[5])+4, $3-l[3], $5-2*(20+l2[5])-l[5]-4 }'); pdfcrop --margins "$cropvals" e.pdf ce.pdf 1>/dev/null
pdftk ce.pdf background scfinal.pdf output scfinal2.pdf && mv -f scfinal2.pdf scfinal.pdf
