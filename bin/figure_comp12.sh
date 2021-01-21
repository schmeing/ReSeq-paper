#!/bin/bash
cd $1
# Get maximum size from all plots
maxsize=$(cat <(pdfinfo figure1.pdf) <(pdfinfo figure2.pdf) <(pdfinfo figure3.pdf) <(pdfinfo figure4.pdf) <(pdfinfo figure5.pdf) <(pdfinfo figure6.pdf) <(pdfinfo figure7.pdf) <(pdfinfo figure8.pdf) <(pdfinfo figure9.pdf) <(pdfinfo figure10.pdf) <(pdfinfo figure11.pdf) <(pdfinfo figure12.pdf) | grep 'Page size:' | awk '{x=x<$3?$3:x; y=y<$5?$5:y}END{print x,y}')

# Extend all plots to maximum size (with whitespace, no scaling: All plots are expected to be at the same scalling and the size difference caused by trimming)
for fi in figure[1-9].pdf figure1[0-2].pdf; do
    cropvals=$(cat <(pdfinfo $fi) | grep 'Page size:' | awk -v maxval="$maxsize" '{split(maxval,mv," ");printf "%i %i %i %i", (mv[1]-$3)/2, (mv[2]-$5)/2, (mv[1]-$3+1)/2, (mv[2]-$5+1)/2}');  pdfcrop --margins "$cropvals" $fi c$fi 1>/dev/null
done

# Merge plots to have final width of Landscape A4: 841.89 x 595.276 pts
pdfjam cfigure[1-9].pdf cfigure1[0-2].pdf --delta '10px 20px' --papersize '{842px,5000px}' --nup 3x4 --outfile final.pdf 2>/dev/null
pdfcrop --margins '0 20 0 0' final.pdf cfinal.pdf 1>/dev/null

# Add letters
echo "a" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > a.pdf 
pdfcrop a.pdf ca.pdf 1>/dev/null #a(11 x 12 pts)
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo ca.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%i %i %i %i", 0, 4, $3-l[3], $5-l[5]-4 }'); pdfcrop --margins "$cropvals" a.pdf ca.pdf 1>/dev/null
pdftk ca.pdf background cfinal.pdf output scfinal.pdf

echo "b" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > b.pdf
pdfcrop b.pdf cb.pdf 1>/dev/null #b(11 x 16 pts)
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo cb.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%i %i %i %i", ($3+10)/3, 0, $3-($3+10)/3-l[3], $5-l[5] }'); pdfcrop --margins "$cropvals" b.pdf cb.pdf 1>/dev/null
pdftk cb.pdf background scfinal.pdf output scfinal2.pdf && mv -f scfinal2.pdf scfinal.pdf

echo "c" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > c.pdf
pdfcrop c.pdf cc.pdf 1>/dev/null #c(11 x 12 pts)
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo cc.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%i %i %i %i", ($3+10)*2/3, 4, $3-($3+10)*2/3-l[3], $5-l[5]-4 }'); pdfcrop --margins "$cropvals" c.pdf cc.pdf 1>/dev/null
pdftk cc.pdf background scfinal.pdf output scfinal2.pdf && mv -f scfinal2.pdf scfinal.pdf

echo "d" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > d.pdf
pdfcrop d.pdf cd.pdf 1>/dev/null #d(11 x 16 pts)
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo cd.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%i %i %i %i", 0, $5/4, $3-l[3], $5-$5/4-l[5] }'); pdfcrop --margins "$cropvals" d.pdf cd.pdf 1>/dev/null
pdftk cd.pdf background scfinal.pdf output scfinal2.pdf && mv -f scfinal2.pdf scfinal.pdf

echo "e" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > e.pdf
pdfcrop e.pdf ce.pdf 1>/dev/null #e(11 x 12 pts)
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo ce.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%i %i %i %i", ($3+10)/3, $5/4+4, $3-($3+10)/3-l[3], $5-$5/4-l[5]-4 }'); pdfcrop --margins "$cropvals" e.pdf ce.pdf 1>/dev/null
pdftk ce.pdf background scfinal.pdf output scfinal2.pdf && mv -f scfinal2.pdf scfinal.pdf

echo "f" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > f.pdf
pdfcrop f.pdf cf.pdf 1>/dev/null #f(7 x 16 pts)
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo cf.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%i %i %i %i", ($3+10)*2/3, $5/4, $3-($3+10)*2/3-l[3], $5-$5/4-l[5] }'); pdfcrop --margins "$cropvals" f.pdf cf.pdf 1>/dev/null
pdftk cf.pdf background scfinal.pdf output scfinal2.pdf && mv -f scfinal2.pdf scfinal.pdf

echo "g" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > g.pdf
pdfcrop g.pdf cg.pdf 1>/dev/null #g(11 x 16 pts)
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo cg.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%i %i %i %i", 0, $5*2/4+4, $3-l[3], $5-$5*2/4-l[5]-4 }'); pdfcrop --margins "$cropvals" g.pdf cg.pdf 1>/dev/null
pdftk cg.pdf background scfinal.pdf output scfinal2.pdf && mv -f scfinal2.pdf scfinal.pdf

echo "h" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > h.pdf
pdfcrop h.pdf ch.pdf 1>/dev/null #h(10 x 16 pts)
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo ch.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%i %i %i %i", ($3+10)/3, $5*2/4, $3-($3+10)/3-l[3], $5-$5*2/4-l[5] }'); pdfcrop --margins "$cropvals" h.pdf ch.pdf 1>/dev/null
pdftk ch.pdf background scfinal.pdf output scfinal2.pdf && mv -f scfinal2.pdf scfinal.pdf

echo "i" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > i.pdf
pdfcrop i.pdf ci.pdf 1>/dev/null #i(4 x 16 pts)
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo ci.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%i %i %i %i", ($3+10)*2/3, $5*2/4, $3-($3+10)*2/3-l[3], $5-$5*2/4-l[5] }'); pdfcrop --margins "$cropvals" i.pdf ci.pdf 1>/dev/null
pdftk ci.pdf background scfinal.pdf output scfinal2.pdf && mv -f scfinal2.pdf scfinal.pdf

echo "j" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > j.pdf
pdfcrop j.pdf cj.pdf 1>/dev/null #j(6 x 20 pts)
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo cj.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%i %i %i %i", 0, $5*3/4, $3-l[3], $5-$5*3/4-l[5] }'); pdfcrop --margins "$cropvals" j.pdf cj.pdf 1>/dev/null
pdftk cj.pdf background scfinal.pdf output scfinal2.pdf && mv -f scfinal2.pdf scfinal.pdf

echo "k" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > k.pdf
pdfcrop k.pdf ck.pdf 1>/dev/null #k(10 x 16 pts)
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo ck.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%i %i %i %i", ($3+10)/3, $5*3/4, $3-($3+10)/3-l[3], $5-$5*3/4-l[5] }'); pdfcrop --margins "$cropvals" k.pdf ck.pdf 1>/dev/null
pdftk ck.pdf background scfinal.pdf output scfinal2.pdf && mv -f scfinal2.pdf scfinal.pdf

echo "l" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > l.pdf
pdfcrop l.pdf cl.pdf 1>/dev/null #l(4 x 16 pts)
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo cl.pdf) | grep 'Page size:' | awk '{getline line; split(line,l," "); printf "%i %i %i %i", ($3+10)*2/3, $5*3/4, $3-($3+10)*2/3-l[3], $5-$5*3/4-l[5] }'); pdfcrop --margins "$cropvals" l.pdf cl.pdf 1>/dev/null
pdftk cl.pdf background scfinal.pdf output scfinal2.pdf && mv -f scfinal2.pdf scfinal.pdf










