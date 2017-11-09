html:index.html
	@echo index.html

index.html:../iteratorCodeLib.md html.tem
	pandoc --mathjax --template html.tem -s -N -o "$@" "$<"

noTocHtml.html:../iteratorCodeLib.md noTocHtml.latex
	pandoc --mathjax --template noTocHtml.latex -s -N -o "$@" "$<"

kkke.me:index.html
	scp "$<" acm@kkke.me:html/index.html

kkke.me/pdf:iteratorCodeLibrary.pdf
	scp "$<" acm@kkke.me:html/ICL.pdf

pdf:iteratorCodeLibrary.pdf
	@echo iteratorCodeLibrary.pdf

iteratorCodeLibrary.pdf:../iteratorCodeLib.md pdf.tex
	pandoc --latex-engine=xelatex --template pdf.tex -s --toc -N -o "$@" "$<"

clearHtml:
	@rm *.html 2>/dev/null && echo clear html ||echo no html file

clearPdf:
	@rm *.pdf 2>/dev/null && echo clear pdf ||echo no pdf file

clear:clearPdf clearHtml
	@echo clear


